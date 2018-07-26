## A variant of rosmap.R that reduces the space of molecular
##   features to protein-coding (pc) genes
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
suppressMessages(library( synapseClient ))
library( stringr )

## Parse local directory specification
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) == 0 )
{
    cat( "NOTE: No directory specified on command line. Using default.\n" )
    local.dir <- "/data/AMP-AD/ROSMAP"
} else { local.dir <- argv[1] }

## Create directory if it doesn't exist
dir.create( local.dir, showWarnings=FALSE )
cat( "Wrangling ROS/MAP dataset to", local.dir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapseLogin( rememberMe=TRUE )

## Read raw expression matrix
cat( "Downloading expression data...\n" )
fnX <- synGet( "syn3505720", downloadLocation = local.dir )@filePath
Xraw <- read_tsv( fnX, col_types = cols() )

## Load biotype annotations and retrieve names of protein-coding genes
cat( "Downloading biotype annotations...\n" )
fnBT <- synGet( "syn14236139", downloadLocation = local.dir )@filePath
BT <- read_csv( fnBT, col_types=cols() ) %>% filter( gene_biotype=="protein_coding" ) %>%
    select( ENSEMBL = gene_id, HUGO = gene_name )

## Map ENSEMBL Gene IDs to HUGO
## Reduce feature space to protein-coding genes
## There is a single duplicate: ENSG00000254093.3 and ENSG00000258724.1 map to the
##   same HUGO ID. However, ENSG00000258724.1 is almost all 0s, so we drop it.
cat( "Mapping gene IDs to HUGO...\n" )
X <- Xraw %>% filter( gene_id != "ENSG00000258724.1" ) %>%
    mutate( ENSEMBL = str_split( gene_id, "\\.", simplify=TRUE )[,1] ) %>%
    inner_join( BT, by="ENSEMBL" ) %>% select( -tracking_id, -gene_id, -ENSEMBL )

## Log-transform the data and combine the replicates
cat( "Additional processing...\n" )
fmed <- function(x) {x %>% as.matrix %>% apply( 1, median )}
XX <- X %>% mutate_at( vars(-HUGO), ~log2(.x+1) ) %>%
    mutate( `492_120515_j` = fmed(select( ., contains("492_120515") )) ) %>%
    select( -`492_120515_0`, -`492_120515_6`, -`492_120515_7` ) %>%
    gather( rnaseq_id, Value, -HUGO ) %>%
    mutate( rnaseq_id = str_sub( rnaseq_id, 0, -3 ) )

## Match sample IDs against individual identifiers
cat( "Matching sample and individual IDs...\n" )
fnZ <- synGet( "syn3382527", downloadLocation = local.dir )@filePath
XZ <- suppressMessages( read_csv(fnZ) ) %>% select( projid, rnaseq_id ) %>% na.omit %>%
    distinct %>% inner_join( XX, ., by="rnaseq_id" )

## Match expression data up against the following clinical covariates:
## ID, PMI, AOD, CDR, Braak, BrodmannArea
cat( "Matching against clinical covariates...\n" )
fnY <- synGet( "syn3191087", downloadLocation = local.dir )@filePath
Y <- suppressWarnings( suppressMessages( read_csv(fnY) ) ) %>%
    select( projid, PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc ) %>%
    mutate( BrodmannArea = "BM9,BM46" )

## Combining everything into a common data frame
cat( "Finalizing...\n" )
XY <- inner_join( Y, XZ, by="projid" ) %>% rename( ID = projid, Barcode = rnaseq_id ) %>%
    spread( HUGO, Value )

## Write out wrangled dataset to file
fnOut <- file.path( local.dir, "rosmap-pc.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
write_tsv( XY, fnOut )
