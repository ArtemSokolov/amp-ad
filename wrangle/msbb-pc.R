## Wrangling of MSBB RNAseq data and matching clinical covariates
##   Reduces gene space to the protein coding regions
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
suppressMessages(library( synapseClient ))

## Parse local directory specification
argv <- commandArgs( trailingOnly = TRUE )
if( length(argv) == 0 )
{
    cat( "NOTE: No directory specified on command line. Using default.\n" )
    local.dir <- "/data/AMP-AD/MSBB"
} else { local.dir <- argv[1] }

## Create directory if it doesn't exist
dir.create( local.dir, showWarnings=FALSE )
cat( "Wrangling MSBB dataset to", local.dir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapseLogin( rememberMe=TRUE )

## Load biotype annotations and retrieve names of protein-coding genes
cat( "Downloading biotype annotations...\n" )
fnBT <- synGet( "syn14236139", downloadLocation = local.dir )@filePath
BT <- read_csv( fnBT, col_types=cols() ) %>% filter( gene_biotype=="protein_coding" ) %>%
    select( ENSEMBL = gene_id, HUGO = gene_name )

## Read raw expression matrix
cat( "Downloading expression data...\n" )
fnX <- synGet( "syn7809023", version=1, downloadLocation = local.dir )@filePath
Xraw <- read.delim( fnX, check.names=FALSE ) %>% rownames_to_column( "ENSEMBL" )

## Map gene IDs to HUGO
## PINX1 is the only gene with duplicate entries, but one of the entries has
##   a higher total count, so we keep it and discard the other entry.
X <- inner_join( BT, Xraw, by="ENSEMBL" ) %>% filter( ENSEMBL != "ENSG00000258724" ) %>%
    select( -ENSEMBL ) %>% gather( barcode, Value, -HUGO )

## Match sample barcodes against individual IDs and brain region information
cat( "Annotating samples with brain region...\n" )
fnZ <- synGet( "syn6100548", downloadLocation = local.dir )@filePath
XZ <- suppressMessages( read_csv(fnZ) ) %>%
    select( BrodmannArea, barcode, individualIdentifier ) %>% distinct %>%
    mutate( barcode = as.character(barcode) ) %>% inner_join( X, by = "barcode" )

## Load clinical covariates and combine with the expression matrix
cat( "Downloading clinical covariates...\n" )
fnY <- synGet( "syn6101474", downloadLocation = local.dir )@filePath
XY <- suppressMessages( read_csv( fnY ) ) %>% select( individualIdentifier, PMI, AOD, CDR, bbscore ) %>%
    inner_join( XZ, by = "individualIdentifier" )

## Flatten the matrix to samples-by-(clin+genes) and save to file
cat( "Finalizing...\n" )
RR <- spread( XY, HUGO, Value ) %>%
    rename( ID = individualIdentifier, Braak = bbscore, Barcode = barcode )
fnOut <- file.path( local.dir, "msbb-pc.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
write_tsv( RR, fnOut )
