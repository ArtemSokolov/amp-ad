## Wrangling of Mayo clinic data and matching clinical annotations
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
    local.dir <- "/data/AMP-AD/Mayo"
} else { local.dir <- argv[1] }

## Create directory if it doesn't exist
dir.create( local.dir, showWarnings=FALSE )
cat( "Wrangling Mayo dataset to", local.dir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapseLogin( rememberMe=TRUE )

## Load biotype annotations and retrieve names of protein-coding genes
cat( "Downloading biotype annotations...\n" )
fnBT <- synGet( "syn14236139", downloadLocation = local.dir )@filePath
BT <- read_csv( fnBT, col_types=cols() ) %>% filter( gene_biotype=="protein_coding" ) %>%
    select( ensembl_id = gene_id, HUGO = gene_name )

## Read raw expression matrices (one for cerebellum, one for temporal cortex)
cat( "Downloading expression data...\n" )
ct <- cols( .default = col_double(), ensembl_id = col_character() )
XCBE <- synGet( "syn5201007", downloadLocation = local.dir )@filePath %>% read_tsv( col_types = ct )
XTCX <- synGet( "syn4650265", downloadLocation = local.dir )@filePath %>% read_tsv( col_types = ct )

## Combine the two expression matrices
## 742_CER has a low number of counts; drop it
## Observation made by the following:
##   XCBE %>% summarize_at( -1, sum ) %>% gather( Sample, Total ) %>% arrange( Total )
Xraw <- XCBE %>% select( -`742_CER` ) %>% inner_join( XTCX, by="ensembl_id" ) %>%
    inner_join( BT, ., by="ensembl_id" )

## The following HUGO IDs are duplicated in the raw data:
##
## Xraw$HUGO %>% keep( duplicated(.) )
##  [1] "TMSB15B" "PINX1"   "EMG1"    "COG8"    "H2BFS"
##
## 1. TMSB15B entries have 0.6 correlation and ~the same [low] total count
## 2. PINX1 entries have almost no correlation. Total #counts is 3,349 and 642.
## 3. EMG1 entries have -0.2 correlation. Total #counts is 4,361 and 1,944.
## 4. COG8 entries have almost no correlation. Total #counts is 19,274 and 2,097.
## 5. H2BFS has an entry with no counts.
##
## Given this, we opt to resolve duplicates in the following way
## 1. Average TMSB15B expression across the two duplicates.
## 2. Drop the lower-count duplicate for PINX1, EMG1, COG8, and H2BFS
cat( "Handling duplicates...\n" )
vDrop <- c("ENSG00000158427", "ENSG00000269226", ## TMSB15B before merging
           "ENSG00000258724",			 ## PINX1 with lower count (642)
           "ENSG00000268439",			 ## EMG1 with lower count (1,944)
           "ENSG00000272617",			 ## COG8 with lower count (2,097)
           "ENSG00000274559")			 ## H2BFS with no counts
X <- Xraw %>% filter( HUGO == "TMSB15B" ) %>% group_by(HUGO) %>%
    summarize_at( vars(-ensembl_id, -HUGO), mean ) %>% bind_rows( Xraw, . ) %>%
    filter( !(ensembl_id %in% vDrop) )

## Finalize the expression matrix
cat( "Log-transforming and finalizing expression matrix...\n" )
XX <- X %>% select( -ensembl_id ) %>% gather( SampleID, Value, -HUGO ) %>%
    mutate( Value = log2( Value+1 ) ) %>% spread( HUGO, Value )

## Load raw clinical covariates files
cat( "Processing clinical covariates...\n" )
YCBE <- synGet( "syn5223705", downloadLocation = local.dir )@filePath %>% read_csv( col_types = cols() )
YTCX <- (synGet( "syn3817650", downloadLocation = local.dir )@filePath) %>%
    read_csv( col_types = cols() ) %>% rename( SampleID = ID, Sex = Gender, Flowcell = FLOWCELL )
Yext <- (synGet( "syn14031984", downloadLocation = local.dir )@filePath) %>%
    read_tsv( col_types = cols() ) %>% rename( SampleID = NETdbID )

## Combine everything into a single clinical matrix
YY <- bind_rows(YCBE,YTCX) %>% select( -Source, -AgeAtDeath, -Tissue, -Sex, -Diagnosis ) %>%
    inner_join( Yext, by="SampleID" ) %>% select( -RIN, -Flowcell ) %>%
    select( SampleID, PMI, AOD = AgeAtDeath, Diagnosis = RNADiagnosis, Braak = `Braak stage`,
           Region, Thal = `Thal phase`, everything() ) %>% filter( !is.na(Braak) )

## Merge everything into a single data frame
cat( "Finalizing...\n" )
XY <- inner_join( YY, XX, by="SampleID" ) %>% rename( ID = SampleID )

## Write to file
fnOut <- file.path( local.dir, "mayo-pc.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
write_tsv( XY, fnOut )
