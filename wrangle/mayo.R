## Wrangling of Mayo clinic data and matching clinical annotations
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
suppressMessages(library( synapseClient ))

## Composes a mapping between ENSEMBL IDs and HUGO names
ens2hugo <- function()
{
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    tx <- ensembldb::transcripts( edb, column=c("gene_id", "gene_name") )
    data_frame( HUGO = tx$gene_name, ENSEMBL = tx$gene_id ) %>% distinct
}

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

## Read raw expression matrices (one for cerebellum, one for temporal cortex)
cat( "Downloading expression data...\n" )
ct <- cols( .default = col_double(), ensembl_id = col_character() )
XCBE <- synGet( "syn5201007", downloadLocation = local.dir )@filePath %>% read_tsv( col_types = ct )
XTCX <- synGet( "syn4650265", downloadLocation = local.dir )@filePath %>% read_tsv( col_types = ct )

## Combine the two expression matrices
## 742_CER has a low number of counts; drop it
## Observation made by the following:
##   XCBE %>% summarize_at( -1, sum ) %>% gather( Sample, Total ) %>% arrange( Total )
cat( "Mapping gene IDs to HUGO...\n" )
Xraw <- XCBE %>% select( -`742_CER` ) %>% inner_join( XTCX, by="ensembl_id" ) %>%
    rename( ENSEMBL = ensembl_id ) %>% inner_join( ens2hugo(), ., by="ENSEMBL" )

## Remove alternative splice forms and non-coding RNA
cat( "Removing alternative splice forms, non-coding RNA and duplicate entries...\n" )
f <- function( x, pattern ) { filter( x, !grepl(pattern, HUGO) ) }
vBL <- c("Y_RNA", "Metazoa_SRP", "Vault", "5S_rRNA", "uc_338", "7SK",
         "mascRNA-menRNA", "U4atac", "U6atac")
X <- Xraw %>% f( "\\." ) %>% f( "^MIR" ) %>% f( "^RNU" ) %>% f( "^SNOR" ) %>%
    f( "^U[1-9]$" ) %>% f( "^SCARNA" ) %>% f( "^sno" ) %>% f( "^LINC" ) %>%
    f( "-AS[1-9]$" ) %>% f( "^ACA[1-9]" ) %>% filter( !(HUGO %in% vBL) )

## Merge duplicate gene entries
cat( "Consolidating duplicate gene entries...\n" )
vdup <- X$HUGO %>% keep( duplicated(.) ) %>% unique
Xdedup <- X %>% filter( HUGO %in% vdup ) %>% select( -ENSEMBL ) %>%
    gather( SampleID, Value, -HUGO ) %>% group_by( HUGO, SampleID ) %>%
    summarize( Value = median(Value) ) %>% ungroup

## Finalize the expression matrix
cat( "Log-transforming and finalizing expression matrix...\n" )
XX <- X %>% filter( !(HUGO %in% vdup) ) %>% select( -ENSEMBL ) %>%
    gather( SampleID, Value, -HUGO ) %>% bind_rows(Xdedup) %>%
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
fnOut <- file.path( local.dir, "mayo-wrangled.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
write_tsv( XY, fnOut )
