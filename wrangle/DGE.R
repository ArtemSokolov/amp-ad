## Wrangles aligned counts table into a genes-by-well matrix
##
## by Artem Sokolov

library( tidyverse )
library( stringr )
library( synapseClient )

synapseLogin()

## Composes a mapping between ENSEMBL IDs and HUGO names
ens2hugo <- function()
{
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    tx <- ensembldb::transcripts( edb, column=c("gene_id", "gene_name") )
    data_frame( HUGO = tx$gene_name, ENSEMBL = tx$gene_id ) %>% distinct
}

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/DGE" )@filePath }

## Consolidates "Raw" concentration values in the origin spreadsheet into a set of common categories
map.conc <- function( v )
{
    setNames( c("0.3uM", "10uM", "1uM", "1uM", "10uM", "3uM", NA),
             c("0.299991", "9.99000001", "1_M", "0.99990001", "10_M", "2.932453597", "2_g/ml") )[v]
}

main <- function()
{
    ## Load the matrix, row names, column names and the well tag-ID map
    X <- syn( "syn11947013" ) %>% read_delim( " ", comment="%", col_types=cols() ) %>%
        select( rowIndex = 1, colIndex = 2, Value = 3 )
    rn <- syn( "syn11947015" ) %>% read_csv( col_names = "ENSEMBL", col_types=cols() ) %>%
        mutate( rowIndex = 1:34947 ) %>% inner_join( ens2hugo() ) %>% filter( !duplicated(HUGO) )
    cn <- syn( "syn11947014" ) %>% read_csv( col_names = "barcode", col_types=cols() ) %>%
        mutate( colIndex = 1:384 )
    wmap <- syn( "syn11947016" ) %>% read_tsv( col_types=cols() )

    ## Pull everything together into a single genes-by-wells matrix
    XX <- inner_join( X, rn ) %>% inner_join( cn ) %>% inner_join( wmap ) %>%
        select( HUGO, well, Value ) %>% spread( well, Value, fill=0L )
    XX %>% write_csv( "wrangled-counts.csv" )

    ## Load the metadata
    Y <- syn( "syn11947020" ) %>% read_tsv( col_types=cols() ) %>%
        select( Well = `Dispensed well`, Drug = `Fluid name`, RawConc = Concentration ) %>%
        mutate( Concentration = map.conc(RawConc) ) %>% select( -RawConc )
    Y %>% write_csv( "wrangled-meta.csv" )
}
