## Wrangling data from DGE1 plate that was re-sequenced
##  with the new protocol
##
## by Artem Sokolov

library( tidyverse )

## Set up the Synapse Downloader
synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/AMP-AD/DGE1_new", ifcollision="overwrite.local" )

preQC <- function()
{
    ## Load the mapping of well barcodes to well IDs
    wmap <- syn( "syn17103737" ) %>% read_tsv( col_types=cols() ) %>% select( Well=well, Barcode=barcode )

    ## Load the mapping of ENSEMBL gene IDs to HUGO gene names
    gmap <- syn( "syn14236139" ) %>% read_csv( col_types=cols() ) %>%
        select( gene_id, gene_name, gene_biotype )

    ## Load the counts table (in long format)
    X <- syn( "syn18143723" ) %>% read_delim( " ", comment="%", col_types=cols() ) %>%
        rename( rowIndex=1, colIndex=2, Value=3 )

    ## Load the matching row and column names
    rn <- syn("syn18143724") %>% read_csv( col_names="gene_id", col_types=cols() ) %>%
        mutate( rowIndex=1:34947 )
    cn <- syn("syn18143725") %>% read_csv( col_names="Barcode", col_types=cols() ) %>%
        mutate( colIndex=1:386 ) %>% left_join( wmap, by="Barcode" ) %>% select( -Barcode )

    ## Data contains several duplicated gene_names
    ## For each gene name, we inspected the corresponding gene IDs and kept the one
    ##   with the highest number of total reads that map to it
    vDrop <- c( "ENSG00000213380", "ENSG00000268439", "ENSG00000234289", "ENSG00000284024",
               "ENSG00000280987", "ENSG00000284741", "ENSG00000258724", "ENSG00000206549",
               "ENSG00000130489", "ENSG00000158427" )
    
    ## Pull everything together into a single genes-by-wells matrix
    ## Map gene IDs to HUGO names and reduce to protein-coding space
    XX <- inner_join( X, rn, by="rowIndex" ) %>% inner_join( cn, by="colIndex" ) %>%
        select( Well, gene_id, Value ) %>% inner_join( gmap, by="gene_id" ) %>%
        filter( gene_biotype == "protein_coding", !(gene_id %in% vDrop) ) %>%
        select( -gene_biotype, -gene_id ) %>% spread( Well, Value, fill=0L ) %>%
        select( -P11_new2, -P11_old ) %>% rename( P11 = P11_new1, HUGO = gene_name )

    XX %>% write_csv( "wrangled-counts.csv" )
}
