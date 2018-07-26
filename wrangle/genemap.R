## Creates a map from ENSEMBL IDs to HUGO for protein-coding regions
##
## by Artem Sokolov

library( tidyverse )

## Location of the latest .gtf release
fnGTF <- "ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz"
fnLocal <- "GRCh38.93.gtf.gz"

main <- function()
{
    ## Download the raw .gtf file
    download.file( fnGTF, fnLocal )

    ## Parse the raw .gtf file
    cn <- c( "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
    ct <- cols( .default = col_character(), start = col_integer(), end = col_integer() )
    X <- read_tsv( fnLocal, comment="#!", col_types=ct, col_names = cn )

    ## Isolate the gene entries and parse key-value pais for each
    vkv <- X %>% filter( feature == "gene" ) %>% .$attribute
    vv <- vkv %>% str_sub(1,-2) %>% str_split( "; " ) %>% map( str_split, " " )
    G <- vv %>% modify_depth( 2, ~set_names(.x[2], .x[1]) ) %>% map( flatten_dfc ) %>%
        bind_rows %>% mutate_all( str_sub, 2, -2 )

    ## Save a copy to a file
    G %>% write_csv( "GRCh38.93.genemap.csv" )
}
