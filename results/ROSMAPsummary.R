## Composes a summary of performance
## NOTE: This is on old (non-pc) ROSMAP data and needs to be updated
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
synapseLogin()

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD" )
{ synGet( id, downloadLocation = dlc )@filePath }

## Summary list of performances for Steve
list4steve <- function()
{
    gm_mean = function(x, na.rm=TRUE){
        exp(mean(log(x)))
    }
    
    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( LINCSID = lincs_id, URL = link, Name = name, Target = target_name )

    ## Load all the relevant stats files
    XX <- map( c( AB = "syn12959508", AC = "syn12958954", BC = "syn12961416" ), syn ) %>%
        map( read_csv, col_types=cols() ) %>% bind_rows( .id = "Task" ) %>%
        select( URL = id, n_genes, p_value, Task )

    ## Join all the data into a single matrix and compute average performance
    RR <- inner_join( XX, MLINCS ) %>% select( -URL ) %>% spread( Task, p_value ) %>%
        rowwise %>% mutate( ave_pval = gm_mean( c(AB,AC,BC) ) ) %>% ungroup %>%
        select( Name, Target, ave_pval, LINCSID, n_genes,
               AB_pval = AB, AC_pval = AC, BC_pval = BC )

    arrange( RR, ave_pval ) %>% write_csv( "perf-summary.csv" )
}
