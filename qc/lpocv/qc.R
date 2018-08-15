## Verification that one-off LPOCV code produces results
##  that are consistent with the BTR gold standard
##
## by Artem Sokolov

source( "../../R/lpocv.R" )

## Evaluates a single task in an LPOCV setting
evalTask <- function( task )
{
    ## Load the data and the pair assignments for a given task
    cat( "Loading data...\n" )
    X <- loadROSMAP() %>% reduceToTask( task )
    lP <- loadPairs( task )

    ## Sample random sets of genes of varying sizes
    cat( "Selecting feature sets...\n" )
    vExclude <- c("ID", "PMI", "AOD", "CDR", "Braak", "BrodmannArea", "Barcode", "Label")
    vFullSet <- setdiff( colnames(X), vExclude )
    XX <- seq( 100, 1000, by=100 ) %>% set_names %>% map( ~sample(vFullSet, .x) ) %>%
        map( ~reduceToGenes(X, .) )
    
    ## Run cross-validation on each subset
    imap( XX, ~lpocv(.x, lP, .y) ) %>% bind_rows %>% mutate( Task = task )
}

## Evaluates all tasks and caches the results to a local file
evalTasks <- function()
{
    ## Use a fixed random seed to allow reproducibility
    set.seed(1)

    ## Evaluate each task separately
    RR <- map( c("AB", "AC", "BC"), evalTask )

    ## Store results to a local file
    bind_rows( RR ) %>% rename( Size = Tag ) %>% write_csv( "lpocv-qc.csv" )
}


## Plots a comparison of BTR results against those from the currect code
plotComparison <- function()
{
    ## Load all the relevant BTR results
    X <- data_frame( Task = c("AB", "AC", "BC"), Framework = "BTR",
                    synid = c("syn15589822", "syn15589816", "syn15589810") ) %>% 
        mutate( Values = map(synid, ~read_csv(syn(.x), col_types = cols())) )

    ## Load the locally-cached results
    Y <- read_csv( "lpocv-qc.csv", col_types = cols() ) %>% mutate( Framework = "One-Off" )

    ## Combine everything into a common data frame
    RR <- X %>% mutate( Values = map( Values, gather, Size, AUC ) ) %>% unnest %>%
        mutate_at( "Size", as.integer ) %>% select( -synid )

    ## Plot the comparison
    ggplot( RR, aes(x=Size, y=AUC) ) + theme_bw() +
        facet_wrap( ~Task ) + geom_smooth( se=FALSE ) +
        geom_point( data = Y )
}

