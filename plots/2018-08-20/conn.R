## Evaluates the relationship between gene set performance and inter-connectivity
##
## by Artem Sokolov

source( "../../R/lpocv.R" )

## Generates a collection of feature sets of the desired size
## Stores the collection to a local file
genFeatSets <- function()
{
    ## Preset parameters
    set.seed(1)
    setSize <- 250
    nSets <- 100
    fnOut <- "featSets250.RData"

    ## Load the data and identify the feature space to sample
    X <- loadROSMAP()
    vExclude <- c("ID", "PMI", "AOD", "CDR", "Braak", "BrodmannArea", "Barcode", "Label")
    vFullSet <- setdiff( colnames(X), vExclude )

    ## Sample nSets of size setSize
    ## Save result to a file
    S <- str_c( "Set", 1:nSets ) %>% set_names %>% map( ~sample(vFullSet, setSize) )
    save( S, file=fnOut )
}

## Applies cross-validation to every feature set using ROSMAP AC
## Stores the results to a local file
evalFeatSets <- function()
{
    ## Load the data, pair assignments and feature sets
    X <- loadROSMAP() %>% reduceToTask( "AC" )
    lP <- loadPairs( "AC" )
    load( "featSets250.RData" )

    ## Apply cross-validation on each set
    XX <- map( S, ~reduceToGenes(X, .) )
    RR <- imap( XX, ~lpocv(.x, lP, .y) )

    bind_rows(RR) %>% write_csv( "fs250-perf.csv" )
}

## Computes average pair-wise correlation for a given subset of columns
avePWCor <- function( Z )
{
    ## Identify and remove columns with no variance
    Z <- as.matrix(Z)
    v <- apply( Z, 2, sd ) %>% keep( ~.x > 1e-5 ) %>% names
    (sum(cor(Z[,v])) - length(v)) / (length(v)^2 - length(v))
}

## Computes pair-wise correlation statistics for the generated feature sets
## Stores the results to a local file
corFeatSets <- function()
{
    ## Load the data and feature sets
    X <- loadROSMAP() %>% reduceToTask( "AC" )
    load( "featSets250.RData" )

    ## Downsample the data to the pre-selected feature sets
    XX <- map( S, ~select(X, one_of(.)) )

    ## Compute average pairwise correlation for each subset
    map_dbl( XX, avePWCor ) %>% data_frame( PWCor = ., Tag = names(.) ) %>% write_csv( "fs250-cor.csv" )
}

main <- function()
{
    ## Load correlation and performance measures
    RR <- inner_join( read_csv( "fs250-cor.csv", col_types = cols() ),
                     read_csv( "fs250-perf.csv", col_types = cols() ) )

    ggplot( RR, aes( x=PWCor, y=AUC ) ) + theme_bw() + geom_point()
}
