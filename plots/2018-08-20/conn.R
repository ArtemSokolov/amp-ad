## Evaluates the relationship between gene set performance and inter-connectivity
##
## by Artem Sokolov

source( "../../R/lpocv.R" )
source( "api.R" )

synapseLogin()

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

## Computes average diffusion score for the generated feature sets
## Stores results to local file
dfFeatSets <- function()
{
    ## Load the feature setss and PathwayCommons diffusion matrix
    load( "featSets250.RData" )
    load( "~/data/PathwayCommons/v10/dfsn.RData" )

    ## Reduce everything to a common set of gene names
    rr <- map( S, intersect, rownames(Df) ) %>% map_dbl( ~sum(Df[.x,.x]) )
    enframe( rr, "Tag", "Dfsn" ) %>% write_csv( "fs250-dfsn.csv" )
}

main <- function()
{
    ## Load all the relevant files
    RR <- c("fs250-cor.csv", "fs250-perf.csv", "fs250-dfsn.csv") %>% str_c( "data/", . ) %>%
        map( ~read_csv(.x, col_types=cols()) ) %>% plyr::join_all(by="Tag")

    gg1 <- ggplot( RR, aes( x=PWCor, y=AUC ) ) + theme_bw() + bold_theme() +
        geom_point() + xlab( "Average Pairwise Correlation" )
    gg2 <- ggplot( RR, aes( x=Dfsn, y=AUC ) ) + theme_bw() + bold_theme() +
        geom_point() + xlab( "Total Diffusion on PathwayCommons Network" )
    gridExtra::arrangeGrob( gg1, gg2, nrow=1 ) %>% ggsave( "plots/conn.png", ., width=10, height=5 )
}
