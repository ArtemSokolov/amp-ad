## Quick-and-dirty plots for one-off gene set runs on ROSMAP
## For large-scale analyses, use the BTR framework: https://github.com/pvtodorov/btr
##
## by Artem Sokolov

library( tidyverse )
library( LiblineaR )
library( synapser )

synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/one-off" )
g_ct <- cols(ID=col_character())

## Train-test for a single pair using liblinear implementation
liblinear <- function( XY, vTest )
{
    ## Argument verification
    stopifnot( all( sort(unique(XY$Label)) == c("neg","pos") ) )
    
    ## Split the data into train and test
    Xte <- XY %>% filter( ID %in% vTest ) %>% select( -ID, -Label ) %>% as.matrix
    Xtr <- XY %>% filter( !(ID %in% vTest) ) %>% select( -ID, -Label ) %>% as.matrix
    ytr <- filter( XY, !(ID %in% vTest) )$Label

    ## Train a model and apply it to test data
    m <- LiblineaR( Xtr, ytr, type=0 )
    ypred <- predict( m, Xte, proba=TRUE )$probabilities[,"pos"]
    XY %>% filter( ID %in% vTest ) %>% select( ID, Label ) %>% mutate( Pred = ypred )
}

## Leave pair out cross-validation for a given dataset
## XY - dataset must contain columns ID (denoting sample identifiers) and Label
## lPairs - list of vectors-of-length-2 specifying IDs to withhold
lpocv <- function( XY, lPairs )
{
    ## Computes AUC from LPOCV
    ## hi - prediction values associated with high-label example in the pair
    ## lo - prediction values associated with low-label example in the pair
    AUC_LPOCV <- function( hi, lo )
    { (0.5*sum(hi==lo) + sum(hi>lo)) / length(hi) }
    
    cat( "." )
    
    ## Ensure that only pairs of samples are withheld
    stopifnot( all( range(map_int( lPairs, length )) == c(2,2) ) )
    stopifnot( all( flatten_chr(lPairs) %in% XY$ID ) )

    ## Traverse the pairs and perform cross-validation
    RR <- map( lPairs, ~liblinear(XY, .x) )

    ## Compute AUC
    bind_rows( RR, .id="index" ) %>% select( -ID ) %>% spread( Label, Pred ) %>%
        with( AUC_LPOCV( pos, neg ) )
}

## Downloads ROSMAP pair assignments used by BTR, associated with a given task
## Returns a list of pairs that can be supplied directly to lpocv()
loadPairs <- function( task )
{
    syn( "syn15665557" ) %>% read_csv( col_types = g_ct ) %>%
        filter( Task == task ) %>% with( split( ID, pair_index ) )
}

## Loads the ROSMAP dataset
loadROSMAP <- function() { syn("syn14306482") %>% read_tsv(col_types = g_ct) }

## Reduces a dataset to the binary task of interest
## X - Dataset, as loaded by loadROSMAP()
## task - one of {"AB", "AC", "BC"}
reduceToTask <- function( X, task )
{
    ## Define the task mapping
    ## Note that the mapping has 1-based indexing
    ##   while Braak score ranges from 0 to 6
    taskMap <- list(
        AB = c("neg", "neg", "neg", "pos", "pos", NA, NA),
        AC = c("neg", "neg", "neg", NA, NA, "pos", "pos"),
        BC = c(NA, NA, NA, "neg", "neg", "pos", "pos")
    )
    
    ## Argument verification
    stopifnot( task %in% names(taskMap) )
    stopifnot( is.integer( X$Braak ) )

    ## Reduce rows accordingly
    ## Note the +1 to align 0-based and 1-based indexing
    X %>% mutate( Label = taskMap[[task]][Braak+1] ) %>% filter( !is.na(Label) )
}

## Reduces a dataset to the gene set of interest
## gs - vector character containing the gene set of interest
## X - Dataset, as loaded by loadROSMAP() or produced by reduceToTask()
reduceToGenes <- function( gs, X )
{
    ## Argument verification
    stopifnot( all(gs %in% colnames(X)) )

    ## Reduce columns accordingly
    X %>% select( ID, Label, one_of(gs) )
}

## Given a dataset and a gene set of interest, generates background sets of equal size
## gsi - gene set of interest
## X - Dataset, as loaded by loadROSMAP() or returned by reduceToTask()
## nBK - number of background sets to generate
## vExclude - set of identifiers to exclude from sampling
genBK <- function( gsi, X, nBK, vExclude=c("ID", "PMI", "AOD", "CDR", "Braak",
                                           "BrodmannArea", "Barcode", "Label") )
{
    ## Intersect the gene set of interest against dataset's feature space
    vFullSet <- setdiff( colnames(X), vExclude )
    vGSI <- intersect( gsi, vFullSet )
    nGSI <- length(vGSI)

    ## Sample random sets of genes of matching size
    seq(1, length.out=nBK) %>% map( ~sample(vFullSet, nGSI) ) %>%
        set_names( rep("BK", nBK) )
}

## Evaluates gene sets of interest in the context of a given dataset / task
## gsi - gene sets of interest (as a named list)
## X - Dataset, as loaded by loadROSMAP()
## task - one or more of {"AB", "AC", "BC"}
evalGeneSets <- function( gsi, X, task=c("AB","BC","AC") )
{
    ## Iterate over tasks if more than one is requested
    if( length( task ) > 1 )
    {
        RR <- map( task, ~evalGeneSets(gsi, X, .x) )
        return( bind_rows(RR) )
    }
    
    ## Reduce the dataset to the requested task and
    ##   identify the corresponding pair assignments
    cat( "Setting up to evaluate task", task, "...\n" )
    XY <- reduceToTask( X, task )
    lP <- loadPairs( task )

    ## Downsample the data according to the requested gene sets
    cat( "Generating data slices...\n" )
    SS <- `if`( is.list( gsi ), gsi, list(gsi) ) %>% enframe( "Name", "Feats" ) %>%
        mutate( Data = map(Feats, reduceToGenes, XY), Task = task )
    
    ## Run LPOCV on each slice of data
    cat( "Running LPOCV" )
    RR <- SS %>% mutate( AUC = map_dbl( Data, lpocv, lP ) )
    cat( "\n" )

    RR %>% select( -Data )
}

## Given a gene set of interest, automatically generates matching background
##   and evaluates the set against it
## gsi - gene set of interest
## nBK - number of background sets to generate
## X - Dataset, as loaded by loadROSMAP()
## task - one or more of {"AB", "AC", "BC"}
evalGeneSetBK <- function( gsi, nBK, X, task=c("AB","BC","AC") )
{
    genBK( gsi, X, nBK ) %>% c( "GSI" = list(gsi), . ) %>%
        evalGeneSets( X, task )
}

## Plots the results of evalGeneSetBK()
## RR - named list of outputs from evalGeneSetBK()
##   Names of list elements are used in the plot and must be specified
plotGeneSetEval <- function( Rl )
{
    ## Argument verification
    stopifnot( is.list(Rl) )
    stopifnot( !is.null(names(Rl)) )
    
    ## Combine everything into a single data.frame
    RR <- bind_rows( Rl, .id = "Tag" )
    
    ## Compound name that is used to separate out individual plots
    RR <- RR %>% mutate( GeneSet = glue::glue("{Tag}({Task})") ) %>%
        mutate_at( "GeneSet", as_factor )

    ## Separate out performance of background and foreground sets
    BK <- RR %>% filter( Name == "BK" )
    FG <- RR %>% filter( Name == "GSI" ) %>% mutate( nGenes = map(Feats, length) ) %>%
        mutate( Label = glue::glue( " {nGenes} genes" ) )

    ## Map gene set labels to the original gene set names
    ylbl <- arrange( FG, GeneSet ) %>% with( set_names(Tag, GeneSet) )
    
    ggplot( BK, aes(x=AUC, y=GeneSet, fill=Task) ) +
        ggridges::theme_ridges( center_axis_labels=TRUE ) +
        ggridges::geom_density_ridges2(scale=1.2, size=1, alpha=0.5) +
        ggthemes::scale_fill_few( "Medium" ) +
        scale_y_discrete( labels = ylbl ) +
        geom_segment( aes(y=as.numeric(GeneSet), yend=as.numeric(GeneSet)+0.9, 
                          x=AUC, xend=AUC), data = FG, color="tomato", lwd=2 ) +
        geom_text( data=FG, aes(x=AUC, y=as.numeric(GeneSet)+0.5, label=Label),
                  hjust=0, fontface="bold", size=4 )
}

## Testing the functionality in this file
mytest <- function()
{
    X <- loadROSMAP()

    ## PIK3CA inhibitor gsk1059615
    rGSK <- c("CSNK1G2", "PIK3CD", "NEK10", "DYRK1A", "PIK3CB", "PIK3CG", "JAK3", "BMP2K",
              "STK10", "MAP3K19", "MTOR", "GAK", "CLK1", "CLK3", "PIK3CA", "MAPK15") %>%
        evalGeneSetBK( 10, X )
    
    ## Apoptosis gene set provided by Kris Sarosiek
    rApoptosis <- c("BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "MCL1") %>%
        evalGeneSetBK( 10, X )

    ## SPARCS gene set provided by Russell Jenkins
    rSPARCS <- c("TRIM22", "TRIM38", "IL32", "SPATS2L", "EPHA3", "HERC3", "ADAM19", "SERPINB9",
                 "IFI44L", "F3", "BEND6", "AIG1", "MSRB2", "TNFRSF9", "ANTXR1" ) %>%
        evalGeneSetBK( 10, X )

    Y <- list( "GSK1059615" = rGSK, "Apoptosis" = rApoptosis, "SPARCS" = rSPARCS )
    plotGeneSetEval( Y ) ##+ ggsave( "mytest.pdf" )
}
