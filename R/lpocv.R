## Leave-pair-out cross-validation for one-off experiments on ROSMAP
## For large-scale analyses, use the BTR framework: https://github.com/pvtodorov/btr
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
library( LiblineaR )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/one-off" )@filePath }

## Downloads ROSMAP pair assignments associated with a given task
## Returns a list of pairs that can be supplied directly to lpocv()
loadPairs <- function( task )
{
    P <- syn( "syn15665557" ) %>% read_csv(col_types = cols(ID=col_character())) %>%
        filter( Task == task )
    with( P, split( ID, pair_index ) )
}

## Loads the ROSMAP dataset
loadROSMAP <- function()
{ Xraw <- syn( "syn14306482" ) %>% read_tsv(col_types = cols(ID = col_character())) }

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
## X - Dataset, as loaded by loadROSMAP() or produced by reduceToTask()
## gs - vector character containing the gene set of interest
reduceToGenes <- function( X, gs )
{
    ## Argument verification
    stopifnot( all(gs %in% colnames(X)) )

    ## Reduce columns accordingly
    X %>% select( ID, Label, one_of(gs) )
}

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

## Computes AUC from LPOCV
## hi - prediction values associated with high-label example in the pair
## lo - prediction values associated with low-label example in the pair
AUC_LPOCV <- function( hi, lo )
{ (0.5*sum(hi==lo) + sum(hi>lo)) / length(hi) }

## Leave pair out cross-validation for a given dataset
## XY - dataset must contain columns ID (denoting sample identifiers) and Label
## lPairs - list of vectors-of-length-2 specifying IDs to withhold
## tag - a tag for the given run
lpocv <- function( XY, lPairs, tag )
{
    cat( "Running LPOCV for", tag, "\n" )

    ## Ensure that only pairs of samples are withheld
    stopifnot( all( range(map_int( lPairs, length )) == c(2,2) ) )
    stopifnot( all( flatten_chr(lPairs) %in% XY$ID ) )

    ## Traverse the pairs and perform cross-validation
    RR <- map( lPairs, ~liblinear(XY, .x) )

    ## Compute AUC
    bind_rows( RR, .id="index" ) %>% select( -ID ) %>% spread( Label, Pred ) %>%
        summarize( AUC = AUC_LPOCV( pos, neg ) ) %>% mutate( Tag = tag )
}

