## Generation of synthetic data
##
## by Artem Sokolov

library( tidyverse )
library( stringr )

## ###### ##
## CONFIG ##
## ###### ##

## Random number generation seed to allow for reproducibility
set.seed( 100 )

## Preset the number of features and samples
nFeats <- list( ABCeasy = 100,		## 10% overlap aligned with A, B, C
            ABChard = 100,		## 50% overlap aligned with A, B, C 
            BACeasy = 100,		## 10% overlap aligned with B, A, C
            NoSignal = 700 )		## No separation
nSamples <- c( A = 200, B = 100, C = 100 )

## ######### ##
## GENERATOR ##
## ######### ##

## Generates a Gaussian cluster of n samples by first sampling (uniformly at random) the mean out of
##   rngMean and st.d. out of rngStD, and then using those parameters for the Gaussian generation
## Final values are offset to ensure non-negativity
rnorm1 <- function( n, rngMean = c(0,10), rngStD = c(0.1,2) )
{
    m <- runif( 1, rngMean[1], rngMean[2] )
    s <- runif( 1, rngStD[1], rngStD[2] )
    x <- rnorm( n, m, s )
    if( min(x) < 0 )
        x <- x - min(x)
    x
}

## Lines up vector x at the desired level of overlap into the upper region of vector y
## fOverlap - fraction of y that is overlapped by x. Must be in [0,1]
vLineup <- function( x, y, fOverlap )
{
    stopifnot( between( fOverlap, 0, 1 ) )
    x - min(x) + quantile( y, 1 - fOverlap )
}

## Generates three Gaussian clusters with the desired level of overlap
## nn - a vector (or list) of length 3, containing the number of samples in each cluster
## fOverlap - desired level of overlap between classes. See vLineup() for more details
## Returns a list of length 3, with each element containing a vector of values for the
##    corresponding cluster
rnormABC <- function( nn, fOverlap )
{
    ## Generate each cluster independently
    xx <- lapply( nn, rnorm1 )

    ## Line up the second cluster at the desired level of overlap with the first
    xx[[2]] <- vLineup( xx[[2]], xx[[1]], fOverlap )

    ## Similarly, line up the third cluster at the desired level of overlap with the second
    xx[[3]] <- vLineup( xx[[3]], xx[[2]], fOverlap )
    xx
}

## Three Gaussian clusters aligned with B, A, C
## nn - a vector (list) that is STILL in the A, B, C order
##  the function uses the second entry of nn to define the first cluster (which is B)
## Same output as rnormABC()
rnormBAC <- function( nn, fOverlap )
{
    xx <- rnormABC( nn[c(2,1,3)], fOverlap )
    xx[c(2,1,3)]
}

## Generates nF features using a predefined generator function
## nF - total number of features
## nn - a vector (or list) of length 3, containing the number of samples in each cluster
## fgen - generator function with same output as what's produced by rnormABC()
## pfx - prefix to use for feature names
genFeats <- function( nF, nn, fgen, pfx )
{
    ## Build vectors for individual features
    vv <- list(nn) %>% rep(nF) %>% lapply( fgen ) %>% lapply( unlist )

    ## Combine everything into a common data.frame and annotate accordingly
    do.call( cbind, vv ) %>% magrittr::set_colnames( str_c( pfx, 1:nF ) ) %>%
        as.data.frame %>% rownames_to_column( "ID" )
}

## Samples from the Braak-CDR bivariate distribution for a given range of Braak scores
sampleBraakCDR <- function( n, brMin, brMax )
{
    source( "braak-cdr.R" )
    B <- distBraakCDR %>% mutate( prob = nSamples / sum(nSamples) ) %>%
        filter( between( Braak, brMin, brMax ) )
    sample( 1:nrow(B), n, replace=TRUE, prob=B$prob ) %>% slice( B, . ) %>% select( CDR, Braak )
}

## Generates a full feature / meta matrix for a given #samples and #features
## nnF - profile for the number of features (see nFeatures global setting at the top)
## nnS - profile for the number of samples (see nSamples global setting at the top)
main <- function( nnF, nnS )
{
    ## Aligned with ABC, 10% overlap (easy separation)
    F1 <- genFeats( nnF$ABCeasy, nnS, function(.) {rnormABC(.,0.1)}, "ABCeasy" )

    ## Aligned with ABC, 50% overlap (hard separation)
    F2 <- genFeats( nnF$ABChard, nnS, function(.) {rnormABC(.,0.5)}, "ABChard" )

    ## Aligned with BAC, 10% overlap
    F3 <- genFeats( nnF$BACeasy, nnS, function(.) {rnormBAC(.,0.1)}, "BACeasy" )

    ## No signal (fNS() has same output as rnormABC(), but uses rnorm1() instead)
    fNS <- function(.) { split( rnorm1(sum(.)), rep(names(.),.) ) }
    F4 <- genFeats( nnF$NoSignal, nnS, fNS, "NS" )

    ## Join all the feature frame together
    FF <- plyr::join_all( list(F1,F2,F3,F4) )

    ## Generate some meta data
    BrCDR <- bind_rows( sampleBraakCDR( nnS[[1]], 0, 2 ),
                       sampleBraakCDR( nnS[[2]], 3, 4 ),
                       sampleBraakCDR( nnS[[3]], 5, 6 ) )
    MM <- select( FF, ID ) %>%
        mutate( PMI = sample( 75:540, sum(nnS), replace=TRUE ),
               AOD = sample( c( as.character(60:89), rep("90+",15) ), sum(nnS), replace=TRUE ) ) %>%
        bind_cols( BrCDR ) %>%
        mutate( BrodmannArea = sample( c("BM22","BM36"), sum(nnS), replace=TRUE ),
               Barcode = sample( 1:sum(nnS) ) + 1000 )

    ## Join everything into a single data.frame
    inner_join( MM, FF )
}

main( nFeats, nSamples ) %>% write_tsv( "synthetic.tsv" )
