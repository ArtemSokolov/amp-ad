## Direct differential expression analysis for A-v-B, A-v-C, and B-v-C
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()

syn <- synExtra::synDownloader( "/data/AMP-AD/ROSMAP" )

## Composes a dichotomy over specified labels v
## The first label in the vector is assumed to correspond to the negative class
makeDichotomy <- function( v, Z )
{
    Z %>% filter( label %in% v ) %>%
        mutate_at( "label", factor, levels=v ) %>%
        select( rnaseq_id, label ) %>% as.data.frame() %>%
        column_to_rownames( "rnaseq_id" )
}

## Applies edgeR to the counts matrix in X, using dichotomy definition Y
myEdgeR <- function( X, Y )
{
    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X, samples = Y )
    cat( "Computing normalization factor...\n" )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~label, data = dl$samples )
    cat( "Estimating dispersion...\n" )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression
    cat( "Fitting a GLM...\n" )
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    RR <- edgeR::topTags( gf, nrow(X) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
    RR
}

main <- function()
{
    ## Define a mapping from Braak score to classification label
    M <- tibble( braaksc = 0:6, label = c("A","A","A","B","B","C","C") )

    ## Download raw (non-log2-normalized) expression data
    ##   and matching metadata matrices
    X <- syn( "syn3505720" ) %>% read_tsv( col_types=cols() ) %>%
        select( -tracking_id, -contains("492_120515") ) %>%
        rename_at( vars(-gene_id), str_sub, 1, -3 )
    Y <- syn( "syn3191087" ) %>% read_csv( col_types=cols() ) %>%
        select( projid, braaksc ) %>% inner_join(M)
    Z <- syn( "syn3382527" ) %>% read_csv( col_types=cols() ) %>%
        select( projid, rnaseq_id ) %>% filter( !is.na(rnaseq_id) ) %>%
        distinct() %>% inner_join(Y) %>% filter( !(rnaseq_id == "492_120515") )

    ## Compose the three dichotomies
    Zs <- list( AvB = c("A", "B"), BvC = c("B", "C"), AvC = c("A", "C") ) %>%
        map( makeDichotomy, Z )
    Xs <- map( Zs, ~select( X, gene_id, one_of(rownames(.)) ) ) %>%
        map( as.data.frame ) %>% map( column_to_rownames, "gene_id" )

    ## Compute differential gene expression
    DFXs <- map2( Xs, Zs, myEdgeR )
    ##    save( DFXs, file="temp.RData" )

    ## Annotate with gene names
    GM <- syn( "syn14236139" ) %>% read_csv( col_types=cols() ) %>% select( gene_id, gene_name )
    RR <- map( DFXs, mutate, gene_id = str_split( Gene, "\\.", simplify=TRUE )[,1] ) %>%
        map( left_join, GM, by="gene_id" ) %>% map( select, -Gene ) %>%
        map( select, gene_id, gene_name, everything() )

    ## Write out individual files
    RR1 <- imap( RR, ~write_csv( .x, str_c("dfx-", .y, ".csv.gz") ) )
}
