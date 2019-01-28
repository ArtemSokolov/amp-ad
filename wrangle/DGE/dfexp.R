## Differential expression computed from DGE datasets
##
## by Artem Sokolov

library( tidyverse )

## Synapse interface
synapser::synLogin()
syn <- synExtra::synDownloader("~/data/AMP-AD/DGE")

## edgeR analysis applied to drug-vs-control dichotomy
dvc_edgeR <- function( X, drugWells, ctrlWells, drugName )
{
    cat( "Comparing", drugName, "against controls\n" )

    Y1 <- bind_rows( data_frame(Drug = drugName, Well = drugWells),
                    data_frame(Drug = "Control", Well = ctrlWells) ) %>%
        mutate_at( "Drug", factor, levels=c("Control", drugName) ) %>%
        as.data.frame %>% column_to_rownames("Well")
    X1 <- X %>% select( HUGO, rownames(Y1) ) %>% as.data.frame %>% column_to_rownames("HUGO")

    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X1, samples = Y1 )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Drug, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    RR <- edgeR::topTags( gf, nrow(X1) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
    RR
}

## Removes wells that have fewer than threshold counts (total, across all genes) in them
## Returns filtered Y data frame
filterWells <- function( Y, X, thresh = 1e5 )
{
    wKeep <- X %>% summarize_at( vars(-HUGO), sum ) %>% gather( Well, TotalCounts ) %>%
        filter( TotalCounts >= thresh ) %>% pull(Well)
    filter( Y, Well %in% wKeep )
}

## Given a drug name, identifies the highest concentration that has at least k replicates
## Returns the associated well IDs
drug2Wells <- function( Y, drugName, k=2 )
{
    ## Determine all concentrations that have at least two replicates
    vConc <- Y %>% filter( Drug == drugName ) %>% group_by( Concentration ) %>%
        summarize( nn = n() ) %>% filter( nn >= k ) %>% pull( Concentration )
    if( length(vConc) == 0 )
        stop( drugName, " has no concentration with at least ", k, " replicates" )

    ## Select the top concentration and retrieve the associated well IDs
    Y %>% filter( Drug == drugName, Concentration %in% vConc ) %>%
        filter( Concentration == max(Concentration) ) %>% pull(Well)
}

## Differential expression on DGE1 experiment
DGE1 <- function()
{
    ## Load the data
    ## Identify drugs that are toxic by looking at the total number of counts in each well
    ## Remove them from consideration
    X <- syn( "syn15673461" ) %>% read_csv( col_types=cols() )
    Y <- syn( "syn15673460" ) %>% read_csv( col_types=cols() ) %>% filterWells( X )

    ## Identify control wells
    wCtrl <- Y %>% filter( Drug == "Drug control" ) %>% pull(Well)
    
    ## Compose the set of drugs to consider and map them to the corresponding wells
    ## Exclude Chitosan as it is highly toxic
    wDrugs <- Y %>% filter( !is.na(Concentration) ) %>% pull(Drug) %>% unique %>%
        setdiff( "Chitosan" ) %>% set_names %>% map( ~drug2Wells(Y, .) )
    
    ## Traverse the drugs and compute differential expression for each
    RR <- imap( wDrugs, ~dvc_edgeR( X, .x, wCtrl, .y ) )
    bind_rows( RR, .id = "Drug" ) %>% write_csv( "DGE1-dfx.csv.gz" )
}

## Differential expression on DGE2 experiment
DGE2 <- function()
{
    ## Load the data
    ## Use a slightly lower threshold for well filtering, since the overall counts distribution
    ##   is shifted to the left (See DGE2 QC vignette)
    X <- syn( "syn17115624" ) %>% read_csv( col_types=cols() )
    Y <- syn( "syn17115625" ) %>% read_csv( col_types=cols() ) %>% filterWells( X, 8.5e4 )

    ## Identify control wells
    wCtrl <- Y %>% filter( Drug == "DMSO" ) %>% pull(Well)

    ## Compose the set of drugs to consider
    wDrugs <- Y %>% filter( !is.na(Concentration) ) %>% pull(Drug) %>% unique %>%
        discard( ~grepl("dsRNAmi|Lipo", .) ) %>% set_names %>% map( ~drug2Wells(Y, .) )

    ## Traverse the drugs and compute differential gene expression for each
    RR <- imap( wDrugs, ~dvc_edgeR( X, .x, wCtrl, .y ) )
    bind_rows( RR, .id = "Drug" ) %>% write_csv( "DGE2-dfx.csv.gz" )
}

