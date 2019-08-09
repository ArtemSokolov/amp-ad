## Additional tests for TF analysis
##
## by Artem Sokolov

library( tidyverse )

## Define Synapse downloader(s)
synapser::synLogin()
syn <- synExtra::synDownloader("./data", ifcollision="overwrite.local")
syn_csv <- function( id ) {syn(id) %>% read_csv(col_types=cols())}
syn_tsv <- function( id ) {syn(id) %>% read_tsv(col_types=cols())}

## Load the full matrix of results
loadResults <- function()
{
    syn_csv( "syn20555046" ) %>%
        rename( DrugTgt = genetarget, TF=pathway ) %>%
        select( -DGE, -fTox, -fConc, -mpi, -leadingEdge ) %>%
        mutate_at( "DrugTgt", str_to_upper )
}

main <- function()
{
    X <- loadResults()
    
    ## Identify the top 10 TF for each drug target
    R <- X %>% nest( -DrugTgt, .key=Res ) %>%
        mutate_at( "Res", map, function( .df, k=10 )
        {
            .df %>% arrange(pval) %>% head(k) %>%
                select(TF) %>% mutate( Rank = 1:k )
        }) %>%
        unnest() %>% spread( DrugTgt, TF )

    ## Check for RFX5 / RFXANK presence (DGE1/2 batch effect)
    R %>% select_if( ~any(grepl("RFX",.x)) )
}

viz <- function()
{
    X <- loadResults()
    SIK3 <- X %>% filter( DrugTgt == "SIK3" ) %>% arrange(pval)
    
    source( "TF_Visualizations_v2.R" )
    fplot <- TF_enriched_targets_Violin_Plot
    fplot( drug_target="SIK3", TF_name="RFX5" )
    fplot( drug_target="SIK3", TF_name="NF1" )
}
