## Scripts that simplify loading of results for all figures
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/figs" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Pre-defined index of DGE results, constructed using the following code:
##   source( "../../R/resmine.R" )
##   synResults( "stats", Region != "cerebellum" ) %>% filter( grepl("DGE", name) ) %>%
##     select( -name, -parentName ) %>% dput()
indexDGE <- function()
{
    structure(list(Dataset = c("ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "MAYOpc", "MAYOpc", "MAYOpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TempCortex", "TempCortex",
                              "TempCortex", "BM10", "BM36", "BM44", "BM22", "BM36", 
                              "BM10", "BM22", "BM10", "BM44", "BM36", "BM44", "BM22"),
                   Task = c("AB", "AC", "BC", "BC", "AB", "AC", "AC", "AB", "BC",
                            "BC", "AC", "AB", "AB", "BC", "AC", "BC", "AB", "AC"),
                   id = c("syn16787402", "syn16787438", "syn16787423", "syn16787399",
                          "syn16787414", "syn16787417", "syn16787435", "syn16787447",
                          "syn16787408", "syn16787411", "syn16787444", "syn16787432", 
                          "syn16787405", "syn16787390", "syn16787420", "syn16787396",
                          "syn16787393", "syn16787387")),
              row.names = c(NA, -18L), class = c("tbl_df", "tbl", "data.frame"))
}

## Pre-defined index of results for Nienke's gene sets, built by the following code:
##   source( "../../R/resmine.R" )
##   synResults( "stats", Region != "cerebellum" ) %>% filter( grepl("Nienke", name) ) %>%
##     select( -name, -parentName ) %>% dput()
indexMined <- function()
{
    structure(list(Dataset = c("ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "MAYOpc", "MAYOpc", "MAYOpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TempCortex", "TempCortex",
                              "TempCortex", "BM10", "BM36", "BM44", "BM22", "BM36", 
                              "BM10", "BM22", "BM10", "BM44", "BM36", "BM44", "BM22"),
                   Task = c("AB", "AC", "BC", "BC", "AB", "AC", "AC", "AB", "BC",
                            "BC", "AC", "AB", "AB", "BC", "AC", "BC", "AB", "AC"),
                   id = c("syn15589825", "syn15589819", "syn15589813", "syn15572216",
                          "syn15571243", "syn15572107", "syn15584832", "syn15585721",
                          "syn15582897", "syn15585424", "syn15581498", "syn15582320", 
                          "syn15586156", "syn15585571", "syn15583851", "syn15584143",
                          "syn15583704", "syn15584689")),
              row.names = c(NA, -18L), class = c("tbl_df", "tbl", "data.frame"))
}

## Pre-defined index of background performances, built by the following code:
##   source( "../../R/resmine.R" )
##   synResults( "score", Region != "cerebellum" ) %>% filter( grepl( "background", name ) ) %>%
##     select( -name, -parentName ) %>% dput()
indexBackground <- function()
{
    structure(list(Dataset = c("ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "MAYOpc", "MAYOpc", "MAYOpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc",
                               "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TempCortex", "TempCortex",
                              "TempCortex", "BM10", "BM36", "BM44", "BM22", "BM36", 
                              "BM10", "BM22", "BM10", "BM44", "BM36", "BM44", "BM22"),
                   Task = c("AB", "AC", "BC", "BC", "AB", "AC", "AC", "AB", "BC",
                            "BC", "AC", "AB", "AB", "BC", "AC", "BC", "AB", "AC"),
                   id = c("syn15589822", "syn15589816", "syn15589810", "syn15572112",
                          "syn15570670", "syn15572002", "syn15584696", "syn15585576",
                          "syn15582756", "syn15585289", "syn15581358", "syn15582188", 
                          "syn15586016", "syn15585432", "syn15583710", "syn15584001",
                          "syn15583564", "syn15584563")),
              row.names = c(NA, -18L), class = c("tbl_df", "tbl", "data.frame"))
}

## Fetches all background matrices provided by the given index
bkFromIndex <- function( IDX = indexBackground() )
{
    IDX %>% mutate( BK = map(id, syn_csv), id=NULL ) %>%
        mutate_at( "BK", map, ~gather(.x, Size, AUC) ) %>%
        mutate_at( "BK", map, ~mutate_at( .x, "Size", as.numeric ) ) %>%
        mutate_at( "Dataset", str_sub, 1, -3 ) %>%
        mutate_at( "Region", recode, TempCortex = "TCX" )
}

## Annotates a results data frame with drug / target information
annotateResults <- function( R )
{
    ## Load the associated LINCS metadata (drug names)
    M <- syn_csv( "syn11801537" ) %>% mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, Drug = name, Target = target_name, URL=link )
    R %>% rename( URL = id ) %>% inner_join( M, by="URL" ) %>%
        select( LINCSID, Drug, Target, Size = intersect, AUC, p_value )
}

## Fetches all results matrices associated with a given results index
##   Annotates results with relevant tags
## IDX - index data frame, as returned by indexDGE() or indexMined()
resFromIndex <- function( IDX )
{
    ## Fetch all the relevant results values
    ## Generate a tag for each Dataset / Region / Task triplet
    csel <- c("URL" = "id","Size" = "intersect", "AUC", "p_value")
    IDX %>% mutate( Results = map(id, syn_csv) ) %>% select( -id ) %>%
        mutate_at( "Results", map, annotateResults ) %>%
        mutate_at( "Dataset", str_sub, 1, -3 ) %>%
        mutate_at( "Region", recode, TempCortex = "TCX" )
}
