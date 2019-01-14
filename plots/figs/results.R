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
##   synResults() %>% filter( Type == "stats", Region != "CBE" ) %>% unnest() %>%
##     filter( grepl("DGE", fileName) ) %>% mutate( Plate = str_sub(fileName, 1, 4) ) %>%
##     select( id = fileId, -fileName ) %>% dput()
indexDGE <- function()
{
    structure(list(Dataset = c("MAYO", "MAYO", "MAYO", "MAYO", "MAYO", 
                               "MAYO", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB" ),
                   Region = c("TCX", "TCX", "TCX", "TCX", "TCX", "TCX", "DLPFC", 
                              "DLPFC", "DLPFC", "DLPFC", "DLPFC", "DLPFC", "BM10", "BM10", 
                              "BM10", "BM10", "BM10", "BM10", "BM22", "BM22", "BM22", "BM22", 
                              "BM22", "BM22", "BM36", "BM36", "BM36", "BM36", "BM36", "BM36", 
                              "BM44", "BM44", "BM44", "BM44", "BM44", "BM44"),
                   Task = c("AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", "BC", "BC", 
                            "AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", "BC", 
                            "BC", "AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", 
                            "BC", "BC"),
                   Type = c("stats", "stats", "stats", "stats", "stats", 
                            "stats", "stats", "stats", "stats", "stats", "stats", "stats", 
                            "stats", "stats", "stats", "stats", "stats", "stats", "stats", 
                            "stats", "stats", "stats", "stats", "stats", "stats", "stats", 
                            "stats", "stats", "stats", "stats", "stats", "stats", "stats", 
                            "stats", "stats", "stats"),
                   id = c("syn17434757", "syn17434764", 
                          "syn17434797", "syn17434803", "syn17434524", "syn17434535", "syn17434578", 
                          "syn17434582", "syn17435113", "syn17435116", "syn17434887", "syn17434893", 
                          "syn17435027", "syn17435035", "syn17435069", "syn17435076", "syn17434396", 
                          "syn17434399", "syn17434620", "syn17434623", "syn17434345", "syn17434352", 
                          "syn17434709", "syn17434712", "syn17435234", "syn17435244", "syn17435191", 
                          "syn17435197", "syn17434485", "syn17434488", "syn17434435", "syn17434445", 
                          "syn17434836", "syn17434848", "syn17434664", "syn17434672"), 
                   Plate = c("DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", 
                             "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", 
                             "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", 
                             "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", 
                             "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2")),
              class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -36L))
}

## Pre-defined index of results for Nienke's gene sets, built by the following code:
##   source( "../../R/resmine.R" )
##   synResults( "stats", Region != "cerebellum" ) %>% filter( grepl("Nienke", name) ) %>%
##     select( -name, -parentName ) %>% dput()
## Dataset and Region names are then adjusted by hand:
## 1. Removal of pc suffix from dataset names
## 2. Recoding TempCortex as TCX
indexMined <- function()
{
    structure(list(Dataset = c("ROSMAP", "ROSMAP", "ROSMAP", "MAYO", "MAYO", "MAYO",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TCX", "TCX",
                              "TCX", "BM10", "BM36", "BM44", "BM22", "BM36", 
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
## Dataset and Region names are then adjusted by hand:
## 1. Removal of pc suffix from dataset names
## 2. Recoding TempCortex as TCX
indexBackground <- function()
{
    structure(list(Dataset = c("ROSMAP", "ROSMAP", "ROSMAP", "MAYO", "MAYO", "MAYO",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TCX", "TCX",
                              "TCX", "BM10", "BM36", "BM44", "BM22", "BM36", 
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
        mutate_at( "BK", map, ~mutate_at( .x, "Size", as.numeric ) )
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
    IDX %>% mutate( Results = map(id, syn_csv) ) %>% select( -id ) %>%
        mutate_at( "Results", map, annotateResults )
}
