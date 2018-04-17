## Mining results files from Synapse
##
## by Artem Sokolov

library( magrittr )
library( tidyverse )
library( synapseClient )

synapseLogin()

## Hard-coded download location for Synapse files used by the code in this file
dlc <- "~/data/AMP-AD"

## Composes a dataset descriptor from dataset_name and processing_scheme_subset annotations
getDatasetDesc <- function( dataset_name, processing_scheme_subset )
{ ifelse( dataset_name=="MSBB_RNAseq", str_c("MSBB_",processing_scheme_subset), dataset_name ) }

## Helper function that applies getDatasetDesc() to synapse objects
## s - synapse object, as loaded by synGet(); annotations(s) must exist
synDatasetDesc <- function( s )
{ getDatasetDesc( annotations(s)$dataset_name, annotations(s)$processing_scheme_subset ) }

## Composes a task name from estimator_name and dataset_filter_name annotations
getTaskDesc <- function( estimator_name, dataset_filter_name )
{ ifelse( estimator_name=="Ordinal", estimator_name, dataset_filter_name ) }

## Helper function that applies getTaskDesc() to synapse objects
## s - synapse object, as loaded by synGet(); annotations(s) must exist
synTaskDesc <- function( s )
{ getTaskDesc( annotations(s)$estimator_name, annotations(s)$dataset_filter_name ) }

## Downloads a mapping of LINCS URLs to drug names
## Returns a mapping that can be fed directly to synStats()
LINCSmap <- function()
{
    s <- synGet( "syn11801537", downloadLocation = dlc )
    M <- read_csv( s@filePath, col_types = cols() )
    setNames( M$name, M$link )
}

## Retrieves a summary table of all *stats files
statsSummary <- function()
{
    ## Identify all the stats files, their corresponding Synapse IDs and
    ##   basic annotations of interest
    aa <- "id,name,dataset_name,processing_scheme_subset,estimator_name,dataset_filter_name"
    str_c( "select ", aa, " from file where btr_file_type==\"stats\"" ) %>% synQuery %>%
        rename_all( ~str_split( ., "\\.", simplify=TRUE )[,2] ) %>%
        mutate( Dataset = getDatasetDesc( dataset_name, processing_scheme_subset ),
               Task = getTaskDesc( estimator_name, dataset_filter_name ) ) %>%
        select( SynID = id, Filename = name, Dataset, Task ) %>% as_data_frame
}

## Loads stats from a given synapseID
## nameMap - optional mapping of names to another space; must be an indexable object
synStats <- function( synid, nameMap=NULL, verbose=FALSE )
{
    if( verbose ) cat( "Loading", synid, "...\n" )
    s <- synGet( synid, downloadLocation = dlc )
    read_csv( s@filePath, col_types = cols() ) %>%
        mutate( Dataset = synDatasetDesc(s), Task = synTaskDesc(s) ) %>%
        select( Name = id, Size = intersect, AUC, pval = p_value, Dataset, Task ) %>%
        mutate( Name = if(is.null(nameMap)) Name else nameMap[Name] )
}

## Loads all the stats files whose filenames grep-match the requested pattern
synStatsMany <- function( fnPattern, nameMap=NULL, verbose=TRUE )
{
    statsSummary() %>% filter( grepl(fnPattern, Filename) ) %$%
        setNames( SynID, Filename ) %>% map( ~synStats(.x, nameMap, verbose) ) %>%
        bind_rows( .id = "Filename" )
}

## Retrieves a summary table of all background score files
bkSummary <- function()
{
    ## Identify all the background score files
    aa <- "id,name,dataset_name,processing_scheme_subset,estimator_name,dataset_filter_name"
    str_c( "select ", aa, " from file where btr_file_type==\"score\"" ) %>% synQuery %>%
        rename_all( ~str_split( ., "\\.", simplify=TRUE )[,2] ) %>%
        mutate( Dataset = getDatasetDesc( dataset_name, processing_scheme_subset ),
               Task = getTaskDesc( estimator_name, dataset_filter_name ) ) %>%
        filter( name == "background_auc.csv" ) %>%
        select( SynID = id, Dataset, Task ) %>% as_data_frame
}

## Loads a background_auc.csv file
synBkAUC <- function( synid, verbose = FALSE )
{
    if( verbose ) cat( "Loading", synid, "...\n" )
    s <- synGet( synid, downloadLocation = dlc )
    read_csv( s@filePath, col_types = cols() ) %>% gather( Size, AUC ) %>%
        mutate( Size = as.integer(Size), Dataset = synDatasetDesc(s), Task = synTaskDesc(s) )
}

## Loads all (non-synthetic) background AUC files
synBkAll <- function( verbose = TRUE )
{
    bkSummary() %>% filter( !grepl("synthetic", Dataset) ) %>% extract2( "SynID" ) %>%
        map( ~synBkAUC(.x, verbose) ) %>% bind_rows
}

