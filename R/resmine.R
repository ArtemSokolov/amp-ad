## API for mining the latest results based on pair selection that minimizes AOD
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

## Lists all settings files available in a given Synapse directory
## parentId may be one of "mayo", "rosmap", "msbb", or a synapse ID
listSettings <- function( parentId )
{
    ## Map parentId to the appropriate synapse ID if provided as a keyword
    vMap <- c( "mayo" = "syn12180241", "rosmap" = "syn15589860", "msbb" = "syn15588043" )
    if( str_to_lower(parentId) %in% names(vMap) ) parentId <- vMap[str_to_lower(parentId)]
    
    synExtra::synq( "id","name","settings_md5", parentId = parentId, Type = "Settings" ) %>% unnest
}

## A clean version of listSettings() that parses the file name into chunks
allSettings <- function( parentId )
{
    listSettings( parentId ) %>% filter(grepl( "^settings", name )) %>%
        mutate( chunks = str_split(str_sub(name, 1, -6), "\\_") ) %>%
        mutate( Dataset = map_chr(chunks, nth, 2),
               Region = map_chr(chunks, nth, 3),
               Method = map_chr(chunks, nth, 4),
               Strategy = map_chr(chunks, nth, 5),
               Task = map_chr(chunks,6) ) %>%
        select( -name, -chunks )
}

## Lists all files associated with a given settings md5 hash
synBySettingsMd5 <- function( md5 )
{
    synExtra::synq( "id","name","parentId", settings_md5 = md5 ) %>%
        mutate( parentName = synExtra::synName(parentId) )
}

