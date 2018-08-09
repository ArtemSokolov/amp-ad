## API for mining the latest results based on pair selection that minimizes AOD
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
library( assertr )

synapseLogin()

## Retrieves file names associated with a synapse IDs
## Works on vectors of ids
synName <- function( ids )
{
    ## Isolate the unique set of ids and retrieve the name for each
    idMap <- unique(ids) %>% purrr::set_names() %>%
        map( ~synGet( .x, downloadFile=FALSE )@properties$name )

    ## Extend the mapping to all the requested values
    unlist( idMap )[ids]
}

## Robust version of synQuery
## Returns a data frame with 0 rows, instead of NULL, if there are no matches
synq <- function( what, ... )
{
    ## Assemble the field request
    fields <- str_flatten(what, ",")

    ## Assemble the requested conditions
    cond <- enquos(...) %>% map( ~rlang::eval_tidy(.x) ) %>%
        imap( ~str_c('"', .y, '"=="', .x, '"') ) %>% str_flatten(" and ")

    ## Compose the query and pass it to Synapse API
    qq <- str_c( "select ", fields, " from file where ", cond )
    cat( qq, "\n" )
    QQ <- synQuery(qq)
    if( is.null(QQ) )
        purrr::set_names(what) %>% map_dfc( ~character() )
    else
        as_data_frame(QQ) %>% rename_all( ~str_split( ., "\\.", simplify=TRUE )[,2] )
}

## Lists all settings files available in a given Synapse directory
## parentId may be one of "mayo", "rosmap", "msbb", or a synapse ID
listSettings <- function( parentId )
{
    ## Map parentId to the appropriate synapse ID if provided as a keyword
    vMap <- c( "mayo" = "syn12180241", "rosmap" = "syn15589860", "msbb" = "syn15588043" )
    if( str_to_lower(parentId) %in% names(vMap) ) parentId <- vMap[str_to_lower(parentId)]
    
    synq( c("id","name","settings_md5"), parentId = parentId, type = "settings" )
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
        select( -id, -name, -chunks )
}

## Lists all files associated with a given settings md5 hash
synBySettingsMd5 <- function( md5 )
{
    synq( c("id","name","parentId"), settings_md5 = md5 ) %>%
        mutate( parentName = synName(parentId) )
}

