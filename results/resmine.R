## API for mining the latest results based on pair selection that minimizes AOD
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
library( assertr )

synapseLogin()

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD" )
{ synGet( id, downloadLocation = dlc )@filePath }

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
listSettings <- function( parentId = "syn12180241" )
{
    synq( c("id","name","settings_md5"), parentId = parentId, type = "settings" )
}

## Lists all files associated with a given settings md5 hash
synBySettingsMd5 <- function( md5 )
{
    synq( c("id","name","parentId"), settings_md5 = md5 ) %>%
        mutate( parentName = synName(parentId) )
}

## Identifies all stats files associated with a given strategy index,
##  task (AB, BC, AC) and linear / nonlinear distinction
listStats <- function( task, stratIndex = 3, lnl = "Linear" )
{
    ## The query should resolve to a unique entry
    ## Retrieve the entry's md5 hash, then compose another query to retrieve
    ##   all files associated with the hash
    ## Finally, reduce the list of files just to those annotated as "stats"
    s <- str_c( "strategy", stratIndex )
    listSettings() %>% filter( grepl(s, name) ) %>%
        filter( grepl(task, name) ) %>% filter( grepl(lnl, name) ) %>%
        verify( nrow(.) == 1 ) %>% .$settings_md5 %>% synBySettingsMd5 %>%
        filter( parentName == "stats" )
}
