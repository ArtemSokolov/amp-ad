## Retrieves the current pairs assignment for a given Dataset / Region / Task
##
## by Artem Sokolov

source( "../R/resmine.R" )

## Downloads a synapse ID to local disk and returns its filename
syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/Pairs" )@filePath }

## Given a data frame of files associated with a given settings_md5
##   retrieves synapse ID of a single background_predictions instance
getBKPid <- function( F )
{
    F %>% filter( parentName == "background_predictions" ) %>% slice(1) %>% .$id
}

## Given the synapse ID of a runs file, retrieves its pair assignment
getPairs <- function( sid )
{
    syn( sid ) %>% read_csv( col_types = cols(ID=col_character()) ) %>%
        select( ID, pair_index )    
}

## Retrieves all pair assignments associated with a given dataset
fetchPairs <- function( ds, fn )
{
    ## Retrieve all the relevant settings files
    X <- allSettings( ds ) %>% filter( Method == "sklLR" ) %>%
        select( -Strategy, -Dataset ) %>%
        mutate( Files = map(settings_md5, synBySettingsMd5) ) %>%
        mutate( bkpid = map_chr(Files, getBKPid) )

    ## Pull one a single instance of Runs associated with each settings file
    XX <- X %>% select( -Files, -settings_md5 ) %>% mutate( Pairs = map(bkpid, getPairs) )

    ## Finalize and write to file
    XX %>% select( -Method, -bkpid ) %>% unnest %>% distinct %>% write_csv( fn )

    ## Upload the file to Synapse
    f <- File( fn, parentId = "syn15590460" )
    synStore(f)
}

main <- function()
{
    fetchPairs( "ROSMAP", "ROSMAPpc-pairs.csv" )
    fetchPairs( "Mayo", "MAYOpc-pairs.csv" )
    fetchPairs( "MSBB", "MSBBpc-pairs.csv" )
}

                
