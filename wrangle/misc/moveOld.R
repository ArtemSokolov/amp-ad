## Re-arranging old files on Synapse

library( tidyverse )

library( synapser )
library( synExtra )
synLogin()

## Given a run ID, identifies all DGE-related files in it
findDGE <- function( runID )
{
    cat( "Retrieving DGE-related files for", runID, "\n" )
    data_frame( resName = c("hypothesis_predictions", "score", "stats") ) %>%
        mutate( resID = map_chr(resName, ~synPluck(runID, .x)), ch = synChildren(resID) ) %>%
        mutate_at( "ch", map, enframe, "name", "id" ) %>% unnest() %>%
        filter( grepl("DGE", name) )
}

## Moves a DGE file to the 2019-01-09 archive folder (syn17232696)
## Creates appropriate directories to avoid name collisions
moveDGE <- function( id, resName, md5sum )
{
    cat( "Moving", id, "to", md5sum, "/", resName, "\n" )
    
    ## Create the md5sum and resName directories (or access existing ones)
    f1 <- Folder( md5sum, parent = "syn17232696" ) %>% synStore()
    f2 <- Folder( resName, parent = f1$properties$id ) %>% synStore()

    ## Move the file to its destination
    synMove( id, f2$properties$id )
}

main <- function()
{
    ## Identify all current Runs directories and DGE-related files in them
    RR <- list( MAYO="syn12180240", MSBB="syn15581352", ROSMAP="syn15589794" ) %>%
        map( synChildren ) %>% map( enframe, "md5sum", "runID" ) %>% bind_rows( .id="Dataset" ) %>%
        mutate( DGEfiles = map(runID, findDGE) ) %>% unnest()

    ## Move each DGE file accordingly
    RR %>% select( id, resName, md5sum ) %>% pmap( moveDGE )
}
