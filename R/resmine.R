## API for mining the latest results based on pair selection that minimizes AOD
##
## by Artem Sokolov

library( tidyverse )
library( synapser )
library( synExtra )

synLogin()

## Returns all available results files
synResults <- function()
{
    ## Given a synapse ID, retrieves settings_md5 annotation associated with it
    getSettingsMd5 <- function( synid )
    {
        s <- synGet( synid, downloadFile=FALSE )
        s$annotations$get("settings_md5")
    }

    ## List of valid types
    vTypes = c("background_predictions", "hypothesis_predictions", "score", "stats")
    
    cat( "Cataloguing available runs...\n" )
    R <- c( "MAYO" = "syn12180241", "ROSMAP" = "syn15589860", "MSBB" = "syn15588043" ) %>%
        map( synChildren ) %>% map( enframe, "name", "synid" ) %>% bind_rows( .id = "Dataset" ) %>%
        mutate( md5 = map_chr( synid, getSettingsMd5 ),
               chunks = str_split(str_sub(name, 1, -6), "\\_"),
               Region = map_chr(chunks, nth, 3),
               Task = map_chr(chunks,6) ) %>% select( -name, -synid, -chunks ) %>%
        mutate_at( "Region", recode, TempCortex = "TCX", cerebellum = "CBE" )

    cat( "Identifying files in each run...\n" )
    R %>% mutate( synRun = map2_chr(Dataset, md5, ~synPluck("syn15590460", str_c(.x, "pc"), .y)),
                 Type = rep( list(vTypes), n() ), md5=NULL ) %>% unnest() %>%
        mutate( synType = map2_chr(synRun, Type, synPluck),
               Files = map(synType, synChildren) ) %>%
        mutate_at( "Files", map, enframe, "fileName", "fileId" ) %>%
        select( -synRun, -synType )
}

