## Interface with Synapse to pull/push files necessary for
##   btr-score and btr-stats scripts
##
## by Artem Sokolov

library( tidyverse )
library( synapser )
library( synExtra )

synLogin()

## Given a synapse ID, retrieves settings_md5 annotation associated with it
getSettingsMd5 <- function( synid )
{
    s <- synGet( synid, downloadFile=FALSE )
    s$annotations$get("settings_md5")
}

## Given a settings_md5 key, finds its Runs folder and returns the corresponding synapse ID
findBySMD5 <- function( smd5 )
{
    c("syn12180240", "syn15581352", "syn15589794" ) %>% map(synChildren) %>%
        lift(c)() %>% pluck( smd5 )
}

## For a given settings_md5, downloads its background_auc.csv to the
##   local Runs/<md5 key>/score directory
getBKAUC <- function( md5 )
{
    ## Ensure the key is present in the local Runs directory
    stopifnot( md5 %in% list.files("Runs") )

    ## Identify the file to download
    sid <- findBySMD5(md5) %>% synPluck( "score", "background_auc.csv" )
    
    ## Retrieve the file using synGet()
    ff <- synGet( sid, downloadLocation = str_c( "Runs/", md5, "/score" ),
                 ifcollision = "overwrite.local" )

    ## Ensure that the output file is named background_auc.csv
    fnOut <- str_c( dirname( ff$path ), "/background_auc.csv" )
    file.rename( ff$path, fnOut )
    fnOut
}

getBKAUC_all <- function()
{
    ## Retrieve all settings md5 hashes
    S <- c( "mayo" = "syn12180241", "rosmap" = "syn15589860", "msbb" = "syn15588043" ) %>%
        map( synChildren ) %>% map( enframe, "name", "id" ) %>%
        bind_rows() %>% mutate( smd5 = map_chr(id, getSettingsMd5) )

    ## Consider only those present in the local Runs directory
    ## Download all the background_auc.csv files
    S %>% filter( smd5 %in% list.files("Runs") ) %>% pull( smd5 ) %>% map( getBKAUC )
}

