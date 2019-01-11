## Uploading results to their correct location on Synapse
##   Marks each file with their corresponding annotations
##
## by Artem Sokolov

library( tidyverse )
library( synapser )
library( synExtra )

synLogin()

## Given a settings_md5 key, finds its Runs folder and returns the corresponding synapse ID
findBySMD5 <- function( smd5 )
{
    c("syn12180240", "syn15581352", "syn15589794" ) %>% map(synChildren) %>%
        lift(c)() %>% pluck( smd5 )
}


## Annotates and uploads a file to its correct directory on Synapse
resUpload <- function( fnRes )
{
    cat( "Annotating and uploading", fnRes, "\n" )
    
    ## Retrieve the file type (score or stats) and the associated md5 key
    fnRes <- normalizePath( fnRes )
    ft <- dirname(fnRes) %>% basename
    md5 <- dirname(fnRes) %>% dirname %>% basename

    ## Identify and parse the associated annotations file
    f0 <- tools::file_path_sans_ext( fnRes )
    fnAnnot <- str_c( dirname(f0), "/.annotations/", basename(f0), ".json" )
    stopifnot( file.exists( fnAnnot ) )
    A <- jsonlite::fromJSON( fnAnnot )

    ## Handle square brackets in annotation names
    names(A) <- names(A) %>% gsub( "\\[", "LFTB", . ) %>% gsub( "\\]", "RGTB", . )

    ## Identify the destination directory on synapse
    synTgt <- findBySMD5(md5) %>% synPluck( ft )

    ## Create a File entity and synStore() it to the identified directory
    File( fnRes, parent = synTgt, annotations = A ) %>% synStore()
}

resUploadAll <- function()
{
    RR <- list.files( pattern = "DGE", recursive=TRUE ) %>%    ## Identify all DGE-related files
        keep( !grepl("hypotheses",.) ) %>%                     ## Exclude .gmt gene set files
        map(resUpload)                                         ## Store each to Synapse
}
