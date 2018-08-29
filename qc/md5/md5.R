## Verification of consistency and discoverability of settings_md5 mappings
##
## by Artem Sokolov

source( "../../R/resmine.R" )

## Retrieves parent ids for a given set of Synapse entities
## Works on vectors of synapse ids
synParentId <- function( ids )
{
    ## Isolate the unique set of ids and retrieve the name for each
    idMap <- unique(ids) %>% purrr::set_names() %>%
        map( ~synGet( .x, downloadFile=FALSE )$properties$parentId )

    ## Extend the mapping to all the requested values
    unlist( idMap )[ids]
}

## Identifies all relevant entities for a given dataset
findEntities <- function( ds )
{
    allSettings( ds ) %>% mutate( Files = map(settings_md5, synBySettingsMd5) ) %>%
        select( settings_md5, Region, Method, Task, Files )    
}

## Count the number and type of files associated with the corresponding settings files
countFiles <- function( X )
{
    unnest(X) %>% group_by( settings_md5, parentName ) %>% summarize( nn=n() ) %>%
        spread( parentName, nn ) %>% ungroup
}

## Identifies "grandparent" folder IDs and names for a given entity matrix
findGParents <- function( X, sExclude )
{
    unnest(X) %>% filter( parentName != sExclude ) %>%
        mutate( gparntId = synParentId(parentId) ) %>%
        mutate( gparntName = synName(gparntId) )
}

## ROSMAP (protein-coding regions)
verifyROSMAP <- function()
{
    ## Identify all the relevant entities
    X <- findEntities("ROSMAP")
    countFiles(X)

    ## Ensure that settings_md5 annotation matches the "grandparent" folder name
    Y <- findGParents(X, "ROSMAPpc")
    all(Y$settings_md5 == Y$gparntName)
}

## Mayo (protein-coding regions)
verifyMayo <- function()
{
    ## Identify all the relevant entities
    X <- findEntities("Mayo")
    countFiles(X)
    
    ## Ensure that settings_md5 annotation matches the "grandparent" folder name
    Y <- findGParents(X, "MAYOpc")
    all(Y$settings_md5 == Y$gparntName)
}

## MSBB (protein-coding regions)
verifyMSBB <- function()
{
    ## Identify all the relevant entities
    X <- findEntities("MSBB")
    countFiles(X)

    ## Examine the problematic slice of data
    ## Conclusion: Missing a score and a stats file (easy to generate if needed)
    ## UPDATE: Removed from consideration, because the focus is on linear models only
    ## Z <- filter( X, settings_md5 == "74d7759158dd5393bfaa4ea1fa24fd2b" )$Files[[1]]
    
    ## Ensure that settings_md5 annotation matches the "grandparent" folder name
    Y <- findGParents(X, "MSBBpc")
    all(Y$settings_md5 == Y$gparntName)
}

## ROSMAP (old rosmap-wrangled.tsv.gz)
verifyROSMAPold <- function()
{
    X <- findEntities("syn15660477")
    countFiles(X)
    Y <- findGParents(X, "ROSMAPall")
    all(Y$settings_md5 == Y$gparntName)
}
