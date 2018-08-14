library( tidyverse )
library( synapser )
synLogin()

## Updates annotation settings_md5 of a given entity with the provided hash
updateAnnot <- function( synID, md5 )
{
    cat( "Updating entity", synID, "\n" )
    v <- synGetAnnotations( synID ) %>% flatten
    v$settings_md5 <- md5
    synSetAnnotations( synID, v )
}

main <- function()
{
    X <- read_csv( "fix1.csv" )
    Y <- X %>% mutate( Annots = map2( id, gparntName, updateAnnot ) )
}
