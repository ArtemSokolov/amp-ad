## Gene set composition for DGE experiments
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/AMP-AD" )

## "Quiet" versions of read_tsv and read_csv
readq_tsv <- function(...) read_tsv( col_types = cols(), ... )
readq_csv <- function(...) read_csv( col_types = cols(), ... )

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
  read_lines(fn) %>% str_split( "\\t" ) %>%
    set_names( map_chr(., nth, iName) ) %>%
    map( ~.x[-2:-1] )
}

## Retrieves the set of all gene names in protein-coding versions of ROSMAP, Mayo and MSBB
genesets_AMPAD <- function()
{
    ## Load the three datasets and extract gene names from each
    vExclude <- c("ID","PMI","AOD","CDR","Braak","BrodmannArea","Barcode",
                  "Diagnosis","Region","Thal","ApoE","Gender","Source","TDP-43","LBD")
    synIDs <- c( ROSMAP = "syn14306482", Mayo = "syn15059396", MSBB = "syn15094062" )
    GG <- map( synIDs, ~readq_tsv(syn(.x)) ) %>% map( colnames ) %>% map( setdiff, vExclude )

    ## Put everything together
    RR <- map_chr( GG, str_flatten, "\t" ) %>% enframe( "Name", "Set" ) %>%
        inner_join( enframe( synIDs, "Name", "Desc" ), ., by="Name" ) %>%
        mutate( Final = str_c(Name, Desc, Set, sep="\t") )

    ## Compose the final strings and write them to file
    cat( RR$Final, sep="\n", file="AMP-AD.gmt" )
}

## Returns the number of genes that Nienke's sets have in common with each AMP-AD dataset
NKsizes <- function()
{
    ## Load gene space of each dataset
    AA <- syn( "syn16204450" ) %>% read_gmt()
    
    ## Load Nienke's mined gene sets and the corresponding metadata
    ##  Count the number of mined associations present
    MGS <- syn( "syn11973633" ) %>% read_gmt() %>% enframe( "URL", "MinedSet" )
    syn( "syn11801537" ) %>% readq_csv() %>% rename( URL = link ) %>%
        mutate( Drug = str_to_lower(name) ) %>% inner_join( MGS, by="URL" ) %>%
        mutate( ROSMAP = map_int( MinedSet, ~length(intersect(., AA$ROSMAP)) ),
               Mayo = map_int( MinedSet, ~length(intersect(., AA$Mayo)) ),
               MSBB = map_int( MinedSet, ~length(intersect(., AA$MSBB)) ) ) %>%
        select( LINCSID = lincs_id, URL, Drug, ROSMAP:MSBB )
}

## Given a differential expression matrix, composes gene sets for each drug that are
##  size-matched against `Nienke_10genes.gmt` in the context of specific datasets
genesets_smatch <- function( DX )
{
    ## Retrieve the number of genes Nienke's sets have in common with each dataset
    NK <- NKsizes()
    vDS <- c("ROSMAP","Mayo","MSBB")	## Dataset-specific column names

    ## Sort drug-specific differential expression by p-value and isolate the top k genes
    DX %>% mutate_at( "Drug", str_to_lower ) %>% nest( -Drug, .key="DFX" ) %>%
        mutate_at( "DFX", map, arrange, sym("PValue") ) %>%
        inner_join( NK, by="Drug" ) %>% mutate_at( vDS, map2, .$DFX, ~slice(.y, 1:.x) ) %>%
        mutate_at( vDS, map, pull, "Gene" ) %>% select(-DFX)
}

## Takes the output of genesets_smatch() and concatenates it into tab-delimited
##   strings in preparation for writing to .gmt files
prep4gmt <- function( GS )
{
    GS %>% mutate( HDR = str_c( LINCSID, "\t", URL ) ) %>%
        mutate_at( vars(ROSMAP,Mayo,MSBB), map_chr, str_flatten, "\t" ) %>%
        transmute_at( vars(ROSMAP,Mayo,MSBB), map2_chr, .$HDR, ~str_c(.y, "\t", .x) )
}

## Processes DGE1 experiment
DGE1 <- function()
{
    ## Compute size-matched gene sets and concatenate into tab-delimited strings
    RR <- syn( "syn15674107" ) %>% readq_csv() %>% genesets_smatch() %>% prep4gmt()

    ## Save everything to .gmt files
    cat( RR$ROSMAP, sep="\n", file="DGE1-ROSMAP-sm.gmt" )
    cat( RR$Mayo, sep="\n", file="DGE1-MAYO-sm.gmt" )
    cat( RR$MSBB, sep="\n", file="DGE1-MSBB-sm.gmt" )
}

## Processes DGE2 experiment
DGE2 <- function()
{
    ## Compute size-matched gene sets and concatenate into tab-delimited strings
    RR <- syn( "syn17167348" ) %>% readq_csv() %>% genesets_smatch() %>% prep4gmt()
    
    ## Save everything to .gmt files
    cat( RR$ROSMAP, sep="\n", file="DGE2-ROSMAP-sm.gmt" )
    cat( RR$Mayo, sep="\n", file="DGE2-MAYO-sm.gmt" )
    cat( RR$MSBB, sep="\n", file="DGE2-MSBB-sm.gmt" )
}

