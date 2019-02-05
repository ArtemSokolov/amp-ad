## Gene set composition for DGE experiments
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/AMP-AD", ifcollision="overwrite.local" )

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

## Given a differential expression matrix, composes gene sets for each drug that
##   consist of [significantly] differentially-expressed genes
## Caps the number of genes at 300 (approximate size of the largest mined set)
genesets_dfexp <- function( DX )
{
    ## Given a differential expression matrix (for a single drug), identifies and returns
    ## genes with FDR < thresh. If there are fewer than atLeast such genes, returns
    ## genes with PValue < thresh instead
    top_genes <- function( .df, thresh = 0.05, atLeast=10 )
    {
        vFDR  <- with( .df, set_names(FDR, Gene) ) %>% keep( ~ .x < thresh ) %>% names
        vpval <- with( .df, set_names(PValue, Gene) ) %>% keep( ~ .x < thresh ) %>% names
        if( length(vFDR) < atLeast ) return(vpval)
        vFDR
    }

    vDS <- c("ROSMAP","Mayo","MSBB")	## Dataset-specific column names
    
    ## Load gene space of each dataset
    AA <- syn( "syn16204450" ) %>% read_gmt()

    ## Subset the gene space to genes in common with each dataset
    ee <- set_names(vDS) %>% map( ~expr( map(Raw, filter, Gene %in% AA[[!!.x]]) ) )
    GS <- DX %>% mutate_at( "Drug", str_to_lower ) %>% nest( -Drug, .key="Raw" ) %>%
        mutate( !!!ee, Raw=NULL ) %>% mutate_at( vDS, map, top_genes ) %>%
        mutate_at( vDS, map_if, ~(length(.x) > 300), ~.x[1:300] )

    ## Annotate with metadata
    M <- syn("syn11801537") %>% read_csv() %>% mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, URL = link, Drug = name )
    stopifnot( all(GS$Drug %in% M$Drug) )

    inner_join( M, GS )
}

## Takes the output of genesets_smatch() and concatenates it into tab-delimited
##   strings in preparation for writing to .gmt files
prep4gmt <- function( GS )
{
    GS %>% mutate( HDR = str_c( LINCSID, "\t", URL ) ) %>%
        mutate_at( vars(ROSMAP,Mayo,MSBB), map_chr, str_flatten, "\t" ) %>%
        transmute_at( vars(ROSMAP,Mayo,MSBB), map2_chr, .$HDR, ~str_c(.y, "\t", .x) ) %>%
        rename( MAYO = Mayo )
}

## Processes DGE1 experiment
DGE1 <- function()
{
    vDS <- c("ROSMAP", "MAYO", "MSBB")
    
    ## Load the raw differential expression matrix
    DX <- syn( "syn15674107" ) %>% readq_csv()
    
    ## Compute size-matched gene sets and concatenate into tab-delimited strings
    GSsm <- DX %>% genesets_smatch() %>% prep4gmt()
    map( vDS, ~cat(GSsm[[.x]], sep="\n", file=str_c("DGE1-", .x, "-sm.gmt")) )

    ## Compute dfx-based gene sets and concatenate into tab-delimited strings
    GSdfx <- DX %>% genesets_dfexp() %>% prep4gmt()
    map( vDS, ~cat(GSdfx[[.x]], sep="\n", file=str_c("DGE1-", .x, "-dfx.gmt")) )
}

## Processes DGE2 experiment
DGE2 <- function()
{
    vDS <- c("ROSMAP", "MAYO", "MSBB")
    
    ## Load the raw differential expression matrix
    DX <- syn( "syn17167348" ) %>% readq_csv()
    
    ## Compute size-matched gene sets and concatenate into tab-delimited strings
    GSsm <- DX %>% genesets_smatch() %>% prep4gmt()
    map( vDS, ~cat(GSsm[[.x]], sep="\n", file=str_c("DGE2-", .x, "-sm.gmt")) )

    ## Compute dfx-based gene sets and concatenate into tab-delimited strings
    GSdfx <- DX %>% genesets_dfexp() %>% prep4gmt()
    map( vDS, ~cat(GSdfx[[.x]], sep="\n", file=str_c("DGE2-", .x, "-dfx.gmt")) )
}

