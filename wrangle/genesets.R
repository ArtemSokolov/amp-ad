## Wrangling gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD", ifcollision = "overwrite.local" )$path }

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
  read_lines(fn) %>% str_split( "\\t" ) %>%
    set_names( map_chr(., nth, iName) ) %>%
    map( ~.x[-2:-1] )
}

## Extracts a gene set associated with dsRNA proteomics dataset
geneset_dsRNAprot <- function()
{
    ## Retrieve the raw data
    X <- syn( "syn11807753" ) %>% read_csv( col_types = cols() ) %>%
        select( Gene=`Gene Symbol`, Control1=`Control Replicate 1`, Control2=`Control Replicate 2`,
               dsRNAmi1 = `DsRNAmi Replicate 1`, dsRNAmi2 = `DsRNAmi Replicate 2` ) %>%
        mutate( Control = map2_dbl(Control1, Control2, ~mean(c(.x,.y),na.rm=TRUE)),
               dsRNAmi = map2_dbl(dsRNAmi1, dsRNAmi2, ~mean(c(.x,.y),na.rm=TRUE)) ) %>%
        mutate( FoldChange = dsRNAmi / Control )

    ## Corner case: protein not expressed in control leading infinite fold change
    v1 <- X %>% filter( Control == 0, dsRNAmi > 10 ) %>% .$Gene

    ## Consider remaining data and compute log2 fold change
    XFC <- X %>% filter( Control != 0 ) %>% mutate( lFC = log2(FoldChange) )

    ## Write the result to file and store it on Synapse
    vup <- XFC %>% filter( lFC > 1 ) %>% .$Gene
    vdn <- XFC %>% filter( lFC < -1 ) %>% .$Gene

    ## Write the results to their respective files and store them to Synapse
    c( v1, vup ) %>% unique %>% cat( file="ReNprot-up.txt", sep="\n" )
    unique( vdn ) %>% cat( file="ReNprot-dn.txt", sep="\n" )
    
    ## Combine the two lists and retrieve unique entries
    ## Write the result to a file and store the file to Synapse
    aa <- list( Type="Gene Set", Category="Interferome", Reference="LSP Experiment" )
    f1 <- File( "ReNprot-up.txt", parentId = "syn11629934" ) %>% synStore() %>% synSetAnnotations(aa)
    f2 <- File( "ReNprot-dn.txt", parentId = "syn11629934" ) %>% synStore() %>% synSetAnnotations(aa)
}

## "Quiet" versions of read_tsv and read_csv
readq_tsv <- function(...) read_tsv( col_types = cols(), ... )
readq_csv <- function(...) read_csv( col_types = cols(), ... )

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

## Composes gene sets associated with DGE data
## These sets are size-matched against `Nienke_10genes.gmt`
##   in the context of specific datasets
genesets_DGE <- function()
{
    ## Load gene space of each dataset
    AA <- syn( "syn16204450" ) %>% read_gmt()
    
    ## Load Nienke's mined gene sets and the corresponding metadata
    ##  Count the number of mined associations present
    MGS <- syn( "syn11973633" ) %>% read_gmt() %>% enframe( "URL", "MinedSet" )
    NK <- syn( "syn11801537" ) %>% readq_csv() %>% rename( URL = link ) %>%
        mutate( Drug = str_to_lower(name) ) %>% inner_join( MGS, by="URL" ) %>%
        mutate( ROSMAP = map_int( MinedSet, ~length(intersect(., AA$ROSMAP)) ),
               Mayo = map_int( MinedSet, ~length(intersect(., AA$Mayo)) ),
               MSBB = map_int( MinedSet, ~length(intersect(., AA$MSBB)) ) ) %>%
        select( LINCSID = lincs_id, URL, Drug, ROSMAP:MSBB )

    ## Load differential gene expression scores for all drugs
    DFE <- syn( "syn15674107" ) %>% read_csv( col_types = cols() ) %>%
        mutate_at( "Drug", str_to_lower ) %>% filter( Drug %in% NK$Drug )

    ## Ensure the genes are sorted by p-value and encapsulate them inside a nested frame
    GG <- DFE %>% group_by(Drug) %>% arrange(PValue) %>% ungroup() %>% select(Drug, Gene) %>%
        nest( Gene ) %>% mutate( Genes = map(data, ~.$Gene) ) %>% select( -data )

    ## NK contains the number of genes Nienke's sets have in common with each dataset
    ##  Use these values to select the top k differentially-expressed genes accordingly
    RR <- inner_join( GG, NK, by="Drug" ) %>% mutate( HDR = str_c( LINCSID, "\t", URL ) ) %>%
        mutate_at( vars(ROSMAP:MSBB), map2, .$Genes, ~.y[1:.x] ) %>%
        select( HDR, ROSMAP, Mayo, MSBB ) %>%
        mutate_at( vars(ROSMAP:MSBB), map_chr, str_flatten, "\t" ) %>%
        transmute_at( vars(ROSMAP:MSBB), map2_chr, .$HDR, ~str_c(.y, "\t", .x) )

    ## Save everything to .gmt files
    cat( RR$ROSMAP, sep="\n", file="DGE-ROSMAP.gmt" )
    cat( RR$Mayo, sep="\n", file="DGE-Mayo.gmt" )
    cat( RR$MSBB, sep="\n", file="DGE-MSBB.gmt" )
}
