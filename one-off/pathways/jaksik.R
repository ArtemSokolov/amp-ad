library( tidyverse )
library( indRa )

## Load relevant Python modules
reticulate::use_virtualenv("~/.virtualenvs/indra/", TRUE)

## Identifies the set of transcription factors that capture interferome-stimulated genes
## External resources:
## ISG list: syn11629935
## TF-target map: https://github.com/bioinfonerd/Transcription-Factor-Databases/
##                                            blob/master/Ttrust_v2/TF_GSEA/all.gmt.gz
## NOTE: The TF-target map is currently NOT in proper .gmt format
##       The code will need to be updated after Nathan fixes the format
TFisg <- function()
{
    ## Setup Synapse downloaders
    synapser::synLogin()
    syn <- synExtra::synDownloader(".", ifcollision="overwrite.local")

    ## Load TF -> Target relationships
    ## Remove TFs that only have a single outgoing edge
    TF <- read_lines("all.gmt.gz") %>% str_split( "\\t" ) %>%
        set_names( map_chr(., first) ) %>% map( ~.x[-1] ) %>%
        keep( ~length(.)>1 ) %>% enframe( "Name", "Targets" )

    ## Load Interferome-Stimulated Genes
    isg <- syn( "syn11629935" ) %>% read_lines()

    ## Determine ISGs that are captured by each TF
    TF <- TF %>% mutate( ISGs = map(Targets, intersect, isg),
                        nISG = map_int(ISGs, length), nTrgt = map_int(Targets, length) ) %>%
        arrange( desc(nISG) ) %>% filter( nISG > 0 )

    ## Population size - total number of black+white balls in an urn
    ## Use total number of genes
    nPop <- 20000

    ## Sample size - number of balls drawn from an urn
    ## Use number of ISGs
    nSampl <- length(isg)

    ## Number of successes in the population - total number of white balls in an urn
    ## Given by number of targets captured by a particular TF
    ## Column nTrgt in the data frame

    ## Number of successes in the sample - number of white balls drawn
    ## Given by number of ISGs captured by a particular TF
    ## Column nISG in the data frame

    ## Perform hypergeometric test for each TF
    fHGT <- ~phyper(.x-1, .y, nPop-.y, nSampl, lower.tail=FALSE)
    TF %>% mutate( pval= map2_dbl(nISG, nTrgt, fHGT),
                  qval = p.adjust(pval,"BH") ) %>% arrange( qval ) %>%
        select( Name, nTrgt, nISG, pval, qval )
}

## Strongly-penalized geometric mean
spgm <- function( v )
{
    mean(log(v)) - length(v)
}

main <- function()
{
    pyLogging <- reticulate::import( "logging" )
    indra()$logger$setLevel( pyLogging$WARNING )
    
    ## Downstream transcription factors, identified by TFisg()
    ## The following line is built through:
    ##   TFisg() %>% pull(Name) %>% magrittr::extract(1:7) %>% dput()
    vDS <- c("NFKB1", "RELA", "STAT1", "IRF1", "STAT3", "STAT2", "IRF3")

    ## SLC2A4 -> CDKN1A, AR, MMP13 -> Exception
    G0 <- dijkstra( "SLC2A4", "IFITM3", fscore=spgm )
    
    ## Find connections from JAK2 and SIK3 to downstream TFs
    G1 <- dijkstra( "JAK2", vDS )
    G2 <- dijkstra( "SIK3", "IFITM3" )
    G3 <- dijkstra( "SIK3", "IFITM3", fscore=spgm )

    ##    save( G1, file="G1.RData" )
    load( "G1.RData" )
    plotPaths( G1$Path ) + ggthemes::scale_color_few()

}
