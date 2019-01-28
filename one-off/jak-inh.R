## Analyses focused on JAK inhibitors
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()
syn <- synExtra::synDownloader( "/data/AMP-AD" )

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
    read_lines(fn) %>% str_split( "\\t" ) %>%
        set_names( map_chr(., nth, iName) ) %>%
        map( ~.x[-2:-1] )
}

## Custom ggplot theme that boldifies text elements
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12, hjust=0),
          strip.background = element_blank() )
}

## Isolates the slice of data associated with a particular drug and applies GSEA to it
## .df - data frame containing differential expression results
## drug - the drug of interest
## lGS - list of gene sets, with each gene set captured by a character vector
sliceGSEA <- function( .df, drug, lGS, lbl=drug )
{
    v <- .df %>% filter( Drug == drug ) %>% with( set_names(logFC, Gene) ) %>% keep( ~.!=0 )
##    print(fgsea::fgsea( lGS, v, 10000 ))
    fgsea::plotEnrichment( lGS[[1]], v ) + bold_theme() + labs( title=lbl ) +
        theme( title=etxt(15) ) + xlab( "Gene Rank" ) + ylab( "Enrichment Score" )
}

main1 <- function()
{
    RR1 <- syn( "syn15674107" ) %>% read_csv( col_types = cols() ) %>%
        mutate_at( "Drug", str_to_lower )
    RR2 <- syn( "syn17167348" ) %>% read_csv( col_types = cols() ) %>%
        mutate_at( "Drug", str_to_lower )

    ## Load the list of Interferome-stimulated genes by Schoggins, et al.
    isg <- syn( "syn11629935" ) %>% scan( what=character() ) %>% list() %>%
        set_names( "Interferome-Stimulated Genes" )

    ## Isolate the Jak inhibitors and compose their differential expression profiles
    gr <- list( sliceGSEA( RR1, "ruxolitinib", isg, "Ruxolitinib (p-value = 1e-4)" ),
               sliceGSEA( RR1, "tofacitinib", isg, "Tofacitinib (p-value = 1e-4)" ),
               sliceGSEA( RR2, "baricitinib", isg, "Baricitinib (p-value = 3e-3)" ) )

    gg <- gridExtra::arrangeGrob( grobs=gr, nrow=1 )
    ggsave( "JakInh-1.pdf", gg, width=12, height=4 )
}

main2 <- function()
{
    source( "../plots/figs/results.R" )

    ## Load the background
    BK <- syn_csv("syn15589816") %>% gather( Size, AUC ) %>% mutate_at( "Size", as.integer ) %>%
        nest( -Size, .key="BK" ) %>% mutate_at( "BK", map, pull, "AUC" )
    
    ## Identify the drugs of interest
    X <- syn_csv( "syn17435116" ) %>% annotateResults %>%
        filter( Drug %in% c("ruxolitinib", "baricitinib", "tofacitinib") ) %>%
        rename( nGenes = Size ) %>% mutate( Size = round(nGenes / 10) * 10 ) %>%
        inner_join( BK, by="Size" ) %>% select( -Size ) %>%
        mutate_at( "Drug", recode,
                  ruxolitinib="Ruxolitinib (p-value = 0.02)",
                  baricitinib="Baricitinib (p-value = 0.02)",
                  tofacitinib="Tofacitinib (p-value < 0.01)" )

    ## Separate the data frame into background and "drug of interest" portion
    RBK <- unnest(X) %>% select(-AUC) %>% rename(AUC=BK)
    
    gg <- ggplot() + theme_bw() + facet_wrap( ~Drug, nrow=1 ) + bold_theme() +
        geom_density( aes(x=AUC), RBK, fill="gray", alpha=0.65, lwd=1 ) +
        xlim( c(0.45, 0.9) ) + ylab("Density") +
        geom_segment( aes(x=AUC, xend=AUC), color="tomato", X, y=0, yend=20, lwd=1 )
    ggsave( "JakInh-2.pdf", gg, width=8, height=3 )
}
