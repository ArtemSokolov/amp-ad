## Figure(s) displaying results of cell viability assays
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/viability" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Short-hand for bold text of desired size s
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Viability from the first experiment
main1 <- function()
{
    XX <- syn_csv( "syn16784373" )
    SS <- XX %>% group_by( Drug ) %>% summarize_at( "Viability", mean ) %>% arrange( Viability )
    vOrder <- c( "Control_Lipo", "dsRNAmi", bind_rows(SS)$Drug ) %>% unique

    RR <- SS %>% mutate_at( "Drug", parse_factor, vOrder )

    ## Compute the reference value
    ref <- filter( RR, Drug == "dsRNAmi" )

    ggplot( RR, aes(x=Drug, y=Viability) ) + theme_bw() +
        geom_bar( stat="identity", fill="#88BDE6" ) +
        geom_hline( data=ref, aes(yintercept = Viability), color="#C84630", lwd=1.5, lty="dashed" ) +
        geom_point( data=bind_rows(XX, .id="XP") ) +
        theme( axis.text.y = etxt(12), axis.title = etxt(14),
              axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
              strip.text = etxt(12), strip.background = element_rect(fill="#f4f9fc") ) +
        ggsave( "viability1.png", width=10, height=4 )
}

main <- function()
{
    ## Load the datasets
    X1 <- syn_csv( "syn16784373" )
    X2 <- syn_csv( "syn16787521" ) %>% select( Drug, Viability = Value )
    XX <- list( `2018-08-18` = X1, `2018-08-31` = X2 )

    ## Compute averages and derive order of drugs
    SS <- map( XX, group_by, Drug ) %>% map( summarize_at, "Viability", mean ) %>%
        map( arrange, Viability )
    vOrder <- c( "Control_Lipo", "dsRNAmi", bind_rows(SS)$Drug ) %>% unique
              
    ## Rearrange all frames based on the derived order
    RR <- map( SS, mutate_at, "Drug", parse_factor, vOrder ) %>% bind_rows( .id = "XP" )

    ## Compute the reference value
    ref <- filter( RR, Drug == "dsRNAmi" )

    ## Identify missing values
    TT <- RR %>% spread( XP, Viability ) %>% gather( XP, Viability, -Drug ) %>%
        filter( is.na(Viability) ) %>% mutate( Viability = 0 )

    ## Highlight G02 well
    G02 <- X2 %>% filter( Drug == "ruxolitinib", Viability < 20000 ) %>%
        mutate( XP = "2018-08-31", Viability = Viability * 0.9 )
        
    ## Plot the results
    gg <- ggplot( RR, aes(x=Drug, y=Viability) ) + theme_bw() +
        facet_wrap( ~XP, ncol=1, scales="free_y", strip.position="right" ) +
        geom_bar( stat="identity", fill="#88BDE6" ) +
        geom_hline( data=ref, aes(yintercept = Viability), color="#C84630", lwd=1.5, lty="dashed" ) +
        geom_point( data=bind_rows(XX, .id="XP") ) +
        geom_text( data=TT, label="NA", vjust=0, fontface="bold", color="#C84630" ) +
        geom_text( data = G02, label = "G02", vjust = 1, fontface="bold" ) +
        theme( axis.text.y = etxt(12), axis.title = etxt(14),
              axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
              strip.text = etxt(12), strip.background = element_rect(fill="#f4f9fc") )
    ggsave( "viability.png", gg, width=12, height=7 )
}
