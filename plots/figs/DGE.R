## Plots comparing performance of DGE-based and mined gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/figs" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Short-hand for bold-face element_text()
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

main <- function()
{
    ## Load the associated LINCS metadata (drug names)
    M <- syn_csv( "syn11801537" ) %>% mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, Drug = name, URL = link, Target = target_name )

    ## Load the performance on mined and DGE sets. The following code was used to identify
    ##   relevant Synapse IDs:
    ## source( "../../R/resmine.R" )
    ## S <- synResults( "stats", Dataset == "ROSMAPpc" )
    X1 <- c(AB = "syn15589825", AC = "syn15589819", BC = "syn15589813") %>%
        map( syn_csv ) %>% bind_rows( .id = "Task" ) %>% mutate( Set = "Mined" )
    X2 <- c(AB = "syn16787402", AC = "syn16787438", BC = "syn16787423") %>%
        map( syn_csv ) %>% bind_rows( .id = "Task" ) %>% mutate( Set = "DGE" )

    ## Combine everything into a common data set and annotate with drug names
    XX <- bind_rows( X1, X2 ) %>% rename( URL = id ) %>% inner_join( M, by="URL" ) %>%
        select( Drug, Task, Set, AUC ) %>% spread( Set, AUC ) %>% na.omit() %>%
        mutate( zoom = FALSE )

    ## Isolate the results for the "zoom" panel
    ZZ <- XX %>% filter( Task == "AC" ) %>% mutate( zoom = TRUE )

    ## Plot the comparison of performance
    ggplot( XX, aes(x=Mined, y=DGE, color=Task) ) + theme_bw() +
        geom_abline( slope=1, intercept=0, lty="dashed", lwd=1.25 ) +
        geom_point( size=2 ) + ggthemes::scale_color_few() +
        geom_point( data=ZZ, size=2 ) +
        xlab( "AUC based on mined sets" ) + ylab( "AUC based on DGE sets" ) +
        ggforce::facet_zoom( y = (DGE > 0.65), zoom.data=zoom ) +
        ggrepel::geom_text_repel( data=ZZ, aes(label=Drug), color="Black", fontface="bold", size=4 ) +
        theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
              legend.title=etxt(14), legend.position = c(0.05,0.85),
              legend.background = element_rect(fill="white", color="#999999") ) +
        ggsave( "DGE-perf.png", width=11.4, height=6 )
}
