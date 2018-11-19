## Plots of top 10 drug candidates for each dataset/region
##
## by Artem Sokolov

source( "results.R" )

library( ggridges )

## Short-hand for bold-face element_text()
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

plotTop10 <- function( Z, myTitle )
{
    ZZ <- Z %>% arrange( AUC ) %>% mutate_at( "Name", as_factor )
    BK <- ZZ %>% select( Name, BK ) %>% unnest

    ggplot( BK, aes(x=AUC, y=Name), fill="steelblue" ) +
        theme_ridges(center_axis_labels=TRUE) +
        geom_density_ridges2(scale=1.2, size=1, alpha=0.5) +
        geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Name), yend=as.numeric(Name)+0.9 ),
                     data=ZZ, color="tomato", lwd=2 ) +
        geom_text( aes(x=AUC, y=as.numeric(Name)+0.5, label=Label, hjust=0),
                  data=ZZ, fontface="bold", size=4 ) +
        theme( axis.text = etxt(12), axis.title = etxt(14) ) +
        ggtitle( myTitle )
}

## Top10 plots for all datasets/regions
main <- function()
{
    ## Identifies the top 10 (by p value) drugs in a results data frame
    ## Breaks ties using AUC
    top10 <- function( R ) {arrange(R, p_value, desc(AUC)) %>% head(10)}

    ## Identify the background performance
    BKall <- bkFromIndex() %>% unnest() %>% nest(AUC, .key="BK")

    ## Load the results
    X <- indexMined() %>% resFromIndex() %>% mutate_at( "Results", map, top10 ) %>% unnest %>%
        mutate( Tag = glue::glue( "{Dataset}_{Region}" ),
               Name = glue::glue( "{Drug} ({Target})" ),
               Label = glue::glue( " {Size}" ), Size = round(Size/10)*10 )

    ## Match up drugs to their backgrounds in the context of dataset / region / task
    XX <- left_join( X, BKall ) %>% select( -Dataset, -Region, -Drug, -Target, -LINCSID, -Size )

    ## Plot everything
    PP <- XX %>% arrange(Task) %>% mutate( Title = glue::glue("{Tag}_{Task}") ) %>%
        nest( -Title, .key=Top10 ) %>%
        transmute( Task = str_sub(Title, -2, -1), Plot = map2(Top10, Title, plotTop10) )

    ## Compose figures for each task
    ## Store each into a separate .png
    FF <- PP %>% group_by(Task) %>%
        summarize( Fig = list(gridExtra::arrangeGrob( grobs=Plot, nrow=2, ncol=3 )) ) %>%
        mutate( Filename = glue::glue("Top10-{Task}.png") )
    with( FF, map2( Filename, Fig, ggsave, width=15, height=10 ) )
}

