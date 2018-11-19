## Volcano plots of performance estimates
##
## by Artem Sokolov

source( "results.R" )

## Short-hand for bold-face element_text()
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Volcano plots for all datasets/tasks
main <- function()
{
    X <- indexMined() %>% resFromIndex() %>% unnest() %>%
        mutate_at( "p_value", recode, `0` = 0.005 ) %>%
        mutate( Tag = glue::glue("{Dataset}_{Region}"), nl10p = -log10(p_value) )

    ggplot( X, aes(x=AUC, y=nl10p, color=Task) ) + theme_bw() + bold_theme() +
        facet_wrap( ~Tag, nrow=2, ncol=3 ) +
        ylab( "-log10( p value )" ) + ggthemes::scale_color_few() + geom_point() +
        ggsave( "volcano.png", width=12, height=7 )
}

