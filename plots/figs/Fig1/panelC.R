## Makes an exemplar gene set vs. background performance plot
##
## by Artem Sokolov

library( tidyverse )

## Produces bold element_text of desired size
etxt <- function(s, ...) element_text(size=s, face="bold", ...)

## Figure shows the intuition of comparing gene set of interest against background performance
panelC <- function()
{
    GSBK_lbl <- data_frame( AUC = c(0.69, 0.55, 0.67), y = c(0.3, 2, 5),
                           lbl=c("p-value", "Random", "Gene Set\nOf Interest") )
    X <- read_csv( "panelC.csv", col_types=cols() )
    ggplot( X, aes(x=AUC) ) + theme_bw() + xlim( c(0.4,0.75) ) + 
        geom_density( fill="steelblue", alpha=0.3, lwd=1.5 ) + xlab("Performance") +
        geom_vline( xintercept = 0.67, color="tomato", lwd=1.5 ) +
        ggrepel::geom_text_repel( data=GSBK_lbl, aes(y=y, label=lbl), size=5, fontface="bold",
                        nudge_y = c(1.5,2.5,5), nudge_x = c(1, -0.1, 1),
                        color=c("black","black","tomato") ) +
        theme( axis.title.x = etxt(14), axis.title.y = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black", fill=NA, size=2) ) +
        ggsave( "panelC.png", width=5, height=2.5 )
}
