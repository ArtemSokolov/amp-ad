## Makes an exemplar gene set vs. background performance plot
##
## by Artem Sokolov

library( tidyverse )

## Produces bold element_text of desired size
etxt <- function(s, ...) element_text(size=s, face="bold", ...)

## Produces the exemplar plot using real background performance values
main <- function()
{
  X <- read_csv( "gs-bk.csv" )
  ggplot( X, aes(x=AUC) ) + theme_bw() + xlim( c(0.4,0.72) ) + 
    geom_density( fill="steelblue", alpha=0.3, lwd=1.5 ) + xlab("Performance") +
    geom_vline( xintercept = 0.68, color="red", lwd=1.5 ) +
    theme( axis.title.x = etxt(14), axis.title.y = element_blank(),
           axis.text = element_blank(), axis.ticks = element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.border = element_rect(colour="black", fill=NA, size=2) )
}

