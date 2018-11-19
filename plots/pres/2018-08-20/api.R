## Functionality common to multiple scripts
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD" )@filePath }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bw() + theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12) )
}

