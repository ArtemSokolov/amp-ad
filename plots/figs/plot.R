## Plotting elements that are shared across multiple figures
##
## by Artem Sokolov

## Short-hand for bold element_text of desired size
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
theme_bold <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Define the dataset palette
dsPal <- function()
{
    c( ROSMAP = "#882256", MAYO = "#cc6677", `MSBB10` = "#89cced",
      `MSBB44` = "#332f85", `MSBB36` = "#13783d", `MSBB22` = "#ddcb76" )
}

