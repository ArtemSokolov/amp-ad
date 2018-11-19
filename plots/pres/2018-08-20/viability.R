## Cell viability under drug treatment
##
## by Artem Sokolov

library( tidyverse )

main <- function()
{
    ## Load the raw data
    X <- read_csv( "180818_Artems_DGE_CTG.csv", col_types = cols() ) %>%
        gather( Drug, Viability ) %>% na.omit()

    ## Compute averages and derive order of drugs
    S <- X %>% group_by( Drug ) %>% summarize( Viability = mean(Viability) ) %>% arrange( Viability )
    vFixed <- c( "Control_Lipo", "dsRNAmi" )
    vOrder <- c( vFixed, setdiff( S$Drug, vFixed ) )
    S <- S %>% mutate( Drug = factor( Drug, vOrder ) )

    ## Compute the reference value
    ref <- S %>% filter( Drug == "dsRNAmi" ) %>% .$Viability
    
    ## Short-hand for bold text of desired size s
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    
    ## Plot the results
    gg <- ggplot( S, aes(x=Drug, y=Viability) ) + theme_bw() +
        geom_bar( stat="identity", fill="steelblue", alpha=0.8 ) +
        geom_hline( yintercept = ref, color="red", lwd=1.5, lty="dashed" ) +
        geom_point( data=X ) +
        theme( axis.text.y = etxt(12), axis.title = etxt(14),
              axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5) )

    ggsave( "viability.png", gg, width=12, height=6.5 )
}
