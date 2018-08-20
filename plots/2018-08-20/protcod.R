## Plots showing the effect of reducing to protein-coding regions
##
## by Artem Sokolov

source( "api.R" )
synapseLogin()

## Identifies the IDs of all relevant
idsBK <- function()
{
    structure(list(Dataset = c("ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "ROSMAP", "ROSMAP", "ROSMAP"),
                   Task = c("AB", "AC", "BC", "AB", "AC", "BC"),
                   BKid = c("syn15589822", "syn15589816", "syn15589810",
                            "syn15661345", "syn15661346", "syn15661344")),
              class = c("tbl_df", "tbl", "data.frame"),
              row.names = c(NA, -6L), .Names = c("Dataset", "Task", "BKid"))
}

mainPlot <- function()
{
    ## Load all the relevant entities
    X <- idsBK() %>% mutate( AUCs = map(BKid, ~read_csv(syn(.x), col_types=cols())) )

    ## Reshape everything into a single data frame
    XX <- X %>% mutate( AUCs = map( AUCs, gather, Size, AUC ) ) %>% unnest %>%
        mutate_at( "Size", as.integer ) %>% select( -BKid )

    ## Tweak the names by hand
    RR <- XX %>% mutate( `Gene Set` = c("ROSMAP" = "28.4k", "ROSMAPpc" = "ProtCod")[Dataset] )

    ## Compute summary distributions at key set sizes
    SS <- RR %>% filter( Size %in% c( 100, 300, 500, 700, 900 ) )
    
    ## Plot the results
    gg <- ggplot( RR, aes( x=Size, y=AUC, color=`Gene Set`) ) + theme_bw() +
        geom_boxplot( aes(group=interaction(Size, `Gene Set`)), data=SS ) +
        geom_smooth(se = FALSE) + facet_wrap( ~Task ) + bold_theme() +
        scale_color_manual( values=c("28.4k"="tomato", "ProtCod"="steelblue") ) +
        theme( legend.position="bottom" )
    ggsave( "plots/protcod.png", gg, width=9, height=4 )
}

