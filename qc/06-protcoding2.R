## Another pass at the analysis of protein coding regions
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
synapseLogin()

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD/QC/06" )
{ synGet( id, downloadLocation = dlc )@filePath }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bw() + theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12) )
}

## Identifies the IDs of all relevant
idsBK <- function()
{
    ## Structure composes through dput() of the following
    ## source( "../R/resmine.R" )
    ## X1 <- allSettings( "rosmap" ) %>% filter( Method == "sklLR" )
    ## X2 <- allSettings( "syn15660477" ) %>% filter( Method == "sklLR" )
    ## XX <- bind_rows( X1, X2 ) %>% select( -Region, -Strategy, -Method ) %>%
    ##     mutate( Files = map( settings_md5, synBySettingsMd5 ) ) %>%
    ##     mutate( BKid = map_chr(Files, ~filter(.x, name=="background_auc.csv")$id) ) %>%
    ##     select( -Files, -settings_md5 )
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
    ggplot( RR, aes( x=Size, y=AUC, color=`Gene Set`) ) + theme_bw() +
        geom_boxplot( aes(group=interaction(Size, `Gene Set`)), data=SS ) +
        geom_smooth(se = FALSE) + facet_wrap( ~Task ) + bold_theme() +
        scale_color_manual( values=c("28.4k"="tomato", "ProtCod"="steelblue") )
}
