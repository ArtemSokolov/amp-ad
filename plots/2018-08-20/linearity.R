## Plots comparing linear vs. non-linear predictors
##
## by Artem Sokolov

source( "../../qc/05-linearity.R" )

compPlot <- function( ids )
{
    ## Load all the files
    XX <- ids %>% mutate( BKPerf = map( synid, ~read_csv(syn(.x), col_types = cols()) ) )
    
    ## Reshape individual score frames into the long format and consolidate everything
    RR <- XX %>% mutate( BKAUC = map( BKPerf, gather, Size, AUC ) ) %>%
        select( -synid, -BKPerf ) %>% unnest %>% mutate_at( "Size", as.integer )

    ## Make a summary plot
    ggplot( RR, aes(x=Size, y=AUC, color=Method) ) + theme_bw() +
        geom_smooth( se=FALSE ) + facet_grid( Region ~ Task ) +
        bold_theme() + scale_color_few( labels=c("Log.Reg.","Rnd.Frst.","XGBoost") )
}

main <- function()
{
    gg1 <- compPlot( idsROSMAP() )
    gg2 <- compPlot( idsMayo() )

    gg <- gridExtra::arrangeGrob( gg1, gg2, ncol=1, heights=c(1,1.5) )
    ggsave( "plots/lin1.png", gg, width=9, height=7 )

    gg3 <- compPlot( idsMSBB() )
    ggsave( "plots/lin2.png", gg3, width=9, height=7 )
}
