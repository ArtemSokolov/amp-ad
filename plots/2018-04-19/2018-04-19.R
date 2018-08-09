## Plots for the 2018-04-19 meeting
##
## by Artem Sokolov

source( "../../R/freezemine.R" )

## Prefixes a filename with the directory corresponding to the analyses in this file
pfn <- function( fn ) {str_c( "2018-04-19/", fn )}

## Bold-face text element of requested size
## Used extensively in ggplot themes throughout the file
etxt <- function(s) {element_text(size = s, face="bold")}

## Makes a background distribution plot
## If provided, places a gene set of interest in its context
## Bk - a slice of relevant background data to plot
##        (same format as the output of synBkAll())
## FgVal - AUC value of the foreground set
##           set to NULL if only background should be plotted
## xrng - (optional) limits on the x axis
bkPlot <- function( Bk, FgVal = NULL, xrng = c(0.4, 0.75) )
{
    gg <- ggplot( Bk, aes(x=AUC) ) + theme_bw() + xlim( xrng ) +
        geom_density( fill="steelblue", alpha=0.3, lwd=1.5 ) +
        theme( axis.text.x = etxt(12), axis.title.x = etxt(14),
              axis.text.y = element_blank(), axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.text = etxt(12), legend.title = etxt(14) )
    if( is.null(FgVal) == FALSE )
        gg <- gg + geom_vline( xintercept = FgVal, color="red", lwd=1.5 )
    gg
}
    

## Plots of performance on background gene sets
bkMain <- function()
{
    ## Download all background results files
    ## Store them locally to save loading time:
    ##  BB <- synBkAll()
    ##  save( BB, file=pfn("bk-all.RData") )
    load( pfn("bk-all.RData") )

    ## Define common theme elements
    thm <- theme( axis.text = etxt(12), axis.title = etxt(14), strip.text = etxt(14),
                 legend.text = etxt(12), legend.title = etxt(14) )
    
    ## Plot an exemplar distribution for a pre-defined gene set size
    gg.s10 <- BB %>% filter( Size == 10, Dataset == "ROSMAP", Task == "AC" ) %>% bkPlot
    ggsave( pfn("bk-s10.png"), gg.s10, width=7, height=3 )

    ## Repeat for size = 20
    gg.s20 <- BB %>% filter( Size == 20, Dataset == "ROSMAP", Task == "AC" ) %>% bkPlot
    ggsave( pfn("bk-s20.png"), gg.s20, width=7, height=3 )
    
    ## Plot all results, showing how performance of background sets varies from
    ##  dataset to dataset and from task to task
    gg1 <- ggplot( BB, aes( x=Size, y=AUC, color=Task ) ) + theme_bw() +
        geom_smooth( lwd = 1.5 ) + facet_wrap( ~Dataset, nrow=2, ncol=3, dir="v" ) +
        thm + theme( legend.position = c(1,0), legend.justification = c(1,0) ) +
        guides( color=guide_legend(override.aes=list(fill="white")) )
    ggsave( pfn("bk1.png"), gg1, width=8, height=5 )

    ## "Zoom in" on MSBB
    MSBB <- BB %>% filter( Dataset != "ROSMAP" )
    gg2 <- ggplot( MSBB, aes( x=Size, y=AUC, color=Task ) ) + theme_bw() +
        geom_smooth( lwd = 1.5 ) + facet_wrap( ~Dataset, nrow=2, ncol=2, dir="v" ) +
        thm + guides( color=guide_legend(override.aes=list(fill="white")) )
    ggsave( pfn("bk2.png"), gg2, width=7, height=5 )
}

## Metformin sets against background on ROSMAP AC
metformin <- function()
{
    ## Load Metformin gene set performance on ROSMAP (AC), mapping gene set names accordingly
    oldNames <- c( "Metformin-upstream.txt", "Metformin-downstream.txt", "Metformin-SPNetwork.txt" )
    newNames <- c( "Mined", "Experimental", "Combined" )
    ##    XX <- synStatsMany( "Metformin", setNames(newNames, oldNames) )
    X <- synStats( "syn11969331", setNames(newNames,oldNames) )

    ## Find matching background
    load( pfn("bk-all.RData") )
    BB <- BB %>% filter( Dataset == "ROSMAP", Task == "AC" )

    ## Plot each of the three sets against background of matching size
    myrng <- c(0.55, 0.8)
    gg1 <- bkPlot( filter(BB, Size == 60), filter(X, Name=="Mined")$AUC, myrng )
    gg2 <- bkPlot( filter(BB, Size == 70), filter(X, Name=="Experimental")$AUC, myrng )
    gg3 <- bkPlot( filter(BB, Size == 470), filter(X, Name=="Combined")$AUC, myrng )
    gg <- gridExtra::arrangeGrob( gg1, gg2, gg3, ncol=1 )
    grid::grid.newpage()
    grid::grid.draw(gg)
    ggsave( pfn("metformin.png"), gg, width = 7, height = 6 )

    ## Look at DGE-derived Metformin set
    Y <- synStats( "syn12046273" )
    ggy <- bkPlot( filter(BB, Size == 900), filter(Y, Name=="metformin")$AUC, myrng )
    gg2 <- gridExtra::arrangeGrob( gg1, gg2, gg3, ggy, ncol=1 )
    ggsave( pfn("metformin2.png"), gg2, width = 7, height = 8 )
}
