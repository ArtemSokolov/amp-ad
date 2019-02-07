## Performance on background gene sets
##
## by Artem Sokolov

source( "../results.R" )

library( ggridges )

## Common mappings for dataset, region and task annotations

## Recodes dataset, region and task annotations
myRecode <- function( .df )
{
    ## Predefine the mappings
    rmap <- c( DLPFC="", TCX="", BM10="10", BM22="22", BM36="36", BM44="44" )
    tmap <- c( AB="A-vs-B", AC="A-vs-C", BC="B-vs-C" )

    .df %>% mutate_at( "Region", recode, !!!rmap ) %>%
        mutate_at( "Task", recode, !!!tmap ) %>%
        mutate( Dataset = ifelse(Dataset == "MSBB", glue::glue("{Dataset}({Region})"), Dataset) ) %>%
        select( -Region )
}

etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Define the dataset palette
dsPal <- function()
{
    c( ROSMAP = "#882256", MAYO = "#cc6677", `MSBB(10)` = "#89cced",
      `MSBB(44)` = "#332f85", `MSBB(36)` = "#13783d", `MSBB(22)` = "#ddcb76" )
}

panelA <- function()
{
    BK <- indexBackground() %>% bkFromIndex() %>% myRecode() %>% unnest()

    ggplot( BK, aes(x=Size, y=AUC, color=Dataset) ) + theme_bw() + bold_theme() +
        facet_wrap( ~Task ) + geom_smooth( se=FALSE ) +
        scale_color_manual( values=dsPal() ) +
        theme( legend.position = "bottom" ) + ggsave( "PanelA.png", width=9, height=4 )
}

panelB <- function()
{
    X <- indexMined() %>% resFromIndex() %>% myRecode() %>% unnest() %>%
        mutate_at( "p_value", recode, `0` = 0.005 ) %>% mutate( nl10p = -log10(p_value) )

    ggplot( X, aes(x=AUC, y=nl10p, color=Dataset) ) + theme_bw() + bold_theme() +
        facet_wrap( ~Task ) + ylab( "-log10( p value )" ) +
        ggthemes::scale_color_few() + geom_point() + guides( color=FALSE ) +
        theme( plot.margin = margin(l=5,r=10) ) + ggsave( "PanelB.png", width=9, height=3 )
}

panelC <- function()
{
    ## Identifies the top 10 (by p value) drugs in a results data frame
    ## Breaks ties using AUC
    top10 <- function( R ) {arrange(R, p_value, desc(AUC)) %>% head(10)}

    ## The matching plotting function
    plotTop10 <- function( Z, myTitle )
    {
        ZZ <- Z %>% arrange( AUC ) %>% mutate_at( "Name", as_factor )
        BK <- ZZ %>% select( Name, Dataset, BK ) %>% unnest

        ggplot( BK, aes(x=AUC, y=Name, fill=Dataset) ) +
            theme_ridges(center_axis_labels=TRUE) +
            geom_density_ridges2(scale=1.2, size=1, alpha=0.5) +
            geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Name), yend=as.numeric(Name)+0.9 ),
                         data=ZZ, color="tomato", lwd=2 ) +
            geom_text( aes(x=AUC, y=as.numeric(Name)+0.5, label=Label, hjust=0),
                      data=ZZ, fontface="bold", size=4 ) +
            theme( axis.text = etxt(12), axis.title = etxt(14),
                  legend.position = "bottom" ) + ggtitle( myTitle )
    }
    
    ## Identify the background performance
    BKall <- bkFromIndex() %>% unnest() %>% nest(AUC, .key="BK") %>% myRecode()

    ## Load all results
    X <- indexMined() %>% resFromIndex() %>% myRecode() %>% unnest() %>%
        nest( -Task, .key="Data" ) %>% mutate_at( "Data", map, top10 ) %>% unnest() %>%
        mutate( Name = glue::glue( "{Drug} ({Target})" ), Label = glue::glue( " {Size}" ),
               Size = round(Size/10)*10 ) %>% select( -LINCSID, -Drug, -Target )

    ## Match the results up against their backgrounds
    XX <- left_join( X, BKall ) %>% nest( -Task, .key=Top10 ) %>%
        mutate( Plot = map2(Top10, Task, plotTop10) )

    ## Compose the joint figure panel for all three tasks
    FF <- gridExtra::arrangeGrob( grobs=XX$Plot, ncol=1 )
    ggsave( "PanelC.png", FF, width=6.5, height=11 )
}
