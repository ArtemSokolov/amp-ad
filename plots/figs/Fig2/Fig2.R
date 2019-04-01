## Performance on background gene sets
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

panelA <- function()
{
    BK <- indexBackground() %>% bkFromIndex() %>% retag() %>% unnest()

    ggplot( BK, aes(x=Size, y=AUC, color=Dataset) ) + theme_bw() + theme_bold() +
        facet_wrap( ~Task ) + geom_smooth( se=FALSE ) +
        scale_color_manual( values=dsPal() ) +
        theme( legend.position = "bottom" ) +
        ggsave( str_c("Fig2A-",Sys.Date(),".pdf"), width=9, height=4 ) +
        ggsave( str_c("Fig2A-",Sys.Date(),".png"), width=9, height=4 )
}

panelB <- function()
{
    library( ggridges )

    ## A closer look some of the following drugs (#non-MAYO datasets where p < 0.05):
    ##  lapatinib (3), baricitinib (MAYO+2), nvp-tae684 (4), torin1 (3)

    ## Identify performance for the drugs of interest
    vDrugs <- c("lapatinib", "nvp-tae684")
    X <- indexDGE() %>% filter( Dataset != "MAYO", Task == "AC" ) %>% resFromIndex() %>%
        retag() %>% unnest() %>% filter( Drug %in% vDrugs, p_value < 0.05 ) %>%
        mutate( Label = ifelse(p_value == 0, " p<0.01", str_c(" p=",p_value)),
               Size = round(Size/10)*10 ) %>%
        select( -Task, -Plate, -IsApproved, -LINCSID )

    ## Get the matching background
    XX <- indexBackground() %>% filter( Dataset != "MAYO", Task == "AC" ) %>% bkFromIndex() %>%
        unnest() %>% nest( AUC, .key="BK" ) %>% retag() %>% select( -Task ) %>% inner_join(X, .) %>%
        mutate( Name = glue::glue("{Drug} ({Target})") ) %>% select( -Drug, -Target ) %>%
        mutate( UID = glue::glue("{Name}-{Dataset}") ) %>%
        arrange( Name, AUC ) %>% mutate_at( "UID", as_factor )
    BK <- XX %>% select( UID, Name, Dataset, BK ) %>% unnest

    ## Generate the ridge plots
    ggplot( BK, aes(x=AUC, y=UID, fill=Dataset) ) +
        theme_ridges(center_axis_labels=TRUE) +
        geom_density_ridges2(scale=1.25, size=1, alpha=0.5) +
        geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(UID), yend=as.numeric(UID)+0.9),
                     data=XX, color="black", lwd=2 ) +
        geom_text( aes(x=AUC, y=as.numeric(UID)+0.7, label=Label, hjust=0),
                  data=XX, fontface="bold", size=4 ) +
        scale_y_discrete( labels = XX$Name, name=NULL ) +
        coord_cartesian(clip="off") +
        scale_fill_manual( values=dsPal() ) +
        theme( axis.text=etxt(12), axis.title=etxt(14), legend.position="bottom",
              legend.text=etxt(12), legend.title=etxt(12) ) +
        ggsave( str_c("Fig2B-",Sys.Date(),".pdf"), width=8, height=5 ) +
        ggsave( str_c("Fig2B-",Sys.Date(),".png"), width=8, height=5 )
}
