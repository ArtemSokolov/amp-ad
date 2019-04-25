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
        scale_color_manual( values=dsPal() )
#        theme( legend.position = "bottom" ) 
##        ggsave( str_c("Fig2A-",Sys.Date(),".pdf"), width=9, height=4 ) +
##        ggsave( str_c("Fig2A-",Sys.Date(),".png"), width=9, height=4 )
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
        theme( axis.text=etxt(12), axis.title=etxt(14), legend.text=etxt(12), legend.title=etxt(12),
              legend.direction="horizontal", legend.position="bottom",
              legend.box.margin=margin(c(0,100,0,0)) )
##        ggsave( str_c("Fig2B-",Sys.Date(),".pdf"), width=8, height=5 ) +
##        ggsave( str_c("Fig2B-",Sys.Date(),".png"), width=8, height=5 )
}

panelC <- function()
{
    load(syn( "syn18565498" ))

    v <- c( AB = "A-vs-B", AC = "A-vs-C", BC = "B-vs-C" )
    RS <- RS %>% mutate_at( "Task", recode, !!!v )
    RBK <- RBK %>% mutate_at( "Task", recode, !!!v )
    
    SBK <- RBK %>% group_by( Task ) %>% summarize_at( "AUC", list(BK=list) ) %>%
        inner_join( RS, ., by="Task" ) %>% mutate( pval = map2_dbl(AUC, BK, ~mean(.x < .y)) ) %>%
        mutate( sig = ifelse( pval < 0.05, "Yes", "No" ) ) %>% arrange( pval ) 
   
    ## Additional plotting elements
    xrng <- bind_rows(RBK, SBK) %>% pull( AUC ) %>% range
    
    ## Plot everything together
    ggplot() + theme_bw() + facet_wrap( ~Task, ncol=1 ) + guides( color=FALSE ) +
        ggthemes::scale_fill_few() + theme_bold() + ylab( "Density" ) +
        scale_x_continuous( limits = c(1,1.1) * xrng ) + ylim( c(0,15) ) +
        ggrepel::geom_text_repel( aes(x=AUC, label=Name, color=sig), SBK, y=3,
                                 nudge_y=100, fontface="bold" ) +
        geom_density( aes(x=AUC, fill=Task), RBK, alpha=0.65, lwd=1 ) +
        geom_segment( aes(x=AUC, xend=AUC, color=sig), SBK, y=0, yend=3, lwd=1 ) +
        scale_color_manual( values=c("Yes" = "red", "No" = "black") ) +
        theme( strip.background = element_blank(), strip.text.x = element_blank(),
              legend.position=c(0.98,0.98), legend.justification=c(1,1),
              legend.background=element_rect(fill=NA), legend.title.align=0.75 )
}

fig2 <- function()
{
    ## Plot individual panels
    pA <- panelA()
    pB <- panelB()
    pC <- panelC()

    ## Place everything onto the same figure
    fBC <- cowplot::plot_grid( pB, NULL, pC, nrow=1, rel_widths=c(1, 0.01, 1.25),
                              labels=c("","","C"), label_size=20,
                              hjust=0.5, vjust=1.5 )
    fABC <- cowplot::plot_grid( NULL, pA, NULL, fBC, nrow=2, rel_widths=c(0.02,1),
                               rel_heights=c(1,2), labels=c("A","","B",""), label_size=20,
                               hjust=-0.2, vjust=1.25 )
    
    ggsave( str_c("Fig2-",Sys.Date(),".pdf"), fABC, width=13, height=7 )
    ggsave( str_c("Fig2-",Sys.Date(),".png"), fABC, width=13, height=7 )
}
