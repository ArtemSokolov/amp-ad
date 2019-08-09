## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

## Returns a score legend grob
legendGrobScore <- function( smx = 2.3, txt="Score", top_title=FALSE,
                            mar = margin(b=0.5, l=0.5, unit="cm"))
{
    X <- tibble( x=c(0,1), y=c(0,1), !!txt:=c(0,smx) )
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    if( top_title )
        gcb <- guide_colorbar(title.position = "top", barwidth=unit(4,"cm"))
    else
        gcb <- guide_colorbar(title.vjust = .85, barwidth=unit(4,"cm"))
    gg <- ggplot( X, aes(x=x, y=y, color=!!sym(txt)) ) + geom_bar(stat="identity") +
        scale_color_gradientn( colors=pal(100), limits=c(0, smx),
                              guide=gcb ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom", legend.margin=mar )
    cowplot::get_legend( gg )
}

Fig3A <- function()
{
    library( cowplot )

    ## Use a single drug for the exemplar panels
    X <- indexDGE() %>% filter( Dataset != "MAYO", Task == "AC" ) %>% resFromIndex() %>%
        retag() %>% unnest() %>% filter( Drug == "nvp-tae684" ) %>%
        mutate( Size = round(Size/10)*10 ) %>% select( Dataset, Size, AUC )

    ## Get the matching background
    XX <- indexBackground() %>% filter( Dataset != "MAYO", Task == "AC" ) %>% bkFromIndex() %>%
        unnest() %>% nest( AUC, .key="BK" ) %>% retag() %>% select( -Task ) %>%
        inner_join( X, ., by=c("Dataset","Size")) %>% arrange( desc(Dataset) ) %>%
        mutate_at( "Dataset", as_factor )
    BK <- XX %>% select( Dataset, BK ) %>% unnest()

    ## Generate faceted plots for the exemplar
    gg <- ggplot( BK, aes(x=AUC) ) + facet_wrap( ~Dataset, nrow=1 ) + geom_density( lwd=1.5 ) +
        theme_bw() + theme_bold() + xlab("") + ylab("") + ylim(c(0,10)) +
        geom_segment( aes(x=AUC, xend=AUC, y=0, yend=8), data=XX, color="red", lwd=1.5 ) +
        theme( axis.text = element_blank(), axis.ticks = element_blank(),
              plot.margin=unit(c(0,0.01,0.5,0), "npc"), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(), panel.border=element_rect(size=1.2) )

    ## Unit helpers
    npc <- function(x) {unit(x,"npc")}
    snpc <- function(x) {unit(x,"snpc")}

    ## Pre-defined coordinates
    x <- c(0.12, 0.24, 0.5)
    y <- c(0.25, 0.48, 0.55)
    z <- c(0.32, 0.52, 0.71, 0.9)
    yadj <- c(0,0.5) * snpc(0.3)
    xadj <- c(0.5,-0.5) * snpc(0.3)

    ## Graphical parameters
    gp <- grid::gpar(lwd=2)
    arw <- grid::arrow( angle=15, length=unit(0.15,"inches"), type="open" )
    
    ## Grob helpers
    rGrob <- function(f, ...) { grid::rectGrob(..., gp=lift(grid::gpar)(c(gp,fill=f))) }
    lGrob <- function(a, ...) { grid::linesGrob(..., gp=gp, arrow=a) }
    
    ## Define additional grobs
    grobs <- list()
    grobs[["r1"]] <- rGrob( "#FFFDDB", npc(x[1]), npc(y[1]), snpc(0.25), snpc(0.25))
    grobs[["r2"]] <- rGrob( "#FEA130", npc(x[2]), npc(y[1]), snpc(0.25), snpc(0.25))
    grobs[["r3"]] <- rGrob( "#FEE28F", npc(x[3]), npc(y[1]), snpc(0.25), snpc(0.25))

    grobs[["l1"]] <- lGrob(arw, c(0.12,x[1]), npc(c(y[3], y[1]))+yadj)
    grobs[["l2"]] <- lGrob(arw, c(x[2],x[2]), npc(c(y[2], y[1]))+yadj)
    grobs[["l2.1"]] <- lGrob(NULL, c(z[1],z[1]), npc(y[3:2]) )
    grobs[["l2.2"]] <- lGrob(NULL, c(z[2],z[2]), npc(y[3:2]) )
    grobs[["l2.3"]] <- lGrob(NULL, c(z[3],z[3]), npc(y[3:2]) )
    grobs[["l2.4"]] <- lGrob(NULL, c(z[4],z[4]), npc(y[3:2]) )
    grobs[["l2.h"]] <- lGrob(NULL, c(x[2],z[4]), npc(c(y[2],y[2])) )
    grobs[["l3"]] <- lGrob(arw, npc(c(x[2],x[3]))+xadj, npc(c(0.25, 0.25)) )

    ## Score legend
    slg <- legendGrobScore( txt="Score: -log10(p)", top_title=TRUE,
                           mar = margin(b=0, l=0, unit="cm") )
    
    ## Combine everything
    reduce( map(grobs, draw_grob), `+`, .init=ggdraw(gg) ) +
        draw_text( "ROSMAP", x[1], 0.07, fontface="bold", size=12 ) +
        draw_text( "MSBB Ave.", x[2], 0.07, fontface="bold", size=12 ) +
        draw_text( "Composite", x[3], 0.07, fontface="bold", size=12 ) +
        ##        draw_text( "Score: -log10(p)", (z[2] + z[4])/2, .1, size=14, fontface="bold" ) +
        draw_grob( slg, (z[2] + z[3])/2, -0.4 ) +
        draw_text( "Geometric Average", (z[2] + z[4])/2, y[2]-0.05, size=12, fontface="bold" ) +
        draw_text( "Geometric\nAverage", (x[2] + x[3])/2, y[1], size=12, fontface="bold" )
##        ggsave( "test.pdf", width=10, height=4 )
}

## Helper function that prepares a data frame for heatmap plotting
fprep <- function( .df ) {
    .df %>% select( Drug, ROSMAP, `MSBB ave` = MSBB, Composite ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" ) %>% as.matrix()
}

## Prepares row annotations (toxicity) for the heatmap
aprep <- function( .df ) {
    .df %>% select( Drug, Toxicity = IsToxic, Approval ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" )
}
        
## Produces a set of values to print over the heatmap
nprep <- function( .df ) {
    .n  <- fprep(.df) %>% round(2)
    storage.mode(.n) <- "character"
    .n[,1:2] <- ""
    .n
}

## Additional MoA information for FDA-approved compounds
FDA_MOA <- function()
{
    tribble( ~Drug, ~MoA,
            "baricitinib", "JAK family",
            "lapatinib", "EGFR, HER2",
            "regorafenib", "RET, VEGFRs, a.o.",
            "tofacitinib", "JAK family",
            "ruxolitinib", "JAK1/2",
            "dasatinib", "BCR-ABL, SRCfamily, a.o.",
            "ponatinib", "BCR-ABL, SRCfamily, a.o.",
            "bortezomib", "26S proteazome",
            "cabozantinib", "MET, VEGFRs, AXL, a.o.",
            "vorinostat", "HDAC1/2/3/6",
            "palbociclib", "CDK4/6",
            "nilotinib", "BCR-ABL",
            "ibrutinib", "BTK" )            
}

## Plots the top FDA-approved candidates ranked by composite score
Fig3B <- function()
{
    ## Fetch the composite score matrix and separate drugs in FDA-approved and non-approved
    XX <- DGEcomposite() %>% select( -(MSBB10:MSBB44), -LINCSID ) %>%
        mutate( IsApproved = ifelse( Approval %in% c("approved","vet_approved"),
                                    "FDA-Approved", "Non-Approved" ) ) %>%
        mutate_at( "IsToxic", recode, `0` = "Non-Toxic", `1` = "Toxic" ) %>%
        mutate_at( "Approval", str_to_title ) %>% left_join( FDA_MOA(), by="Drug" ) %>%
        mutate( Target = ifelse(is.na(MoA), Target, MoA) ) %>% select( -MoA ) %>%
        mutate( Drug = str_c(Drug, " (", Target, ") [", str_sub(Plate, 4, 4), "]"),
               Plate=NULL, Target=NULL ) %>% split( ., .$IsApproved ) %>%
        map( arrange, desc(Composite) ) %>% map( head, 15 )

    ## Set up the plotting mechanism
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    palA <- list( Toxicity = c("Toxic"="tomato", "Non-Toxic"="steelblue"),
                 Approval = set_names(ggthemes::few_pal()(3),
                                      c("Approved", "Experimental", "Investigational")) )
    fplot <- partial( pheatmap::pheatmap, cluster_rows=FALSE, cluster_cols=FALSE,
                     color=pal(100), fontsize=11, fontface="bold", gaps_col=2,
                     number_color="black", width=6.5, height=6, silent=TRUE,
                     annotation_colors=palA, legend=FALSE, annotation_legend=FALSE )
    
    ## Plot approved and non-approved drugs separately
    ggh <- map2( XX, c("FDA-approved\n", "Experimental and Investigational\n"),
                ~fplot(fprep(.x), display_numbers=nprep(.x), annotation_row=aprep(.x), main=.y) ) %>%
        map( pluck, "gtable" )

    ## Faux plot to generate the Score legend
    gl1 <- legendGrobScore( max(XX[[1]]$ROSMAP) )
    
    ## Faux plot to generate the Toxicity legend
    X2 <- XX[[1]] %>% select( Drug, Toxicity=IsToxic )
    gg2 <- ggplot( X2, aes(x=Drug, fill=Toxicity) ) + geom_bar(aes(y=1), stat="identity") +
        scale_fill_manual( values=palA$Toxicity,
                          guide=guide_legend(title.position="top") ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom", legend.margin=margin(b=0.5, l=0.5, unit="cm") )
    gl2 <- cowplot::get_legend( gg2 )
    
    ## Another faux plot to generate the Approval legend
    X3 <- bind_rows(XX) %>% select( Drug, Approval )
    gg3 <- ggplot( X3, aes(x=Drug, fill=Approval) ) +
        geom_bar( aes(y=1), stat="identity" ) +
        scale_fill_manual( values=palA$Approval,
                          guide=guide_legend(title.position="top") ) +
        theme( legend.title=etxt(12), legend.text=etxt(10),
              legend.position="bottom", legend.margin=margin(l=0, b=0.5, unit="cm"))
    gl3 <- cowplot::get_legend( gg3 )

    ## Put everything together
    gleg <- gridExtra::arrangeGrob( grobs = list(gl3,gl2,gl1), nrow=1,
                                   widths=c(4,2,3) )
    ly <- matrix( c(1,3,2,3), 2, 2 )
    gg <- gridExtra::arrangeGrob( grobs = c(ggh,list(Legend=gleg)), layout_matrix=ly,
                                 heights=c(0.9,0.1), widths=c(1.27,1) )

    ## Add custom annotations
    cowplot::ggdraw(gg) +
        cowplot::draw_text( "Drug (FDA-proposed MoA) [Plate Index]",
                          0.35, 0.93, fontface="bold", size=11 ) +
            cowplot::draw_text( "Drug (Vendor Target) [Plate Idx]",
                               0.85, 0.93, fontface="bold", size=11 )
##        cowplot::ggsave( str_c("Fig3-", Sys.Date(), ".pdf"), width=10, height=9 )
}

Fig3 <- function()
{
    F3A <- Fig3A()
    F3B <- Fig3B()

    ff <- cowplot::plot_grid( NULL, F3A, NULL, NULL, NULL, F3B, ncol=2,
                             rel_heights=c(0.25, 0.02, 0.75), rel_widths=c(0.05,1),
                             labels=c("A","","","","B",""), label_size=24 )
    cowplot::ggsave( str_c("Fig3-", Sys.Date(), ".pdf"), width=9, height=11 )
    cowplot::ggsave( str_c("Fig3-", Sys.Date(), ".png"), width=9, height=11 )
}
