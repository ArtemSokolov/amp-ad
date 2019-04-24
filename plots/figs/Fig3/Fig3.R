## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

## Helper function that prepares a data frame for heatmap plotting
fprep <- function( .df ) {
    .df %>% select( Drug, MSBB10, MSBB22, MSBB36, MSBB44, `MSBB ave` = MSBB,
                   ROSMAP, Composite ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" ) %>% as.matrix()
}

## Prepares row annotations (toxicity) for the heatmap
aprep <- function( .df ) {
    .df %>% select( Drug, Toxicity = IsToxic ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" )
}
        
## Produces a set of values to print over the heatmap
nprep <- function( .df ) {
    .n  <- fprep(.df) %>% round(2)
    storage.mode(.n) <- "character"
    .n[,1:6] <- ""
    .n
}

## Plots the top FDA-approved candidates ranked by composite score
fig3 <- function()
{
    ## Fetch the composite score matrix and separate drugs in FDA-approved and non-approved
    XX <- DGEcomposite() %>%
        mutate_at( "IsApproved", recode, `0` = "Non-Approved", `1` = "FDA-Approved" ) %>%
        mutate_at( "IsToxic", recode, `0` = "Non-Toxic", `1` = "Toxic" ) %>%
        mutate( Drug = str_c(Drug, " (", Target, ") [", str_sub(Plate, 4, 4), "]"),
               Plate=NULL, Target=NULL ) %>% split( ., .$IsApproved ) %>%
        map( arrange, desc(Composite) ) %>% map( head, 15 )

    ## Set up the plotting mechanism
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    palA <- list( Toxicity = c("Toxic"="tomato", "Non-Toxic"="steelblue") )
    fplot <- partial( pheatmap::pheatmap, cluster_rows=FALSE, cluster_cols=FALSE,
                     color=pal(100), fontsize=11, fontface="bold", gaps_col=c(4,6),
                     number_color="black", width=6.5, height=6, silent=TRUE,
                     annotation_colors=palA, legend=FALSE, annotation_legend=FALSE )
    
    ## Plot approved and non-approved drugs separately
    ggh <- map( XX, ~fplot(fprep(.x), display_numbers=nprep(.x), annotation_row=aprep(.x)) ) %>%
        map( pluck, "gtable" )

    ## Faux plot to generate the legend
    X1 <- XX[[1]] %>% select( Drug, Score=ROSMAP, Toxicity=IsToxic )
    smx <- max(X1$Score)
    gg0 <- ggplot( X1, aes(x=Drug, color=Score, fill=Toxicity) ) +
        geom_bar(aes(y=1), stat="identity") +
        scale_color_gradientn( colors=pal(100), limits=c(0, smx) ) +
        scale_fill_manual( values=palA$Toxicity ) +
        guides( color=guide_colorbar(title.vjust = .85) ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom" )
    ggl <- cowplot::get_legend( gg0 )

    ## Put everything together
    ly <- matrix(c(1,3,2,3),2,2)
    gg <- gridExtra::arrangeGrob( grobs = c(ggh,list(ggl)), layout_matrix=ly,
                                 heights=c(0.9,0.1), widths=c(0.96,1) )
    
    ## Write out to file
    ggsave( str_c("Fig3-", Sys.Date(), ".pdf"), gg, width=10.4, height=7 )
    ggsave( str_c("Fig3-", Sys.Date(), ".png"), gg, width=10.4, height=7 )
}
