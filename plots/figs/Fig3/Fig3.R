## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

panelA <- function()
{
    X <- indexDGE() %>% resFromIndex() %>% retag() %>% unnest() %>%
        mutate_at( "p_value", recode, `0` = 0.005 ) %>% mutate( nl10p = -log10(p_value) )

    ggplot( X, aes(x=AUC, y=nl10p, color=Dataset) ) + theme_bw() + theme_bold() +
        facet_wrap( ~Task, scales="free" ) + ylab( "-log10( p value )" ) +
        scale_color_manual( values=dsPal() ) +
        geom_point() + guides( color=FALSE ) +
        theme( plot.margin = margin(l=5,r=10) ) +
        ggsave( str_c("Fig3A-",Sys.Date(),".pdf"), width=12, height=4 ) +
        ggsave( str_c("Fig3A-",Sys.Date(),".png"), width=12, height=4 )
}


## Plots the top FDA-approved candidates ranked by composite score
panelsBC <- function()
{
    ## Helper function that prepares a data frame for heatmap plotting
    fprep <- function( .df ) {
        .df %>% arrange(desc(Composite)) %>% head(15) %>%
            mutate( Drug = str_c(Drug, " (", Target, ")",
                                 " [", str_sub(Plate, 4, 4), "]") ) %>%
            select( Drug, MSBB10, MSBB22, MSBB36, MSBB44, `MSBB ave` = MSBB,
                   ROSMAP, Composite ) %>%
            as.data.frame() %>% column_to_rownames( "Drug" ) %>% as.matrix()
    }

    ## Produces a set of values to print over the heatmap
    nprep <- function( .p ) {
        .n  <- round( .p, 2 )
        storage.mode(.n) <- "character"
        .n[,1:6] <- ""
        .n
    }
    
    ## Fetch the composite score matrix and separate drugs in FDA-approved and non-approved
    XX <- DGEcomposite() %>%
        mutate_at( "IsApproved", recode, `0` = "Non-Approved", `1` = "FDA-Approved" ) %>%
        split( ., .$IsApproved ) %>% map( fprep )

    ## Compose output filenames
    fns <- c( "FDA-Approved" = "Fig3B", "Non-Approved" = "Fig3C" ) %>%
        str_c( "-", Sys.Date(), ".pdf" )
    fnsPng <- str_replace( fns, "pdf", "png" )

    ## Set up the plotting mechanism
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    fplot <- partial( pheatmap::pheatmap, cluster_rows=FALSE, cluster_cols=FALSE,
                     color=pal(100), fontsize=11, fontface="bold", gaps_col=c(4,6),
                     number_color="black", width=5, height=6, silent=TRUE )

    ## Plot the heatmaps
    map2( XX, fns, ~fplot(.x, display_numbers=nprep(.x), file=.y) )
    map2( XX, fnsPng, ~fplot(.x, display_numbers=nprep(.x), file=.y) )
}
