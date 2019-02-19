## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "results.R" )

## Plots the top FDA-approved candidates ranked by composite score
main <- function()
{
    ## Helper function that prepares a data frame for heatmap plotting
    fprep <- function( .df ) {
        .df %>% replace_na( list(Target="") ) %>%
            mutate( Drug = str_c(Drug, " (", Target, ")",
                                 " [", str_sub(Plate, 4, 4), "]") ) %>%
            select( Drug, MSBB10:Composite ) %>% as.data.frame() %>% column_to_rownames( "Drug" ) %>%
            head( 20 ) %>% as.matrix()
    }

    ## Produces a set of values to print over the heatmap
    nprep <- function( .p ) {
        .n  <- round( .p, 2 )
        storage.mode(.n) <- "character"
        .n[,1:6] <- ""
        .n
    }
    
    ## Define the palette
    flog10 <- function(n) {log10( 10 / seq( 10, 1, length.out=n ))}
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    
    ## Fetch the composite score matrix
    X <- DGE1new_composite() %>% select(Plate, Drug, Target, starts_with("MSBB"), ROSMAP, Composite) %>%
        rename( `MSBB ave` = MSBB )

    ## Separate drugs into FDA-approved and non-approved
    DB <- syn("syn16932412") %>% read_csv() %>% transmute( Drug = str_to_lower(Name) )
    P1 <- X %>% filter( Drug %in% DB$Drug ) %>% fprep()
    P2 <- X %>% filter( !(Drug %in% DB$Drug) ) %>% fprep()

    pheatmap::pheatmap( P1, cluster_rows=FALSE, cluster_cols=FALSE, color=pal(100),
                       fontsize=12, fontface="bold", gaps_col=c(4,6), display_numbers=nprep(P1),
                       number_color="black",
                       file="composite_FDA.png", width=5, height=7.5 )

    pheatmap::pheatmap( P2, cluster_rows=FALSE, cluster_cols=FALSE, color=pal(100),
                       fontsize=12, fontface="bold", gaps_col=5, display_numbers=nprep(P2),
                       number_color="black",
                       file="composite_nonapprv.png", width=5, height=7.5 )
}
