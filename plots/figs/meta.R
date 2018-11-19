## Meta-analysis of performances across datasets/tasks
##
## by Artem Sokolov

source( "results.R" )

## Converts a tibble to matrix, using the specified column as rownames
tb2m <- function( TB, rnCol )
{ TB %>% as.data.frame() %>% column_to_rownames( rnCol ) %>% as.matrix() }

## Plots a similarity heatmap, where the similarity is computed as spearman correlation
##  between by-p-value rankings of drugs
pvalHeatmap <- function( XX, mainTitle = NA )
{
    SS <- XX %>% unnest %>% select( Tag, Drug, p_value ) %>% spread( Drug, p_value ) %>%
        tb2m( "Tag" ) %>% t() %>% cor( method="spearman" )
    pheatmap::pheatmap( SS, silent=TRUE, main=mainTitle,
                       breaks=seq( -0.5, 1, length.out=101 ) )$gtable
}

## Plots a similarity heatmap across each task individually
pvalHeatmap3 <- function( IDX, fnOut )
{
    XX <- resFromIndex( IDX ) %>%
        mutate( Tag = pmap_chr(list(Dataset, Region, Task), str_c, sep="_") )
    H3 <- split( XX, XX$Task ) %>% map( pvalHeatmap )
    HH <- gridExtra::arrangeGrob( grobs = H3, ncol=2 )
    ggsave( fnOut, HH, width=10, height=9 )
}

## General concordance in drug rankings across datasets
main <- function()
{
    indexMined() %>% pvalHeatmap3( "mined-concordance.png" )
    indexDGE() %>% pvalHeatmap3( "DGE-concordance.png" )
}

