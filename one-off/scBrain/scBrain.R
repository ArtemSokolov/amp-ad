## Out-of-scope analysis: single-cell data from the mouse brain
## External data provided by Scott Lipnick
##
## by Artem Sokolov

source( "../R/lpocv.R" )

## Short-hand for bold-face element_text()
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = element_blank(),
#          strip.text = etxt(12, hjust=0, color="gray40"),
          strip.background = element_blank() )
}

## Parses one of the supplemental .xlsx files
##   and retrieves tables associated with the cell types of interest
parseTable <- function( fn, vCTI = c("OPC", "OLG", "ASC", "mNEUR", "EC", "PC", "MG", "MAC") )
{
    ## Read gene names for the cell types of interest
    openxlsx::getSheetNames( fn ) %>% match( vCTI, . ) %>%
        set_names( vCTI ) %>% map( openxlsx::read.xlsx, xlsxFile = fn )
}

## Generates a background distribution for 200 genes for each task
BK200 <- function()
{
    X <- loadROSMAP()
    vBK <- colnames(X) %>% tail(200) %>% genBK( X, 1000 )
    RBK <- evalGeneSets( vBK, X )
    save( RBK, file="BK200.RData" )
}

## Plots a result in the context of its background
plotResults <- function( RS, RBK )
{
    SBK <- RBK %>% nest( -Task, .key = "BK" ) %>% mutate_at( "BK", map, pull, "AUC" )
    RS <- inner_join( RS, SBK ) %>% mutate( pval = map2_dbl( AUC, BK, ~mean(.x < .y) ) ) %>%
        mutate( sig = ifelse( pval < 0.05, "Yes", "No" ) ) %>% arrange( pval )
    
    ## Additional plotting elements
    xrng <- bind_rows(RBK, RS) %>% pull( AUC ) %>% range
    
    ## Plot everything together
    ggplot() + theme_bw() + facet_wrap( ~Task, ncol=1 ) + guides( color=FALSE ) +
        ggthemes::scale_fill_few() + bold_theme() + ylab( "Density" ) +
        scale_x_continuous( limits = c(0.9,1.1) * xrng ) + ylim( c(0,13) ) +
        ggrepel::geom_text_repel( aes(x=AUC, label=Name), RS, y=3,
                                 nudge_y=100, fontface="bold" ) +
        geom_density( aes(x=AUC, fill=Task), RBK, alpha=0.65, lwd=1 ) +
        geom_segment( aes(x=AUC, xend=AUC, color=sig), RS, y=0, yend=3, lwd=1 ) +
        scale_color_manual( values=c("Yes" = "red", "No" = "black") )
}

## Plots a slice of gene expression for a given task
## X - dataset, as loaded by loadROSMAP()
## task - one of {"AB", "AC", "BC"}
## vgs - vector for gene names to look at
plotSlice <- function( X, task, vgs )
{
    ## Isolate the slice of interest
    X1 <- X %>% reduceToTask( task ) %>% reduceToGenes( vgs, . ) %>% arrange( Label )

    ## Median-center each gene
    X2 <- X1 %>% mutate_at( vars(-ID, -Label), ~(.x - median(.x)) )
    
    ## Separate into plotting elements
    A <- X2 %>% select( ID, Label ) %>% as.data.frame %>% column_to_rownames("ID")
    P <- X2 %>% select( -Label ) %>% as.data.frame %>% column_to_rownames( "ID" ) %>% as.matrix

    ## Split the matrix by label, cluster each portion and recombine
    vv <- split( X2, X2$Label ) %>% map( pull, "ID" )
    PP <- map(vv, ~P[.x,])
    SS <- map(PP, ~hclust(dist(.x))$order) %>% map2( PP, ~.y[.x,] )
    P2 <- lift(rbind)(SS)

    ## Plot the results
    pheatmap::pheatmap( P2, annotation_row = A, cluster_rows=FALSE, silent=TRUE,
                       color = colorRampPalette( list("blue","white","red") )(100),
                       show_rownames=FALSE, fontsize_col = 7,
                       annotation_colors = list( Label = c("pos"="tomato","neg"="steelblue") ) )
}

rankSlice <- function( X, task, vgs )
{
    ## Isolate the slice of interest
    X1 <- X %>% reduceToTask( task ) %>% reduceToGenes( vgs, . ) %>%
        gather( Gene, Value, -ID, -Label )

    ## Compute the t.test statistic for each gene
    S <- X1 %>% select( -ID ) %>% nest( Value, .key="Vals" ) %>%
        mutate_at( "Vals", map, pull, "Value" ) %>% spread( Label, Vals ) %>%
        mutate( S = map2_dbl( pos, neg, ~t.test(.x,.y)$statistic ) ) %>%
        arrange( S )

    ## Identify the top and bottom gene
    v <- S %>% slice( c(1,n()) ) %>% pull( "Gene" )
    X2 <- X1 %>% filter( Gene %in% v )

    ## Plot expression distributions for the selected genes
    gg <- ggplot( X2, aes( x = Value, color = Label ) ) +
        geom_density() + theme_bw() +
        scale_color_manual( values = c("pos"="tomato", "neg"="steelblue") ) +
        facet_wrap( ~Gene, ncol=2 )
    ggsave( "MG-dens.png", gg, width=9, height=3 )
}

## Cell-type specific gene sets (independent of age)
TableS3 <- function()
{
    X <- loadROSMAP()
    
    ## Load the raw sheets and reduce to the set of genes that occurs in ROSMAP
    Graw <- parseTable( "440032-3.xlsx" )
    G <- map( Graw, mutate_at, "Gene", str_to_upper ) %>% map( filter, Gene %in% colnames(X) )

    ## Evaluate the joint set, as well as gene sets associated with each cell type
    vJoint <- map( G, head, 27 ) %>% bind_rows %>% pull( Gene ) %>% unique %>% setdiff( "TMSB4X" )
    vGSI <- map( G, head, 200 ) %>% map( pull, Gene ) %>% c( Joint = list(vJoint) )
    RS3 <- evalGeneSets( vGSI, X )
    save( RS3, file="TableS3.RData" )

    ## Plot the performance against BK200
    load( "BK200.RData" )
    plotResults( RS3, RBK ) + ggsave( "TableS3.png", width=8, height=6 )

    ## Show mircoglia expression across AC
    pp <- plotSlice( X, "AC", vGSI$MG )
    ggsave( "MG-AC.png", pp, width = 16, height = 6 )
}

## Cell-type specific age-associated gene sets
TableS5 <- function()
{
    X <- loadROSMAP()
    
    ## Load the raw sheets and reduce to the set of genes that occurs in ROSMAP
    Graw <- parseTable( "440032-5.xlsx" )
    G <- map( Graw, mutate_at, "Gene", str_to_upper ) %>% map( filter, Gene %in% colnames(X) )

    ## Select the top genes from each cell type
    vJoint <- map( G, head, 41 ) %>% bind_rows %>% pull( Gene ) %>% unique
    vGSI <- map( G, head, 200 ) %>% map( pull, Gene ) %>% c( Joint = list(vJoint) )
    RS5 <- evalGeneSets( vGSI, X )
    save( RS5, file="TableS5.RData" )
    
    ## Plot the performance against BK200
    load( "BK200.RData" )
    plotResults( RS5, RBK ) + ggsave( "TableS5.png", width=8, height=6 )
}

## Common age-associated gene sets across multiple cell types in the mouse brain
TableS6 <- function()
{
    X <- loadROSMAP()

    ## Compose gene sets
    Zct <- cols(.default=col_double(), Gene=col_character())
    Z <- read_csv( "440032-6.csv", col_types=Zct) %>% mutate_at( "Gene", str_to_upper ) %>%
        select( -NRP, -ImmN, -TNC, -VSMC ) %>% filter( Gene %in% colnames(X) ) %>%
        gather( CellType, Value, -Gene ) %>% na.omit() %>% group_by( CellType ) %>%
        summarize( GSI = list(Gene) )

    ## ** The rest needs revision **
    
    ## Evaluate each gene set
    ##    RR <- Z %>% mutate( Results = map2(GSI, CellType, ~evalGeneSet(X, .x, 30, .y)) )
    ##    save( RR, file="TableS6.RData" )
    ## load( "TableS6.RData" )
    
    ## Plot the results
    ## gg <- RR %>% select( Results ) %>% unnest() %>% plotGeneSetEval()
    ## ggsave( "TableS6.pdf", gg, width=6, height=12 )
}

