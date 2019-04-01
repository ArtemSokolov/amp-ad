## Out-of-scope analysis: single-cell data from the mouse brain
## External data provided by Scott Lipnick
##
## by Artem Sokolov

source( "../lpocv.R" )

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

## Generates a background distribution for 200 genes for each task for ROSMAP
BK200_ROSMAP <- function()
{
    X <- loadROSMAP()
    P <- loadPairs( "ROSMAP", "DLPFC" )

    ## Same random 200-gene sets 1000 times
    ## Use the last 200 genes names from the data as the length seed
    vBK <- colnames(X) %>% tail(200) %>% genBK( X, 1000 )
    RBK <- evalGeneSets( vBK, X, P )
    save( RBK, file="BK200-ROSMAP.RData" )
}

## Generates a background distribution for 200 genes for each task for ROSMAP
BK200_MSBB10 <- function()
{
    X <- loadMSBB( "BM10" )
    P <- loadPairs( "MSBB", "BM10" )

    ## Same random 200-gene sets 1000 times
    ## Use the last 200 genes names from the data as the length seed
    vBK <- colnames(X) %>% tail(200) %>% genBK( X, 1000 )
    RBK <- evalGeneSets( vBK, X, P )
    save( RBK, file="BK200-MSBB10.RData" )
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

## Evaluates a single set of cell-type speicifc gene sets against a single dataset
## fnXLS - filename of the Excel table containing gene sets
## X - dataset
## P - matching pairs for LPOCV
## RBK - performance on background sets, as computed by BK200_*()
## pfx - prefix to be used in the output filenames
eval1 <- function( fnXLS, X, P, RBK, pfx )
{
    ## Load the raw sheets and reduce to the set of genes that occurs in the data
    vGSI <- parseTable( fnXLS ) %>% map( mutate_at, "Gene", str_to_upper ) %>%
        map( filter, Gene %in% colnames(X) ) %>% map( head, 200 ) %>% map( pull, Gene )

    ## Evaluate the joint set, as well as gene sets associated with each cell type
    RS <- evalGeneSets( vGSI, X, P )
    save( RS, file=str_c(pfx, ".RData") )

    ## Plot the performance against BK200
    plotResults( RS, RBK ) + ggsave( str_c(pfx,".png"), width=8, height=6 )

    ## Show mircoglia expression across AC
    ##pp <- plotSlice( X, "AC", vGSI$MG )
    ##ggsave( "MG-AC.png", pp, width = 16, height = 6 )
}

main <- function()
{
    ## Cell-type specific gene sets
    fnS3 <- "440032-3.xlsx"     ## independent of age
    fnS5 <- "440032-5.xlsx"     ## age-associated gene sets
    
    ## ROSMAP
    XROSMAP <- loadROSMAP()
    PROSMAP <- loadPairs( "ROSMAP", "DLPFC" )
    load( "BK200-ROSMAP.RData" )
    eval1( fnS3, XROSMAP, PROSMAP, RBK, "TableS3-ROSMAP" )
    eval1( fnS5, XROSMAP, PROSMAP, RBK, "TableS5-ROSMAP" )

    ## MSBB10
    XMSBB10 <- loadMSBB("BM10")
    PMSBB10 <- loadPairs("MSBB", "BM10")
    load( "BK200-MSBB10.RData" )
    eval1( fnS3, XMSBB10, PMSBB10, RBK, "TableS3-MSBB10" )
    eval1( fnS5, XMSBB10, PMSBB10, RBK, "TableS5-MSBB10" )
}

## Trains a single model to predict A vs C in ROSMAP data
## Uses Microglia gene set from Table S3
MGsig <- function()
{
    ## Load ROSMAP data
    X <- loadROSMAP()
    
    ## Load the raw sheets and reduce to the set of genes that occurs in ROSMAP
    ## Isolate the set of genes associated with microglia
    Graw <- parseTable( "440032-3.xlsx" )
    gsi <- Graw[["MG"]] %>% mutate_at( "Gene", str_to_upper ) %>%
        filter( Gene %in% colnames(X) ) %>% head(200) %>% pull( Gene )

    ## Isolate the relevant data slice
    XY <- X %>% reduceToTask( "AC" ) %>% select( Label, one_of(gsi) )

    ## Split into features and labels
    XX <- XY %>% select( -Label ) %>% as.matrix
    yy <- XY$Label

    ## Train a model
    mdl <- LiblineaR( XX, yy, type=0 )

    ## Retrieve the weights
    W <- t(mdl$W) %>% as.data.frame %>% rownames_to_column( "Gene" ) %>% rename( Weight=V1 ) %>%
        filter( Gene != "Bias" ) %>% mutate_at( "Weight", ~-.x )

    arrange(W, desc(Weight)) %>% write_csv( "TableS3-MG-model.csv" )
}

## Trains a single model to predict A vs C in ROSMAP data
## Uses Endothelial Cell gene set from Table S5
ECsig <- function()
{
    ## Load ROSMAP data
    X <- loadROSMAP()
    P <- loadPairs( "ROSMAP", "DLPFC" )

    ## Load the relevant set of genes
    gsi <- read_csv( "2018-09-26_Seur_YX_vs_OX_EC_MAST_logfc_thresh_0_min_percent_0.csv" ) %>%
        select( -X1 ) %>% mutate_at( "Gene", str_to_upper ) %>% filter( Gene %in% colnames(X) ) %>%
        arrange( p_val_adj ) %>% head(200) %>% pull(Gene)

    ## Evaluate the gene set
    RS <- evalGeneSets( gsi, X, P )
    load( "BK200-ROSMAP.RData" )
    RBK %>% nest( -Task, .key="BK" ) %>% mutate_at( "BK", map, pull, "AUC" ) %>% inner_join( RS ) %>%
        select( -Name, -Feats ) %>% mutate( pval = map2_dbl( AUC, BK, ~mean(.x < .y) ) ) %>%
        select( -BK )
    
    ## Isolate the relevant data slice
    XY <- X %>% reduceToTask( "AB" ) %>% select( Label, one_of(gsi) )

    ## Split into features and labels
    XX <- XY %>% select( -Label ) %>% as.matrix
    yy <- XY$Label

    ## Train a model
    mdl <- LiblineaR( XX, yy, type=0 )
    
    ## Retrieve the weights
    W <- t(mdl$W) %>% as.data.frame %>% rownames_to_column( "Gene" ) %>% rename( Weight=V1 ) %>%
        filter( Gene != "Bias" )
    
    arrange(W, desc(Weight)) %>% write_csv( "2018-09-26-AvB-model.csv" )
}
