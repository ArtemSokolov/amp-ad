## Figures focusing on DGE data results
##
## by Artem Sokolov

source( "results.R" )

## Retrieves drug names associated with Nienke's gene sets that were used
##  in the prediction setting
minedDrugs <- function()
{
    ## An exemplar results file that shows which drugs have predictions available
    R <- syn_csv( "syn15572216" ) %>% rename( url = id )

    ## The associated metadata (in particular, the drug names)
    M <- syn_csv( "syn11801537" ) %>% rename( url = link )

    inner_join( M, R, by="url" ) %>% .$name %>% str_to_lower
}

DGE1 <- function()
{
    ## Load the data
    X <- syn_csv( "syn15673461" )
    Y <- syn_csv( "syn15673460" ) %>%
        mutate( Conc = as.numeric( str_split( Concentration, "u", simplify=TRUE )[,1] ) ) %>%
        select( -Concentration )

    ## Set of evaluated drugs for which mined associations are available
    vMined <- minedDrugs() %>% c( "metformin" )

    ## Set of FDA approved Drugs
    vApprv <- syn_csv( "syn16932412" )$Name %>% str_to_lower
    
    ## Identify drugs that are toxic by looking at the total number of counts in each well
    ## Remove them from consideration
    wKeep <- X %>% summarize_at( -1, sum ) %>% gather( Well, TotalCounts ) %>%
        filter( TotalCounts > 1e5 ) %>% .$Well
    Y <- Y %>% filter( Well %in% wKeep )

    ## Tally the number of replicates for each drug / concentration combo
    Z <- Y %>% filter( !is.na( Conc ) ) %>% count( Drug, Conc ) %>%
        group_by( Drug ) %>%
        mutate( dfexp = ifelse( Conc == max(Conc) & n > 1, 1, 0 ) ) %>%
        ungroup() %>% 
        mutate( dfexp = ifelse( Conc == 1 & Drug == "ldn-193189", 1, dfexp ) )

    ## Prepare the matrices to plot
    P <- Z %>% select( -n ) %>% spread( Conc, dfexp, fill=0 ) %>% as.data.frame %>%
        column_to_rownames( "Drug" )
    N <- Z %>% select( -dfexp ) %>% spread( Conc, n, fill=0 ) %>% as.data.frame %>%
        column_to_rownames( "Drug" )

    ## Prepare additional annotations
    A <- Z %>% select( Drug ) %>% distinct() %>%
        mutate( Mined = ifelse( Drug %in% vMined, "Yes", "No" ),
               Approved = ifelse( Drug %in% vApprv, "Yes", "No" ) ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" )

    ## Prepare the palette
    pal <- list( Mined = c("Yes" = "steelblue", "No" = "gray"),
                Approved = c("Yes" = "steelblue", "No" = "gray") )

    ## Plot the summary
    hh <- pheatmap::pheatmap( P, cluster_cols=FALSE, cluster_rows=FALSE, color=c("white","tomato"),
                       legend=FALSE, display_numbers = N, number_color = "black", silent=TRUE,
                       fontsize_number=12, fontsize_row=12, fontsize_col=12,
                       annotation_row = A, annotation_colors=pal )

    ## Summarize the controls (all wells with NA concentration)
    S <- Y %>% filter( is.na(Conc) ) %>% count( Drug ) %>% rename( Treatment = Drug )
    g1 <- gridExtra::tableGrob( S, row=NULL )

    ## Put everything together into a single plot
    lyt <- rbind( c(1,2), c(1,NA) )
    gg <- gridExtra::arrangeGrob( grobs=list(hh$gtable, g1), layout_matrix=lyt,
                                 widths=c(2,1), heights=c(1,3))

    ggsave( "DGE1.pdf", gg, width=5, height=9 )
}

## Generates labels for Dataset/Region combo
## Inputs can be vectors
DRlbl <- function( ds, rg )
{
    ## Compose the mappings
    dsm <- c( "ROSMAPpc" = "ROSMAP", "MAYOpc" = "MAYO", "MSBBpc" = "MSBB" )
    rgm <- c("BM10", "BM22", "BM36", "BM44", "DLPFC") %>% set_names %>%
        c( "TempCortex" = "TCX", "cerebellum" = "CB" )

    ## Compose the labels
    str_c( dsm[ds], ".", rgm[rg] )
}

DGE2 <- function()
{
    ## Load the performance of mined sets
    X <- indexMined() %>% resFromIndex()

    ## Set of FDA approved Drugs
    vApprv <- syn_csv( "syn16932412" )$Name %>% str_to_lower
    
    ## List of drugs for which DGE was already performed
    vDGE <- syn_csv( "syn15673460" ) %>% .$Drug %>% unique

    ## Concentrate on the AC task
    ## Consider drugs we haven't profiled yet
    XX <- X %>% filter( Task == "AC" ) %>% unnest %>%
        filter( !(Drug %in% vDGE) ) %>%
        mutate( Tag = glue::glue("{Dataset}_{Region}"),
               Approved = ifelse( Drug %in% vApprv, "Yes", "No" ) ) %>%
        select( LINCSID, Drug, p_value, Tag, Approved ) %>%
        spread( Tag, p_value ) %>%
        mutate( MSBB = pmin(MSBB_BM10, MSBB_BM22, MSBB_BM36, MSBB_BM44) ) %>%
        select( -contains("MSBB_") ) %>%
        rename( MAYO = MAYO_TCX, ROSMAP = ROSMAP_DLPFC )

    ## Sort by (geometric) average p-value
    P <- XX %>% mutate( p_ave = pmap_dbl( list(MAYO, ROSMAP, MSBB), ~..1 * ..2 * ..3 ) ) %>%
        arrange( p_ave ) %>% select( -p_ave )

    ## Write to .csv
    P %>% filter( Approved == "Yes" ) %>% select( -Approved ) %>% write_csv( "DGE2-approved.csv" )
    P %>% filter( Approved == "No" ) %>% select(-Approved) %>% head(30) %>%
        write_csv( "DGE2-top30unapprv.csv" )
}
