source( "../lpocv.R" )

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bw() + theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12), strip.background=element_blank() )
}

main <- function()
{
    X <- loadROSMAP()
    P <- loadPairs("ROSMAP", "DLPFC")

    ## New gene set related to apoptosis
    vApt <- c( "APAF1", "BAD", "BAK1", "BAX", "BCL2", "BCL2A1",
              "BCL2L1", "BCL2L11", "BID", "BIK", "CASP3", "CASP7",
              "CASP8", "CASP9", "HRK", "PMAIP1", "BBC3", "BMF",
              "BCL2L2", "CYCS", "MCL1" )
    
    vSet1 <- c( "BCL2", "MCL1", "BAX", "BCL2L1", "BMF", "CASP8" )
    vSet2 <- c( "FADD", "FAS", "TNF", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C",
               "TNFRSF10D", "TNFRSF1A", "TNFRSF1B", "TNFSF10", "FASLG", "TNFSF9", "TRADD" )
    lSets <- list( Set1 = vSet1, Set2 = vSet2 )

    ## rApt <- evalGeneSetBK( vApt, 100, X, P )
    ## save( rApt, file="apoptosis-lpocv.RData" )
    ## load( "apoptosis-lpocv.RData" )
    
    ## R <- map( lSets, evalGeneSetBK, 100, X, P )
    ## save( R, file="sets-lpocv.RData" )
    load( "sets-lpocv.RData" )
    R <- map( R, mutate_at, "Task", recode, AB="A-vs-B", BC="B-vs-C", AC="A-vs-C" )

    ## Plot the performance across all three tasks
    plotGeneSetEval( R ) + ylab( "Gene Set" ) +
        ggsave( "apoptosis-performance.pdf", width=8, height=4.5 )

    ## Plot raw expression of all genes across the Braak stages in ROSMAP
    X1 <- X %>% select( Braak, one_of(vApt) ) %>% gather( Gene, Expression, -Braak ) %>%
        mutate_at( "Braak", as.factor )
    ggplot( X1, aes(x=Braak, y=Expression) ) + theme_bw() + bold_theme() +
        geom_boxplot() + facet_wrap( ~Gene, nrow=3, ncol=7, scales="free" ) +
        ggsave( "apoptosis-expression.pdf", width=12, height=6 )
}


