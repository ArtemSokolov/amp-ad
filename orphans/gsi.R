## Evaluation of several specific gene sets
##
## by Artem Sokolov

## ############################################## ##
## NOTE: THE CODE IS OLD AND NEEDS TO BE UPDATED  ##
##   TO ACCOUNT FOR THE LATEST CHANGES TO lpocv.R ##
## ############################################## ##

## Evaluates a task and a gene set of interest in the context of background
##   sets of the same size
eval.GSI <- function( Xfull, vGS, vTask, nBK = 100, MainTag = "Foreground" )
{
    ## Generate background sets
    ## Append them to the true gene set of interest
    ##   and evaluate everything
    vGenes <- colnames(Xfull) %>% setdiff( c("ID", "Stage") )
    RR <- data_frame( i = 1:nBK ) %>%
        transmute( Feats = map( i, ~sample(vGenes, length(vGS)) ),
                  Tag = str_c("Background",i) ) %>%
        bind_rows( data_frame( Feats = list(vGS), Tag = MainTag ), . ) %>%
        mutate( Data = map( Feats, ~reduce.ds(Xfull, .x, vTask) ) ) %>%
        select( -Feats )
    
    ## Run LPOCV on background sets
    ## Combine everything into a single data frame and return
    RR %>% transmute( LPOCV = map2( Data, Tag, lpocv ) ) %>% unnest
}

main <- function()
{
    ## Load ROSMAP data
    RSM <- loadROSMAP()

    ## Define the tasks of interest
    vTask1 <- c("late","early")
    vTask2 <- c("late","intermediate")
    vTask3 <- c("intermediate","early")
    
    ## Define the gene sets of interest

    ## Apoptosis gene set provided by Kris Sarosiek
    vApoptosis <- c("BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "MCL1")

    ## SPARCS gene set provided by Russell Jenkins
    vSPARCS <- c("TRIM22", "TRIM38", "IL32", "SPATS2L", "EPHA3", "HERC3",
                 "ADAM19", "SERPINB9", "IFI44L", "F3", "BEND6", "AIG1",
                 "MSRB2", "TNFRSF9", "ANTXR1" )

    ## Traverse the tasks and evaluate each one
    eval.GSI(Xfull, vGS, vTask1, 100, "Apoptosis") %>% write_csv("Apoptosis-AC.csv")
    eval.GSI(Xfull, vGS, vTask2, 100, "Apoptosis") %>% write_csv("Apoptosis-BC.csv")
    eval.GSI(Xfull, vGS, vTask3, 100, "Apoptosis") %>% write_csv("Apoptosis-AB.csv")
}

## Plots the results of main()
plotResults <- function()
{
    ## Task names
    vTasks <- c( "Early-v-Late", "Interm-v-Late" )
    
    ## Load individual results files
    X1 <- read_csv("Apoptosis-AC.csv", col_types=cols()) %>% mutate( Task=vTasks[1] )
    X2 <- read_csv("Apoptosis-BC.csv", col_types=cols()) %>% mutate( Task=vTasks[2] )
    XX <- bind_rows(X1,X2) %>% mutate( Task = factor(Task,vTasks) )

    ## Separate GSI and Background performance values
    GSI <- XX %>% filter( Tag == "Apoptosis" ) %>% select( -Tag )
    BK <- XX %>% filter( grepl( "Background", Tag ) ) %>% select( -Tag )

    ## Plot background distributions with overlayed GSI performance
    library( ggridges )
    gg <- ggplot( BK, aes(x=AUC, y=Task) ) + theme_ridges(center_axis_labels=TRUE) +
        geom_density_ridges2( scale=0.9, size=1, alpha=0.5, fill="steelblue" ) +
        geom_segment( aes(x=AUC, xend=AUC, y=Task, yend=as.numeric(Task)+0.9),
                     data=GSI, color="red", lwd=2 ) +
        scale_y_discrete(expand=c(0.1,0)) +
        theme( axis.text = element_text(size=12,face="bold"),
              axis.title = element_text(size=14,face="bold") )

    ggsave( "apoptosis-Predict.pdf", gg, width=7, height=3.5 )
}

## Plots the distribution of apoptosis gene expression across Braak stages in ROSMAP
plotApoptosisROSMAP <- function()
{
    ## Load ROSMAP data and reduce to the gene set of interest
    vGS <- c("BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "MCL1")
    Xraw <- syn( "syn11738012" ) %>% readq_tsv
    X <- Xraw %>% select( Braak, one_of(vGS) ) %>%
        gather( Gene, Expression, -Braak )

    ## Show the distribution of each gene across the Braak stages
    gg <- ggplot( X, aes(x=Braak, y=Expression) ) + theme_bw() +
        stat_summary( fun.y="mean", geom="bar", alpha=0.75, fill="steelblue" ) +
        geom_jitter( width=0.05 ) + facet_wrap( ~Gene ) +
        theme( axis.text = element_text(size=12,face="bold"),
              axis.title = element_text(size=14,face="bold"),
              strip.text = element_text(size=14,face="bold"),
              strip.background = element_blank() )
    ggsave( "apoptosis-ROSMAP.pdf", gg, width=7, height=7 )
}

## Plots the distribution of apoptosis gene expression from DGE data
plotApoptosisDGE <- function()
{
    ## Load DGE counts and meta data
    vGS <- c("BAX", "BAK1", "BCL2", "BCL2L1", "BCL2L11", "MCL1")
    X <- syn("syn11948496") %>% read_csv( col_types=cols() )
    Y <- syn("syn11948497") %>% read_csv( col_types=cols() )

    ## Identify the wells of interest
    vControls <- c( "Drug control", "Lipo control" )
    vTreats <- c( "naked dsRNA", "dsRNA +lipo" )
    W <- Y %>% filter( Drug %in% c(vControls, vTreats) ) %>%
        select( Well, Treatment = Drug ) %>%
        filter( !(Well %in% c("L09","O14")) )

    ## Combine everything into a single data frame of interest
    XY <- X %>% filter( HUGO %in% vGS ) %>% select( Gene=HUGO, one_of(W$Well) ) %>%
        gather( Well, Expression, -Gene ) %>% inner_join(W) %>%
        mutate_at( vars(Expression), ~log2(.+1) ) %>%
        mutate( Treatment = factor(Treatment, c(vControls,vTreats)) )

    ## Show the distribution
    gg <- ggplot( XY, aes(x=Treatment, y=Expression) ) + theme_bw() +
        stat_summary( fun.y="mean", geom="bar", alpha=0.75, fill="steelblue" ) +
        geom_jitter( width=0.05 ) + facet_wrap( ~Gene ) +
        theme( axis.text.x = element_text(size=12,face="bold",angle=90,vjust=0.5,hjust=1),
              axis.text.y = element_text(size=12,face="bold"),
              axis.title = element_text(size=14,face="bold"),
              strip.text = element_text(size=14,face="bold"),
              strip.background = element_blank() )
    ggsave( "apoptosis-DGE.pdf", gg, width=7, height=7 )
}
