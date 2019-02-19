## Identifies wells treated with dsRNA on the DGE2 plate
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()

syn <- synExtra::synDownloader("~/data/AMP-AD/dsRNA/")

## edgeR analysis applied to drug-vs-control dichotomy
dvc_edgeR <- function( X, drugWells, ctrlWells, drugName )
{
    cat( "Comparing", drugName, "against controls\n" )
    Y1 <- bind_rows( tibble(Drug = drugName, Well = drugWells),
                    tibble(Drug = "Control", Well = ctrlWells) ) %>%
        mutate_at( "Drug", factor, levels=c("Control", drugName) ) %>%
        as.data.frame %>% column_to_rownames("Well")
    X1 <- X %>% select( HUGO, rownames(Y1) ) %>% as.data.frame %>% column_to_rownames("HUGO")

    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X1, samples = Y1 )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Drug, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    RR <- edgeR::topTags( gf, nrow(X1) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
    RR
}

main <- function()
{
    ## Load all the relevant files
    X <- syn( "syn17114495" ) %>% read_csv( col_types=cols() )
    Y <- syn( "syn17114496" ) %>% read_csv( col_types=cols() )

    ## Load the set of ISGs
    vISG <- syn( "syn11629935" ) %>% scan( what=character() )
    
    ## Isolate all wells currently marked as DMSO
    vDMSO <- Y %>% filter( Drug == "DMSO" ) %>% pull(Well)

    ## Downsample to the appropriate slice of the expression matrix
    X1 <- X %>% filter( HUGO %in% vISG ) %>% select( one_of(vDMSO) )

    ## Compute and plot the correlation among all well marked as DMSO
    S <- as.matrix(X1) %>% cor( method="sp" )
    pheatmap::pheatmap( S, filename="dsRNA.pdf" )

    ## ## ##
    
    ## Identify drug wells
    wDrug <- c("Ruxolitinib", "Tofacitinib") %>% set_names %>%
        map( ~filter(Y, Drug == .x, Concentration == 10) ) %>% map( pull, "Well" )

    ## Compute differential gene expression
    DFX <- wDrug %>% imap( ~dvc_edgeR(X, .x, vDMSO, .y) )

    ## Query for enrichment of ISGs
    GSEA <- DFX %>% map(with, set_names(logFC, Gene)) %>%
        map( ~fgsea::fgsea(list(ISG=vISG), .x, 1e4) )
}
