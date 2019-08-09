## Testing DGE2: dsRNA vs. DMSO with edgeR
## Positive control: ISGs should be enriched
##
## by Artem Sokolov

library( tidyverse )

## Set up Synapse downloaders
synapser::synLogin()
syn <- synExtra::synDownloader(".", ifcollision="overwrite.local")

## Downloads counts data and matching metadata
## Reduces to the dichotomy of interest and performs differential gene expression
## synCounts - Synase ID of the counts matrix
## synMeta - Synapse ID of the metadata
## dich - dichotomy of interest, assumed to be in -/+ order
DGE_DFX <- function( synCounts, synMeta, dich )
{
    ## Load the metadata and identify wells that contain dsRNA and DMSO
    Y <- syn(synMeta) %>% read_csv() %>% filter( Drug %in% dich  )%>%
        filter( is.na(Concentration) | Concentration == max(Concentration, na.rm=TRUE) ) %>%
        select(Well, Drug)

    ## Download the counts data and reduce to the wells of interest
    X <- syn(synCounts) %>% read_csv() %>% select( HUGO, one_of(Y$Well) )
    
    ## Prepare the metadata for EdgeR
    ## Impose class ordering through the factor levels
    YER <- Y %>% mutate_at( "Drug", factor, dich ) %>% as.data.frame %>% column_to_rownames("Well")

    ## Prepare the counts data for EdgeR
    XER <- X %>% as.data.frame %>% column_to_rownames("HUGO")

    ## Create the design matrix and estimate dispersion
    DGE <- edgeR::DGEList( counts = XER, samples = YER )
    DGE <- edgeR::calcNormFactors( DGE )
    mmx <- model.matrix( ~Drug, data = DGE$samples )
    DGE <- edgeR::estimateDisp( DGE, mmx )

    ## Compute differential expression
    gf <- edgeR::glmFit( DGE, mmx ) %>% edgeR::glmLRT( coef = 2 )
    edgeR::topTags( gf, nrow(DGE$counts) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
}

main <- function()
{
    DGE_DFX1 <- partial( DGE_DFX, "syn18145755", "syn15673460" )
    DGE_DFX2 <- partial( DGE_DFX, "syn17115624", "syn17115625" ) 
    
    ## Compute differential gene expression for DGE1 and DGE2 sets separately
    DFX <- list( dsRNA1 = DGE_DFX1(c("Drug control", "dsRNA +lipo")),
                Ruxo1 = DGE_DFX1(c("Drug control", "ruxolitinib")),
                Tofa1 = DGE_DFX1(c("Drug control", "tofacitinib")),
                dsRNA2 = DGE_DFX2(c("DMSO", "dsRNAmi")),
                Ruxo2 = DGE_DFX2(c("DMSO", "Ruxolitinib")),
                Tofa2 = DGE_DFX2(c("DMSO", "Tofacitinib")),
                `Ruxo+Lipo` = DGE_DFX2(c("DMSO", "Ruxolitinib (Lipo)")),
                `Tofa+Lipo` = DGE_DFX2(c("DMSO", "Tofacitinib (Lipo)"))
                )
    
    ## Compose logFC signatures
    dfx <- map(DFX, with, set_names(logFC, Gene))
    
    ## Load the set of ISGs
    ISG <- syn("syn11629935") %>% read_lines()

    ## Gene set enrichment analysis
    fgsea <- partial( fgsea::fgsea, pathways=list(ISG), nperm=1e4 )
    map( dfx, fgsea ) %>% bind_rows( .id="Dich" ) %>% select( Dich, pval, NES, leadingEdge )

}
