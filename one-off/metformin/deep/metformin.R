library( tidyverse )

syn <- synExtra::synDownloader(".")

## Relevant bit from the metadata table:
## LIB039987_TRA00149187_S1	B2	10µM Metformin -1	Metformin
## LIB039987_TRA00149188_S2	B3	10µM Metformin -2	Metformin
## LIB039987_TRA00149189_S3	B4	10µM Mefformin -3	Metformin
## LIB039987_TRA00149190_S4	B5	10µM Glyburide -1	Glyburide
## LIB039987_TRA00149191_S5	B6	10µM Glyburide -2	Glyburide
## LIB039987_TRA00149192_S6	B7	10µM Glyburide -3	Glyburide
## LIB039987_TRA00149193_S7	B8	DMSO -1	DMSO
## LIB039987_TRA00149196_S10	C2	DMSO -2	DMSO
## LIB039987_TRA00149197_S11	C3	DMSO -3	DMSO

## Metadata annotations for Metformin vs. DMSO
metaMetformin <- function()
{
    tibble::tribble(
                ~Drug, ~Sample,
                "Metformin", "S1",
                "Metformin", "S2",
                "Metformin", "S3",
                "DMSO", "S7",
                "DMSO", "S10",
                "DMSO", "S11" ) %>%
        mutate_at( "Drug", factor, levels=c("DMSO","Metformin") ) %>%
        as.data.frame %>% column_to_rownames( "Sample" )
}

## Metadata annotations for Glyburide vs. DMSO
metaGlyburide <- function()
{
    tibble::tribble(
                ~Drug, ~Sample,
                "Glyburide", "S4",
                "Glyburide", "S5",
                "Glyburide", "S6",
                "DMSO", "S7",
                "DMSO", "S10",
                "DMSO", "S11" ) %>%
        mutate_at( "Drug", factor, levels=c("DMSO","Glyburide") ) %>%
        as.data.frame %>% column_to_rownames( "Sample" )
}


getCounts <- function(M, ...)
{
    M %>% filter( !duplicated(gene_name) ) %>%
        select( gene_name, ... ) %>%
        as.data.frame() %>% column_to_rownames("gene_name")
}

myEdgeR <- function( X, Y )
{
    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = X, samples = Y )
    dl <- edgeR::calcNormFactors( dl )
    mmx <- model.matrix( ~Drug, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential expression
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    RR <- edgeR::topTags( gf, nrow(X) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
    as_tibble(RR)
}

main <- function()
{
    X <- read_tsv( "combined-counts.tsv", col_types=cols() ) %>%
        filter( !is.na(gene_name) )

    ## Metformin vs. DMSO
    R1 <- myEdgeR( getCounts(X, S1, S2, S3, S7, S10, S11), metaMetformin() )
    vMet <- R1 %>% filter( logCPM > 0 ) %>% head(20) %>% pull(Gene)

    ## Compare to pre-existing gene sets
    ## synapser::synLogin()
    ## v1 <- syn("syn11617551") %>% scan( what=character() )
    ## v2 <- syn("syn11617549") %>% scan( what=character() )
    ## R1 %>% filter( Gene %in% v1 )
    ## R1 %>% filter( Gene %in% v2 )

    ## Glyburide vs. DMSO
    R2 <- myEdgeR( getCounts(X, S4, S5, S6, S7, S10, S11), metaGlyburide() )
    vGly <- R2 %>% head(16) %>% pull(Gene)

    ## Save everything into a common .gmt file
    cat( file="diabetes.gmt", sep="\n",
        str_flatten( c("Metformin","Metformin",vMet), "\t" ),
        str_flatten( c("Glyburide","Glyburide",vGly), "\t" ) )
}
