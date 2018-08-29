## Wrangles aligned counts table into a genes-by-well matrix
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/DGE" )$path }

## Consolidates "Raw" concentration values in the origin spreadsheet into a set of common categories
map.conc <- function( v )
{
    setNames( c("0.3uM", "10uM", "1uM", "1uM", "10uM", "3uM", NA),
             c("0.299991", "9.99000001", "1_M", "0.99990001", "10_M", "2.932453597", "2_g/ml") )[v]
}

## Function allows for examination of duplicate gene_name entries
## Used for exploration only; not used in the final wrangling
showDupe <- function( Z, g )
{
    Z %>% filter( gene_name == g ) %>% gather( Well, Value, -gene_id, -gene_name ) %>%
        spread( gene_id, Value )
}

## Pre-QC wrangling of the data
## Joins all the relevant matrices from the "Aligned" directory (syn11947012)
preQC <- function()
{
    ## Load ENSEMBL to HUGO map
    GM <- syn( "syn14236139" ) %>% read_csv( col_types=cols() ) %>%
        select( gene_id, gene_name, gene_biotype )
    
    ## Load the matrix, row names, column names and the well tag-ID map
    X <- syn( "syn11947013" ) %>% read_delim( " ", comment="%", col_types=cols() ) %>%
        select( rowIndex = 1, colIndex = 2, Value = 3 )
    rn <- syn( "syn11947015" ) %>% read_csv( col_names = "gene_id", col_types=cols() ) %>%
        mutate( rowIndex = 1:34947 )
    cn <- syn( "syn11947014" ) %>% read_csv( col_names = "barcode", col_types=cols() ) %>%
        mutate( colIndex = 1:384 )
    wmap <- syn( "syn11947016" ) %>% read_tsv( col_types=cols() )

    ## Data contains several duplicated gene_names
    ## Using showDupe(), we examined the duplicates and identified the gene_id with the highest
    ##   total number of transcripts. We keep this gene_id and drop all others in each duplicated case
    vDrop <- c( "ENSG00000268439", "ENSG00000274559", "ENSG00000280987", "ENSG00000283706",
               "ENSG00000284024", "ENSG00000130489", "ENSG00000158427", "ENSG00000213380" )
    
    ## Pull everything together into a single genes-by-wells matrix
    ## Map gene IDs to HUGO names and reduce to protein-coding space
    XX <- inner_join( X, rn, by="rowIndex" ) %>% inner_join( cn, by="colIndex" ) %>%
        inner_join( wmap, by="barcode" ) %>% select( gene_id, well, Value ) %>%
        spread( well, Value, fill=0L ) %>% inner_join( GM, by="gene_id" ) %>%
        filter( !(gene_id %in% vDrop), gene_biotype == "protein_coding" ) %>%
        select( -gene_biotype, -gene_id ) %>% select( HUGO = gene_name, everything() )

    ## Load the metadata and adjust a couple of things by hand:
    ## - ruxolitinib+dsRNA is just ruxolitinib
    ## - metformin concentration in K16 is 0.3
    ## - paropanib is a typo; drug name is pazopanib
    Y <- syn( "syn11947020" ) %>% read_tsv( col_types=cols() ) %>%
        select( Well = `Dispensed well`, Drug = `Fluid name`, RawConc = Concentration ) %>%
        mutate( Concentration = map.conc(RawConc) ) %>% select( -RawConc ) %>%
        mutate( Drug = ifelse( Drug == "ruxolitinib+dsRNA", "ruxolitinib", Drug ) ) %>%
        mutate( Drug = ifelse( Drug == "paropanib", "pazopanib", Drug ) ) %>%
        mutate( Concentration = ifelse( Well == "K16", "0.3uM", Concentration ) )

    ## Write everything to files
    XX %>% write_csv( "wrangled-counts.csv" )
    Y %>% write_csv( "wrangled-meta.csv" )
}

## Additional modifications to the data based on the output of QC analyses
postQC <- function()
{
    ## Load pre-QC matrices
    X <- syn( "syn11948496" ) %>% read_csv( col_types = cols() )
    Y <- syn( "syn11948497" ) %>% read_csv( col_types = cols() )

    ## Remove well P11 due to primer sequence issue
    ## Remove DMSO wells L09 and O14 because they cluster with Lipo controls
    vRemove <- c( "P11", "L09", "O14" )
    Y %>% filter( !(Well %in% vRemove) ) %>% write_csv( "postqc-meta.csv" )
    X %>% select( -one_of(vRemove) ) %>% write_csv( "postqc-counts.csv" )
}

## edgeR analysis applied to drug-vs-control dichotomy
dvc_edgeR <- function( drugName, X, Y )
{
    cat( "Comparing", drugName, "against controls\n" )
    
    Y1 <- filter( Y, Drug %in% c( drugName, "Drug control" ) ) %>% replace_na( list(Conc=0) ) %>%
        group_by( Drug ) %>% top_n( 1, Conc ) %>% ungroup %>% select( Well, Drug ) %>%
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

diffexp <- function()
{
    ## Load the data
    X <- syn( "syn15673461" ) %>% read_csv( col_types=cols() )
    Y <- syn( "syn15673460" ) %>% read_csv( col_types=cols() ) %>%
        mutate( Conc = as.numeric( str_split( Concentration, "u", simplify=TRUE )[,1] ) ) %>%
        select( -Concentration )

    ## Identify drugs that are toxic by looking at the total number of counts in each well
    ## Remove them from consideration
    wKeep <- X %>% summarize_at( -1, sum ) %>% gather( Well, TotalCounts ) %>%
        filter( TotalCounts > 1e5 ) %>% .$Well
    Y <- Y %>% filter( Well %in% wKeep )

    ## After the previous step, some of the drugs only have a single replicate at the top concentration
    ## Go down to the next highest concentration by removing the single-repcliate wells
    WW <- Y %>% group_by( Drug ) %>% filter( Conc == max(Conc) ) %>% mutate(nn=n()) %>% filter( nn==1 )
    Y <- Y %>% filter( !(Well %in% WW$Well) )

    ## Ensure that all top concentrations now have at least 2 replicates
    nn <- Y %>% group_by( Drug ) %>% filter( Conc == max(Conc) ) %>% mutate(nn=n()) %>% .$nn
    stopifnot( min(nn) > 1 )

    ## Clip the counts matrix to match the remaining wells in the metadata
    X <- X %>% select( HUGO, Y$Well )

    ## Compose the set of drugs to consider
    vDrugs <- c( "dsRNA +lipo", "Drug control", "Lipo control", "naked dsRNA", "LPS" ) %>%
        setdiff( unique(Y$Drug), . )
    
    ## Traverse the drugs and compute differential expression for each
    RR <- map( vDrugs, dvc_edgeR, X, Y ) %>% set_names( vDrugs )
    bind_rows( RR, .id = "Drug" ) %>% write_csv( "diffexp.csv.gz" )
}
