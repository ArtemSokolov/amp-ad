library( tidyverse )

## Define Synapse downloader(s)
synapser::synLogin()
syn <- synExtra::synDownloader("./data", ifcollision="overwrite.local")
syn_csv <- function( id ) {syn(id) %>% read_csv(col_types=cols())}

## Removes wells that have fewer than threshold counts (total, across all genes) in them
## Returns filtered Y data frame
filterWells <- function( Y, X, thresh = 1e5 )
{
    wKeep <- X %>% summarize_at( vars(-HUGO), sum ) %>% gather( Well, TotalCounts ) %>%
        filter( TotalCounts >= thresh ) %>% pull(Well)
    filter( Y, Well %in% wKeep )
}

## Custom parser for TF->Target relationships
parseTF <- function(fn)
{
    read_lines( fn ) %>% str_split( "\\t" ) %>%
        set_names( map_chr(., first) ) %>% map( ~.x[-1:-2] ) %>%
        keep( ~length(.)>=5 )
}

## Wrangles DGE data into common expression and metadata matrices
## fTox - whether to filter wells by toxicity (total counts)
## fConc - whether to filter wells to max concentration
loadDGE <- function( fTox=FALSE, fConc=FALSE )
{
    ## DGE1
    X1 <- syn_csv( "syn18145755" ) %>% rename_at( vars(-HUGO), ~str_c("DGE1_",.) )
    Y1 <- syn_csv( "syn15673460" ) %>% mutate_at( "Well", ~str_c("DGE1_",.) )

    ## DGE2
    X2 <- syn_csv( "syn17115624" ) %>% rename_at( vars(-HUGO), ~str_c("DGE2_",.) )
    Y2 <- syn_csv( "syn17115625" ) %>% mutate_at( "Well", ~str_c("DGE2_",.) )

    ## Filter by toxicity (total counts), if applicable
    if( fTox )
    {
        Y1 <- filterWells( Y1, X1, 3e4 )
        Y2 <- filterWells( Y2, X2, 8.5e4 )
    }

    ## Combine everything into a single set of matrices
    ## Filter by concentration, if applicable
    Y <- bind_rows(Y1, Y2) %>% mutate_at( "Drug", str_to_lower )
    if( fConc ) Y <- filter( Y, Concentration == 10 )
    X <- inner_join(X1, X2, by="HUGO") %>% select( HUGO, Y$Well )

    list( X=X, Y=Y )
}

## Given a drug target, marks each well with whether it contains a binder or a non-binder
idBinders <- function( target, Wells, TAS )
{
    TAS %>% select( Drug = name, Target = symbol, tas ) %>%
        inner_join( Wells, by="Drug" ) %>% filter( Target == target ) %>%
        mutate( Binder = ifelse(tas %in% c(1,2,3), "Yes", "No") ) %>%
        select( -tas )
}

## Computes differential gene expression for binders vs. non-binders using edgeR
## mpi - explicitly model plate index
myEdgeR <- function( target, DGE, TAS, mpi=FALSE )
{
    ## Compose inputs for edgeR
    eRY <- idBinders(target, DGE$Y, TAS) %>%
        mutate_at("Binder", factor, levels=c("No","Yes")) %>%
        mutate(Plate = str_sub(Well, 1, 4)) %>%
        as.data.frame %>% column_to_rownames("Well")
    eRX <- DGE$X %>% select(HUGO, rownames(eRY)) %>%
        as.data.frame %>% column_to_rownames("HUGO")

    ## Create the design matrix and estimate dispersion
    dl <- edgeR::DGEList( counts = eRX, samples = eRY )
    dl <- edgeR::calcNormFactors( dl )
    if( mpi ) mmx <- model.matrix( ~Binder + Plate, data = dl$samples )
    else mmx <- model.matrix( ~Binder, data = dl$samples )
    dl <- edgeR::estimateDisp( dl, mmx )

    ## Compute differential gene expression
    gf <- edgeR::glmFit( dl, mmx ) %>% edgeR::glmLRT( coef = 2 )
    edgeR::topTags( gf, nrow(eRX) ) %>% as.data.frame %>% rownames_to_column( "Gene" )
}

## Enrichment analysis to infer TF activity based on target expression
myGSEA <- function( DFX, TF )
{
    fgsea <- function(v) { suppressWarnings(fgsea::fgseaMultilevel(TF, v)) }
##    fgsea <- function(v) { suppressWarnings(fgsea::fgsea(TF, v, 10000)) }
    with( DFX, set_names(logFC, Gene) ) %>% fgsea() %>% as_tibble() %>% arrange(pval)
}

jak2 <- function()
{
    ## TAS scores
    TAS <- syn_csv("syn18268627") %>% mutate_at("name", str_to_lower)
    
    ## Transcription Factor -> Target relationships
    ## Include the ISG list as a positive control
    ISG <- syn("syn11629935") %>% scan( what=character() )
    fnTF <- str_c("https://raw.githubusercontent.com/bioinfonerd/",
                  "Transcription-Factor-Databases/master/Ttrust_v2/human_TF_GSEA/all.gmt")
    TF <- parseTF(fnTF) %>% c( list(ISG=ISG) )

    ## Consider three different parameter settings
    ## fTox - Filtering wells by drug toxicity (i.e., total mRNA counts)
    ## fConc - Filtering wells to those with max drug concentration
    ## mpi - Explicitly model the plate index in edgeR

    ## Compose a global data matrix for each parameter setting
    ## Load and wrangle DGE data for each scenario
    X <- crossing( fTox=c(FALSE,TRUE), fConc=c(FALSE,TRUE) ) %>%
        mutate( DGE = map2(fTox, fConc, loadDGE) ) %>%
        crossing( mpi=c(FALSE,TRUE) )

    ## In each scenario, compute differential expression between
    ##   drugs that bind to JAK2 and those that don't
    Y <- mutate( X, DFX = map2(DGE, mpi, ~myEdgeR("JAK2", .x, TAS, .y)) )

    ## Perform GSEA to generate a ranking of TFs
    Z <- select( Y, -DGE ) %>% mutate( TF = map(DFX, myGSEA, TF) )

    ## Save all results to a single .csv file
    Z %>% select( -DFX ) %>% mutate_at( "TF", map, select, expr(-leadingEdge) ) %>%
        unnest() %>% write_csv( "jak2-TFs-mlvl.csv" )

    ## Short test for SIK3
#    X1 <- X %>% filter( fTox, fConc )
#    Y1 <- mutate( X1, DFX = map2(DGE, mpi, ~myEdgeR("SIK3", .x, TAS, .y)) )
#    Z1 <- select( Y1, -DGE ) %>% mutate( TF = map(DFX, myGSEA, TF) )
#    R1 <- Z1 %>% select( -DFX, -fTox, -fConc ) %>%
#        mutate_at( "TF", map, select, expr(-leadingEdge) ) %>% unnest()
}

mostToxic <- function()
{
    ## Load relevant files
    DGE <- loadDGE()
    DGEf <- loadDGE(TRUE,TRUE)
    TAS <- syn_csv("syn18268627") %>% mutate_at("name", str_to_lower)

    ## Identify binders and non-binders for each target
    v <- TAS$symbol %>% unique %>% set_names
    Z <- v %>% map( idBinders, DGE$Y, TAS )
    Zf <- v %>% map( idBinders, DGEf$Y, TAS )

    ## Compute total counts in each well and average the values across wells
    f <- function( X, W )
    {
        if( length(W) == 0 ) return(NA)
        X %>% select(W) %>% summarize_all(sum) %>%
            gather(Well, TotalCnt) %>% pull(TotalCnt) %>% mean
    }
    TC <- map(Z, ~f(DGE$X, .x$Well))
    TCf <- map(Zf, ~f(DGEf$X, .x$Well))

    ## Pull everything into a common dataframe
    R <- list( enframe(Z, "Target", "Full"),
              enframe(Zf, "Target", "Filtered"),
              enframe(TC, "Target", "AWCfull"),
              enframe(TCf, "Target", "AWCfilt") ) %>%
        reduce( inner_join, by="Target" ) %>%
        mutate_at( c("Full","Filtered"), map_int, nrow ) %>% unnest() %>%
        mutate( AWCdiff = AWCfilt - AWCfull )

    ## Identify a toxicity-related target that is well-represented
    R %>% arrange( desc(AWCdiff) ) %>% filter( Filtered >= 10 )
}
