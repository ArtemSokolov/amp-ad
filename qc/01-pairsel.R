## Functionality used by 01-pairsel.Rmd
##
## by Artem Sokolov

## Downloads a synapse ID to local disk and returns its filename
syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/QC" )@filePath }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12) )
}

## Load relevant ROSMAP metadata
loadROSMAPMeta <- function()
{
    syn( "syn12299431" ) %>% read_tsv( col_types=cols(ID=col_character()) ) %>%
        select( ID, AOD, Braak ) %>% mutate( AOD = ifelse(AOD=="90+", "90", AOD) ) %>%
        mutate_at( "AOD", as.numeric )
}

## Download all AUC estimates for background sets
loadBKAUCs <- function()
{
    c( Current="syn11960715", `All-pairs`="syn12299006", `AOD-all`="syn12299007",
      `AOD-one`="syn12299008", `ID-reuse`="syn12299009" ) %>%
        map( ~read_csv( syn(.x), col_types=cols() ) ) %>% map( gather, Size, AUC ) %>%
        bind_rows( .id = "Estimate" ) %>% mutate_at( vars(Size), as.integer) %>%
        filter( Size <= 100 )
}

## Loads all the relevant pairs files
loadPairs <- function()
{
    c( Current="syn11942625", `All-pairs`="syn12299014", `AOD-all`="syn12299015",
      `AOD-one`="syn12299016", `ID-reuse`="syn12299017" ) %>%
        map( ~read_csv( syn(.x), col_types = cols(ID=col_character()) ) ) %>%
        bind_rows( .id = "Estimate" ) %>% select( ID, pair_index, Estimate )
}

## Loads a run file using the provided Synapse ID
loadRunFile <- function( synid )
{
    cat( "Loading", synid, "...\n" )
    syn( synid ) %>% read_csv( col_types = cols(ID=col_character()) )
}

## Loads all the runs files from a given Synapse directory
## Tags the result with the provided AUC estimate description tag
loadRunsDir <- function( dirID, tagName )
{
    ## Identify all the files present in a given directory
    str_c( "select id from file where parentId==\"", dirID, "\"" ) %>% synQuery %>%
        as_tibble %>% mutate( Pred = map( file.id, ~loadRunFile(.x) ) ) %>% unnest %>%
        mutate( Estimate = tagName )
}

## Loads all the runs files from all relevant Synapse directories
## Caches the results to . (current working directory)
## Uploads the caches to Synapse
cacheRunsAll <- function()
{
    ## Current AUC estimates
    loadRunsDir( "syn11942410", "Current" ) %>% select( -(`110`:`1000`) ) %>%
        write_csv( "0_runs_cache.csv.gz" )

    ## Alternative schemes of pair selection
    loadRunsDir( "syn12299034", "All-pairs" ) %>% write_csv( "1_runs_cache.csv.gz" )
    loadRunsDir( "syn12299035", "AOD-all" ) %>% write_csv( "2_runs_cache.csv.gz" )
    loadRunsDir( "syn12299036", "AOD-one" ) %>% write_csv( "3_runs_cache.csv.gz" )
    loadRunsDir( "syn12299037", "ID-reuse" ) %>% write_csv( "4_runs_cache.csv.gz" )

    ## Upload each cached file to Synapse
    str_c(0:4, "_runs_cache.csv.gz") %>% map(File, parentId="syn12299033") %>% map(synStore) 
}

## Computes AUC from LPOCV
## hi - prediction values associated with high-label example in the pair
## lo - prediction values associated with low-label example in the pair
auc.lpocv <- function( hi, lo )
{
    ( 0.5*sum(hi==lo) + sum(hi>lo) ) / length(hi)
}

## Computes leave-pair-out AUC values for each (Estimate, file.id, Size) tuple
##   in R using .lbl as the column of labels
computeAUC <- function( R, .lbl = "Braak" )
{
    ## Load the relevant ROSMAP metadata and match it up against the runs frame
    RR <- R %>% gather( Size, Pred, `10`:`100` ) %>% inner_join( loadROSMAPMeta() ) %>%
        mutate_at( "Size", as.integer )

    ## Order each pair based on the labels column
    cat( "Ordering pairs based on", .lbl, "...\n" )
    X <- RR %>% group_by( Estimate, file.id, Size, pair_index ) %>%
        mutate( Ord = c("Lo","Hi")[order(!!rlang::sym(.lbl))] ) %>% ungroup

    ## Separate out predictions associated with the lowest and highest samples in each pair
    ##  Compute AUC from the side-by-side comparison
    cat( "Computing AUC for", .lbl, "...\n" )
    X %>% select( Estimate, file.id, Size, pair_index, Pred, Ord ) %>%
        spread( Ord, Pred ) %>% group_by( Estimate, file.id, Size ) %>%
        summarize( AUC = auc.lpocv( Hi, Lo ) ) %>% ungroup
}

## Computes AUC values for AOD and Braak for each relevant scenario using
##   cached runs data frames. Caches the result to disk.
computeAllAUCs <- function()
{
    ## Load the cached versions of runs data
    RR <- c( "syn14716431", "syn14716433", "syn14716434", "syn14716436", "syn14716438" ) %>%
        map( ~read_csv(syn(.x), col_types=cols(ID=col_character())) )
    
    ## Compute Braak and AOD AUCs for each
    A1 <- map( RR, computeAUC, "Braak" )
    A2 <- map( RR, computeAUC, "AOD" )

    ## Combine all the AUC values and cache the combination to disk
    AA1 <- bind_rows(A1) %>% mutate( Labels = "Braak" )
    AA2 <- bind_rows(A2) %>% mutate( Labels = "AOD" )
    bind_rows(AA1, AA2) %>% write_csv( "cache-AUCs.csv.gz" )

    ## Store the cache to Synapse
    File( "cache-AUCs.csv.gz", parentId = "syn12299033" ) %>% synStore()
}

## Evaluates the impact of AOD as a confounding factor
mainAODConf <- function()
{
    ## Load cached AUC estimates
    AA <- syn( "syn14716566" ) %>% read_csv( col_types=cols() )

    ## Positive control:
    ## Compare Braak AUC estimates to the BTR output
    BB <- loadBKAUCs() %>% group_by( Estimate, Size ) %>%
        summarize( BTR = mean(AUC) ) %>% ungroup
    AA %>% filter( Labels == "Braak" ) %>% group_by( Estimate, Size ) %>%
        summarize( REst = mean(AUC) ) %>% ungroup() %>% inner_join( BB ) %>%
        group_by(Estimate) %>% summarize( Cor = cor(REst, BTR) )

    ## Create a summary plot
    gg <- ggplot( AA, aes( x = Size, y = AUC, color = Estimate ) ) +
        theme_bw() + geom_smooth( se=FALSE ) + bold_theme() +
        facet_wrap( ~Labels, nrow=1 )
    ggsave( "AOD-confounder.pdf", gg, width=9, height=4 )
}
