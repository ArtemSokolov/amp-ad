## Performance summary plot
##
## by Artem Sokolov

source( "../R/resmine.R" )

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD/results" )
{ synGet( id, downloadLocation = dlc )@filePath }

## Given a settings md5, identifies and load the corresponding stats file
##   for Nienke's gene sets
statsByMd5 <- function( md5 )
{
    id <- synBySettingsMd5( md5 ) %>% filter( grepl("Nienke", name), parentName == "stats" ) %>%
        verify( nrow(.) == 1 ) %>% .$id
    syn(id) %>% read_csv( col_types=cols() )
}

## Plots the top candidates for a given task
plotTopCand <- function( X1 )
{
    etxt <- function( s, ... ) {element_text( size = s, face = "bold", ... )}
    
    stopifnot( length(unique(X1$Task)) == 1 )

    ## Isolate the relevant portion of the matrix
    ## Cap the AUC and p values, compute the log score
    X2 <- X1 %>% select( Tag, Dataset, AUC, p_value, Task, DGEprof ) %>%
        mutate( AUC = ifelse(AUC < 0.5, 0.5, AUC),
               p_value = ifelse(p_value == 0, 0.01, p_value) ) %>%
        mutate( nlp = -log10(p_value) )

    ## Compute the average performance across all datasets
    ## Rank the drugs based on this average and keep the top 20
    v20 <- X2 %>% group_by( Tag ) %>% summarize( Score = sum(nlp) ) %>%
        arrange( desc(Score) ) %>% head( 20 ) %>% .$Tag

    ## Plot the top 20 drugs
    X3 <- X2 %>% filter( Tag %in% v20 ) %>% mutate( Tag = factor(Tag, rev(v20)) ) %>% arrange( Tag )
    vCol <- X3 %>% select( Tag, DGEprof ) %>% distinct %>% .$DGEprof %>% ifelse( "red", "black" )
    ggplot( X3, aes( x=Dataset, y=Tag, size=AUC, color=p_value ) ) + facet_wrap( ~Task ) +
        geom_point() + theme_bw() + scale_radius( range=c(2,16) ) +
        scale_color_gradient( low="red", high="white", trans = "log", breaks = c(0.05, 0.3) ) +
        xlab("") + ylab("") +
        theme( axis.text.x = etxt(10), axis.text.y = etxt(10, color=vCol),
              axis.title = etxt(12), strip.text = etxt(12),
              legend.text = etxt(10), legend.title = etxt(12) )
}

main <- function()
{
    ## Identify md5 keys of all relevant results files
    dsMap <- c( "ROSMAPpc" = "ROSMAP", "MAYOpc" = "Mayo" )
    rgMap <- c( "DLPFC" = "", "cerebellum" = "CER", "TempCortex" = "TCX" )
    MD5 <- allSettings() %>% filter( Strategy=="strategy3", Method=="sklLR", Dataset=="MAYOpc" ) %>%
        mutate( Dataset = str_c( dsMap[Dataset], rgMap[Region] ) ) %>%
        select( -Strategy, -Method, -Region )

    ## Load the corresponding stats files
    Xraw <- MD5 %>% mutate( Stats = map( settings_md5, statsByMd5 ) ) %>%
        select( -settings_md5 ) %>% unnest( Stats )

    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( id = link, Drug = name, Target = target_name )

    ## Load the list of drugs that have been profiled by DGE
    vDGE <- scan( "DGE-profiled.txt", what=character(), sep='\n' ) %>% str_to_lower()
    
    ## Match the LINCS IDs against Drug names and targets
    X <- inner_join( MLINCS, Xraw, by="id" ) %>% select( -id, -description, -n_genes) %>%
        rename( nGenes = intersect ) %>%
        mutate( Tag = ifelse( is.na(Target), Drug, str_c( Drug,"(",Target,")" ) ),
               DGEprof = str_to_lower(Drug) %in% vDGE )

    ## Plot each task individually
    PP <- split( X, X$Task ) %>% map( plotTopCand )
    gg <- gridExtra::arrangeGrob( grobs = PP, nrow=1 )
    ggsave( "summary.pdf", gg, width = 14, height=11 )
}
