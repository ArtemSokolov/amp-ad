## Comparison of results from Freeze against the new pair selection scheme that
##   minimizes AOD distance
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
library( ggthemes )
library( plotly )

synapseLogin()

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD" )
{ synGet( id, downloadLocation = dlc )@filePath }

## Downloads a mapping of LINCS URLs to drug names
## Returns a mapping that can be fed directly to synStats()
LINCSmap <- function()
{
    M <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        with( setNames(name, link) )
    M <- read_csv( s@filePath, col_types = cols() )
    setNames( M$name, M$link )
}

## Loads the requested stats file, and pads it with additional columns
synStats <- function( id, ... )
{
    cat( "Loading", id, "...\n" )
    syn( id ) %>% read_csv( col_types=cols() ) %>%
        select( URL = id, AUC ) %>% mutate( ... )
}

## Generates a plotly associated with a given data frame
## The data frame must contain columns x, y, Task, Name and Target
myplot <- function( X, xlabel, ylabel )
{
    ## Display the results
    etxt <- function(s, ...) {element_text(size=s, face="bold", ...)}

    gg <- ggplot( X, aes(color=Task, label = Name, label2 = Target ) ) +
        theme_bw() + scale_color_few( "medium" ) +
        theme( axis.title=etxt(14), axis.text=etxt(12),
              legend.text=etxt(12), legend.title=element_blank() ) +
        geom_point( aes(x=x, y=y) ) + xlab( xlabel ) + ylab( ylabel ) +
        geom_abline( slope=1, color="gray", lwd=1.25, lty="dashed" )
    
    ggplotly( gg, tooltip=c("label","label2"), width=800, height=600 ) %>%
        add_annotations( text="Task", xref="paper", yref="paper",
                        font = list( size=20, face="bold" ),
                        x=1.02, xanchor="left", y=0.8, yanchor="bottom",
                        legendtitle=TRUE, showarrow=FALSE ) %>%
        layout( legend=list(y=0.8, yanchor="top" ) )
}

main <- function()
{
    ## Identify results associated with each pair selection scheme
    ids <- list(
        ## Class calls, no AOD minimization
        CCno = c( AB = "syn11961620", AC = "syn11961532", BC = "syn11961676" ),
        ## Probabilities, no AOD minimization
        PRno = c( AB = "syn12977808", AC = "syn12978062", BC = "syn12978675" ),
        ## Probabilities, AOD minimization
        PRyes = c( AB = "syn12959508", AC = "syn12958954", BC = "syn12961416" )
    )

    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( URL = link, Name = name, Target = target_name )

    ## Load all the relevant stats files
    SS <- map( ids, map, synStats ) %>% map( bind_rows, .id="Task" ) %>%
        bind_rows( .id="PairSel" ) %>% inner_join( MLINCS ) %>% select( -URL ) %>%
        spread( PairSel, AUC )

    ## Generate all pair-wise plots
    SS %>% rename( x = CCno, y = PRno ) %>%
        myplot( "Class calls, no AOD minimization", "Probabilities, no AOD minimization" )
    SS %>% rename( x = PRno, y = PRyes ) %>%
        myplot( "Probabilities, no AOD minimization", "Probabilities, w/ AOD minimization" )
    SS %>% rename( x = CCno, y = PRyes ) %>%
        myplot( "Class calls, no AOD minimization", "Probabilities, w/ AOD minimization" )

    ## Store the plots to files
    ##    htmlwidgets::saveWidget( pp, "minAOD-comparison.html" )
}

## Summary list of performances for Steve
list4steve <- function()
{
    gm_mean = function(x, na.rm=TRUE){
        exp(mean(log(x)))
    }
    
    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( LINCSID = lincs_id, URL = link, Name = name, Target = target_name )

    ## Load all the relevant stats files
    XX <- map( c( AB = "syn12959508", AC = "syn12958954", BC = "syn12961416" ), syn ) %>%
        map( read_csv, col_types=cols() ) %>% bind_rows( .id = "Task" ) %>%
        select( URL = id, n_genes, p_value, Task )

    ## Join all the data into a single matrix and compute average performance
    RR <- inner_join( XX, MLINCS ) %>% select( -URL ) %>% spread( Task, p_value ) %>%
        rowwise %>% mutate( ave_pval = gm_mean( c(AB,AC,BC) ) ) %>% ungroup %>%
        select( Name, Target, ave_pval, LINCSID, n_genes,
               AB_pval = AB, AC_pval = AC, BC_pval = BC )

    arrange( RR, ave_pval ) %>% write_csv( "perf-summary.csv" )
}
