## Performance on background gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/figs" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD/QC/05" )
{ synGet( id, downloadLocation = dlc )$path }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Performance as a function of gene set size (ROSMAP only)
bkfig1 <- function()
{
    ## The results files are identified by the following code:
    ## source( "../../R/resmine.R" )
    ## allSettings( "rosmap" ) %>% mutate( Files = map( settings_md5, synBySettingsMd5 ) ) %>%
    ##   select( Task, Files ) %>% unnest %>% filter( name=="background_auc.csv" )
    X <- c( AB = "syn15589822", AC = "syn15589816", BC = "syn15589810" ) %>% map( syn_csv ) %>%
        map( gather, Size, AUC ) %>% bind_rows( .id = "Task" ) %>% mutate_at( "Size", as.integer )

    gg <- ggplot( X, aes( x = Size, y = AUC, color=Task ) ) + theme_bw() + bold_theme() +
        ggthemes::scale_color_few() + stat_smooth( span = 0.25, method="loess", se=FALSE, size=1.5 )
    ggsave( "bk-rosmap.png", gg, width=6.5, height=4.5 )
}

## Performance as a function of gene set size (add datasets)
bkfig2 <- function()
{
    ## The results files are identified by the following code:
    ## source( "../../R/resmine.R" )
    ## SS <- map( c("rosmap","msbb","mayo"), allSettings ) %>% bind_rows %>%
    ##  filter( Region != "cerebellum" ) %>% mutate( Files = map(settings_md5, synBySettingsMd5) )
    ## SS %>% mutate( Files = map( Files, filter, name == "background_auc.csv" ) ) %>%
    ##  select( -id ) %>% unnest() %>% select( Dataset, Region, Task, synid = id )
    S <- structure(list(Dataset = c("ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "MSBBpc", "MSBBpc", "MSBBpc",
                                    "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc", "MSBBpc",
                                    "MSBBpc", "MSBBpc", "MSBBpc", "MAYOpc", "MAYOpc", "MAYOpc"),
                        Region = c("DLPFC", "DLPFC", "DLPFC", "BM10", "BM36", "BM44", "BM22", "BM36",
                                   "BM10", "BM22", "BM10", "BM44", "BM36", "BM44", "BM22",
                                   "TempCortex", "TempCortex", "TempCortex"),
                        Task = c("AB", "AC", "BC", "AC", "AB", "BC", "BC", "AC", "AB", "AB", "BC",
                                 "AC", "BC", "AB", "AC", "BC", "AB", "AC"),
                        synid = c("syn15589822", "syn15589816", "syn15589810", "syn15584696",
                                  "syn15585576", "syn15582756", "syn15585289", "syn15581358",
                                  "syn15582188", "syn15586016", "syn15585432", "syn15583710",
                                  "syn15584001", "syn15583564", "syn15584563", "syn15572112", 
                                  "syn15570670", "syn15572002")), row.names = c(NA, -18L),
                   class = c("tbl_df", "tbl", "data.frame"))

    ## Load all the data
    X <- S %>% mutate( Data = map( synid, syn_csv ) )

    ## Prepare the data for plotting
    XX <- X %>% mutate_at( "Region", recode, DLPFC = "", TempCortex="" ) %>%
        mutate_at( "Dataset", str_sub, 1, -3 ) %>% mutate_at( "Dataset", recode, MSBB = "MSBB." ) %>%
        mutate( Dataset = str_c( Dataset, Region ) ) %>% select( -Region, -synid ) %>%
        mutate( Data = map( Data, gather, Size, AUC ) ) %>% unnest %>%
        mutate_at( "Size", as.integer )

    gg <- ggplot( XX, aes( x=Size, y=AUC, color=Task ) ) + theme_bw() + bold_theme() +
        facet_wrap( ~Dataset ) + ggthemes::scale_color_few() + geom_smooth( se = FALSE )
    ggsave( "bk-all.png", gg, width=9, height=5 )
}
