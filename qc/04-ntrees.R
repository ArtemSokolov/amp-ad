## Evaluates the impact of #trees on performance of nonlinear predictors
##
## by Artem Sokolov

## Downloads a synapse ID to local disk and returns its filename
syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/QC" )@filePath }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bw() + theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12) )
}

getData <- function()
{
    data_frame( Forest = c( "Linear", "10trees", "100trees", "100trees_bal" ),
               synID = c( "syn14715369", "syn14715304", "syn14877617", "syn14880031" ) ) %>%
        mutate( fn = map_chr(synID, syn) ) %>%
        mutate( X = map(fn, read_csv, col_types=cols()) ) %>% select( -fn, -synID ) %>%
        mutate( X2 = map(X, gather, Size, AUC) ) %>% select(-X) %>% unnest(X2) %>%
        mutate_at( "Size", as.integer )
}
