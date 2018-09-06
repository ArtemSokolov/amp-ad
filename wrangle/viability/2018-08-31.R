## Wrangles 2018-08-31 viability data
##
## Plates 1 & 2 contain dsRNA + drug data
## Plate 3 contains "drug only" data
## Plate 4 is a separate experiment and should be omitted from wrangling
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

## Loaders
syn <- synExtra::synDownloader( "/data/AMP-AD/viability" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }
load_plate <- function( synid )
{
    suppressWarnings( syn_csv( synid ) ) %>% rename( Row = X1 ) %>%
        gather( Col, Value, -Row ) %>% mutate_at( "Col", str_pad, 2, pad="0" ) %>%
        mutate( Well = str_c(Row, Col) ) %>% select( Well, Value )
}

## Annotates raw luminescene data with the provided metadata
## iPlate is the plate index
annot_plate <- function( P, M, iPlate )
{ filter( M, Plate == iPlate ) %>% inner_join( P, by="Well" ) %>% na.omit }

main <- function()
{
    ## Load the well metadata
    M <- syn_csv( "syn16784399" ) %>%
        select( Plate, Well = `Dispensed well`, Drug = `Fluid name`, Concentration )

    ## Load luminescence values
    P1 <- load_plate( "syn16784395" )
    P23 <- load_plate( "syn16784396" )

    ## Annotate each plate
    MP1 <- annot_plate( P1, M, 1 )
    MP2 <- annot_plate( P23, M, 2 )
    MP3 <- annot_plate( P23, M, 3 )

    ## Plates 1 & 2 contain dsRNA + drug data
    bind_rows( MP1, MP2 ) %>% write_csv( "Drug-dsRNA.csv" )

    ## Plate 3 contains "drug only" data
    write_csv( MP3, "Drug-only.csv" )

    ## Store to Synapse
    File( "Drug-dsRNA.csv", parent = "syn16784374") %>% synStore()
    File( "Drug-only.csv", parent = "syn16784374") %>% synStore()
}
