## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "results.R" )

## Geometric average
gavg <- function(v) {exp(mean(log(v)))}

## Geometric average of non-zero values in the vector v
gavg_nz <- function( v ) {keep( v, ~.!=0 ) %>% gavg}

## Retrieves a slice of DGE results relevant to the requested task
## Removes MAYO results, due to the observed batch effect
DGEslice <- function( task="AC" )
{
    indexDGE() %>% filter( Task == task, Dataset != "MAYO" ) %>%
        resFromIndex() %>% mutate_at("Region", recode, DLPFC="" ) %>%
        mutate_at( "Region", str_sub, 3, 5 ) %>%
        mutate( Dataset = str_c( Dataset, Region ) ) %>% unnest() %>%
        select( Dataset, Plate, Drug, Target, p_value )    
}

## The composite score is defined by three values (in the order of importance)
## 1. Geometric average of p-values
## 2. The total number of p-values < 0.01 across all datasets / regions
## 3. [Log] Geometric average of all p-values that are >= 0.01
DGEcomposite <- function( task="AC" )
{
    ## Compute a composite score
    DGEslice(task) %>% nest( Dataset, p_value, .key="Perf" ) %>%
        mutate( Individual = map(Perf, spread, Dataset, p_value),
               Composite = map(Perf, summarize_at, "p_value",
                               funs(CS1 = gavg, CS2 = sum(.==0), CS3 = gavg_nz)) ) %>%
        select(-Perf) %>% unnest() %>% arrange( CS1, desc(CS2), CS3 )
}

