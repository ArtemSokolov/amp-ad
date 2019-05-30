## Disease-activated microglia
##
## by Artem Sokolov

source( "../lpocv.R" )

## Disease-associated microglia gene sets
damSets <- function()
{
    ## Fetch the raw supplemental files from the corresponding papers
    download.file( "https://ars.els-cdn.com/content/image/1-s2.0-S0092867417305780-mmc3.xlsx",
                  "28602351.xlsx" )
    download.file( str_c("https://www.ncbi.nlm.nih.gov/pmc/articles/",
                         "PMC5719893/bin/NIHMS901481-supplement-2.xlsx"), "28930663.xlsx" )
    download.file( str_c("https://onlinelibrary.wiley.com/action/downloadSupplement?",
                         "doi=10.1002%2Fglia.23572&file=GLIA23572-sup-0001-TableS6.xlsx"),
                  "30758077.xlsx" )

    ## Retrieve the signatures
    v1 <- openxlsx::read.xlsx("28602351.xlsx", na.strings="NaN") %>%
        set_names( c("Gene", "nUMI", "logFC", "nl10p", "nl10q") ) %>%
        as_tibble() %>% drop_na() %>% arrange( desc(nl10q) ) %>% head( 100 ) %>%
        pull(Gene) %>% str_to_upper()
    v2 <- openxlsx::read.xlsx("28930663.xlsx", sheet=9) %>% pull(1) %>% str_to_upper()
    v3 <- openxlsx::read.xlsx("30758077.xlsx", startRow=7, colNames=FALSE) %>% select(X1, X4) %>%
        drop_na() %>% pull(X1)

    ## Put everything into a single .gmt
    list( "PMID:28602351" = v1, "PMID:28930663" = v2, "PMID:30758077" = v3 ) %>%
        map( str_flatten, "\t" ) %>% imap( ~str_c(.y, "\tDisease-Activated-Microglia\t", .x) ) %>%
        lift_dl(cat)( file="DAM.gmt", sep="\n" )
}

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
  read_lines(fn) %>% str_split( "\\t" ) %>%
    set_names( map_chr(., nth, iName) ) %>%
    map( ~.x[-2:-1] )
}

## Spot-check performance of sets before feeding them to BTR
main <- function()
{
    X <- loadROSMAP()
    P <- loadPairs( "ROSMAP", "DLPFC" )
    DAM <- read_gmt( "DAM.gmt" ) %>% map( keep, ~.x %in% colnames(X) )

    RR <- map( DAM, evalGeneSetBK, 10, X, P, "AC" )
    save( RR, file="res-DAM.RData" )
}

