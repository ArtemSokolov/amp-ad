library( tidyverse )

## Define Synapse downloader(s)
synapser::synLogin()
syn <- synExtra::synDownloader(".", ifcollision="overwrite.local")

## Download differential gene expression
DFX1 <- syn("syn18145776") %>% read_csv(col_types=cols())
DFX2 <- syn("syn17167348") %>% read_csv(col_types=cols())

## Download ISG list (positive control)
ISG <- syn("syn11629935") %>% scan( what=character() )

## Load TF-target associations
## Source: https://github.com/bioinfonerd/Transcription-Factor-Databases/
##              blob/master/Ttrust_v2/TF_GSEA/all.gmt.gz
TF <- read_lines( "all.gmt.gz" ) %>% str_split( "\\t" ) %>%
        set_names( map_chr(., first) ) %>% map( ~.x[-1] ) %>%
        keep( ~length(.)>5 ) %>% c( list(ISG=ISG) )

## Enrichment of gene sets GS for drug d in differential expression matrix DFX
mygsea <- function( DFX, d, GS )
{
    DFX %>% filter( Drug==d ) %>% with( set_names(logFC, Gene) ) %>%
        fgsea::fgsea( GS, ., 10000 ) %>% as_tibble() %>% arrange(pval)
}

## Tofacitinib (DGE1)
Tofa1 <- mygsea( DFX1, "tofacitinib", TF )
## # A tibble: 286 x 8
##    pathway     pval   padj     ES   NES nMoreExtreme  size leadingEdge
##    <chr>      <dbl>  <dbl>  <dbl> <dbl>        <dbl> <int> <list>     
##  1 ISG     0.000100 0.0286 -0.501 -1.52            0   321 <chr [101]>
##  2 MSX1    0.000582 0.0550  0.930  1.82            1     6 <chr [4]>  
##  3 CIITA   0.000741 0.0550 -0.712 -1.74            5    28 <chr [12]> 
##  4 RFXANK  0.000961 0.0550 -0.822 -1.73            6    13 <chr [8]>  
##  5 RFXAP   0.000961 0.0550 -0.822 -1.73            6    13 <chr [8]>  
##  6 TBP     0.00161  0.0767 -0.764 -1.66           11    15 <chr [5]>  
##  7 RFX5    0.00241  0.0986 -0.758 -1.64           17    15 <chr [8]>  
##  8 PAX8    0.00454  0.162  -0.855 -1.61           30     8 <chr [3]>  
##  9 ELF1    0.00556  0.167  -0.848 -1.60           37     8 <chr [2]>  
## 10 MAZ     0.00583  0.167  -0.869 -1.59           38     7 <chr [2]>  

## Ruxolitinib (DGE1)
Ruxo1 <- mygsea( DFX1, "ruxolitinib", TF )
## # A tibble: 286 x 8
##    pathway     pval   padj     ES   NES nMoreExtreme  size leadingEdge
##    <chr>      <dbl>  <dbl>  <dbl> <dbl>        <dbl> <int> <list>     
##  1 CIITA   0.000143 0.0225 -0.766 -1.98            0    28 <chr [12]> 
##  2 RFXANK  0.000309 0.0225 -0.844 -1.85            1    13 <chr [9]>  
##  3 RFXAP   0.000309 0.0225 -0.844 -1.85            1    13 <chr [9]>  
##  4 ISG     0.000315 0.0225 -0.423 -1.47            2   321 <chr [98]> 
##  5 RFX5    0.000610 0.0349 -0.795 -1.80            3    15 <chr [9]>  
##  6 ZEB1    0.00390  0.186  -0.787 -1.69           24    12 <chr [3]>  
##  7 MYBL2   0.00604  0.247   0.829  1.75           22     8 <chr [4]>  
##  8 HNF4A   0.00764  0.266   0.501  1.61           19    41 <chr [11]> 
##  9 RUNX2   0.00837  0.266  -0.702 -1.61           54    16 <chr [3]>  
## 10 TFAP2C  0.0191   0.400  -0.774 -1.55          118     9 <chr [2]>  

## Tofacitinib (DGE2)
Tofa2 <- mygsea( DFX2, "Tofacitinib", TF )
## # A tibble: 286 x 8
##    pathway    pval  padj     ES   NES nMoreExtreme  size leadingEdge
##    <chr>     <dbl> <dbl>  <dbl> <dbl>        <dbl> <int> <list>     
##  1 MSC     0.00226 0.264 -0.823 -1.69           13    10 <chr [5]>  
##  2 ISG     0.00307 0.264 -0.397 -1.37           28   317 <chr [91]> 
##  3 TP53    0.00390 0.264 -0.450 -1.48           33   157 <chr [47]> 
##  4 HIF1A   0.00523 0.264 -0.501 -1.53           41    81 <chr [18]> 
##  5 PITX1   0.00548 0.264 -0.745 -1.66           34    14 <chr [5]>  
##  6 NFKB1   0.00555 0.264 -0.394 -1.36           51   287 <chr [71]> 
##  7 TRERF1  0.00729 0.264  0.870  1.66           29     6 <chr [3]>  
##  8 RELA    0.00739 0.264 -0.390 -1.34           68   283 <chr [70]> 
##  9 APC     0.0130  0.311  0.731  1.71           47    12 <chr [3]>  
## 10 ELF3    0.0134  0.311  0.885  1.61           55     5 <chr [2]>  

## Ruxolitinib (DGE2)
Ruxo2 <- mygsea( DFX2, "Ruxolitinib", TF )
## # A tibble: 286 x 8
##    pathway    pval  padj     ES   NES nMoreExtreme  size leadingEdge
##    <chr>     <dbl> <dbl>  <dbl> <dbl>        <dbl> <int> <list>     
##  1 RFXANK  0.00104 0.148 -0.814 -1.71            6    12 <chr [10]> 
##  2 RFXAP   0.00104 0.148 -0.814 -1.71            6    12 <chr [10]> 
##  3 RFX5    0.00519 0.495 -0.759 -1.66           35    14 <chr [10]> 
##  4 AHR     0.00977 0.550 -0.663 -1.59           71    22 <chr [10]> 
##  5 CTNNB1  0.0176  0.550 -0.643 -1.54          129    22 <chr [8]>  
##  6 CEBPE   0.0187  0.550 -0.841 -1.50          114     6 <chr [2]>  
##  7 PROX1   0.0204  0.550  0.724  1.65           68    10 <chr [7]>  
##  8 GATA1   0.0217  0.550 -0.540 -1.47          175    46 <chr [13]> 
##  9 MBD1    0.0219  0.550 -0.901 -1.46          127     4 <chr [1]>  
## 10 SRY     0.0228  0.550 -0.835 -1.49          139     6 <chr [2]>
## ...
## 16 ISG     0.0336  0.550 -0.381 -1.24          331   317 <chr [103]>

## Combine everything into a common data frame and examine TFs
##   that best explain ISGs
TF_ISG <- c( "NFKB1", "RELA", "STAT1", "STAT2", "STAT3", "IRF1", "IRF3" )
rlang::exprs( Tofa1, Tofa2, Ruxo1, Ruxo2 ) %>% set_names() %>% map(eval) %>%
    enframe( "Drug", "GSEA" ) %>% unnest() %>% filter( pathway %in% TF_ISG ) %>%
    select( -padj ) %>% arrange( pval )
## # A tibble: 28 x 8
##    Drug  pathway    pval     ES   NES nMoreExtreme  size leadingEdge
##    <chr> <chr>     <dbl>  <dbl> <dbl>        <dbl> <int> <list>     
##  1 Tofa2 NFKB1   0.00555 -0.394 -1.36           51   287 <chr [71]> 
##  2 Tofa2 RELA    0.00739 -0.390 -1.34           68   283 <chr [70]> 
##  3 Ruxo1 NFKB1   0.0199  -0.374 -1.29          186   282 <chr [60]> 
##  4 Ruxo1 RELA    0.0212  -0.373 -1.29          198   281 <chr [63]> 
##  5 Ruxo1 STAT1   0.0235  -0.464 -1.42          188    77 <chr [21]> 
##  6 Tofa1 IRF1    0.0241  -0.549 -1.45          209    47 <chr [13]> 
##  7 Tofa1 STAT3   0.0416  -0.442 -1.29          404   135 <chr [36]> 
##  8 Ruxo2 NFKB1   0.0435  -0.381 -1.23          428   287 <chr [74]> 
##  9 Tofa1 RELA    0.0472  -0.400 -1.21          469   281 <chr [76]> 
## 10 Tofa1 NFKB1   0.0526  -0.398 -1.20          523   282 <chr [75]> 

