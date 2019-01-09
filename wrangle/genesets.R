## Wrangling gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD", ifcollision = "overwrite.local" )$path }

## Extracts a gene set associated with dsRNA proteomics dataset
geneset_dsRNAprot <- function()
{
    ## Retrieve the raw data
    X <- syn( "syn11807753" ) %>% read_csv( col_types = cols() ) %>%
        select( Gene=`Gene Symbol`, Control1=`Control Replicate 1`, Control2=`Control Replicate 2`,
               dsRNAmi1 = `DsRNAmi Replicate 1`, dsRNAmi2 = `DsRNAmi Replicate 2` ) %>%
        mutate( Control = map2_dbl(Control1, Control2, ~mean(c(.x,.y),na.rm=TRUE)),
               dsRNAmi = map2_dbl(dsRNAmi1, dsRNAmi2, ~mean(c(.x,.y),na.rm=TRUE)) ) %>%
        mutate( FoldChange = dsRNAmi / Control )

    ## Corner case: protein not expressed in control leading infinite fold change
    v1 <- X %>% filter( Control == 0, dsRNAmi > 10 ) %>% .$Gene

    ## Consider remaining data and compute log2 fold change
    XFC <- X %>% filter( Control != 0 ) %>% mutate( lFC = log2(FoldChange) )

    ## Write the result to file and store it on Synapse
    vup <- XFC %>% filter( lFC > 1 ) %>% .$Gene
    vdn <- XFC %>% filter( lFC < -1 ) %>% .$Gene

    ## Write the results to their respective files and store them to Synapse
    c( v1, vup ) %>% unique %>% cat( file="ReNprot-up.txt", sep="\n" )
    unique( vdn ) %>% cat( file="ReNprot-dn.txt", sep="\n" )
    
    ## Combine the two lists and retrieve unique entries
    ## Write the result to a file and store the file to Synapse
    aa <- list( Type="Gene Set", Category="Interferome", Reference="LSP Experiment" )
    f1 <- File( "ReNprot-up.txt", parentId = "syn11629934" ) %>% synStore() %>% synSetAnnotations(aa)
    f2 <- File( "ReNprot-dn.txt", parentId = "syn11629934" ) %>% synStore() %>% synSetAnnotations(aa)
}

