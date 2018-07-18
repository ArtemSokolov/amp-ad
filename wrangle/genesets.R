## Wrangling gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/" )@filePath }

## Extracts a gene set associated with dsRNA
geneset.dsRNA <- function()
{
    ## Retrieve the raw data
    X <- syn( "syn11807753" ) %>% read_csv( col_types = cols() ) %>%
        select( Gene=`Gene Symbol`, Control1=`Control Replicate 1`, Control2=`Control Replicate 2`,
               dsRNAmi1 = `DsRNAmi Replicate 1`, dsRNAmi2 = `DsRNAmi Replicate 2` ) %>%
        rowwise %>% mutate( Control = mean( c(Control1, Control2), na.rm=TRUE ),
                           dsRNAmi = mean( c(dsRNAmi1, dsRNAmi2), na.rm=TRUE ) ) %>% ungroup %>%
        mutate( FoldChange = dsRNAmi / Control )

    ## Corner case: protein not expressed in control leading infinite fold change
    v1 <- X %>% filter( Control == 0, dsRNAmi > 10 ) %>% magrittr::extract2( "Gene" )

    ## Consider remaining data and split the genes into up- and down- sets
    XFC <- X %>% filter( Control != 0 ) %>% select( Gene, FoldChange ) %>%
        mutate( lFC = log2( FoldChange ) )
    vup <- XFC %>% filter( lFC > 1 ) %>% magrittr::extract2( "Gene" )
    vdn <- XFC %>% filter( lFC < -1 ) %>% magrittr::extract2( "Gene" )

    ## Write the results to their respective files and store them to Synapse
    c( v1, vup ) %>% unique %>% cat( file="ReNprot-up.txt", sep="\n" )
    unique( vdn ) %>% cat( file="ReNprot-dn.txt", sep="\n" )
    
    ## Combine the two lists and retrieve unique entries
    ## Write the result to a file and store the file to Synapse
    f1 <- File( "ReNprot-up.txt", parentId = "syn11629934" )
    annotations(f1) <- list( Type="Gene Set", Category="Interferome", Reference="LSP Experiment" )
    synStore(f1)
    f2 <- File( "ReNprot-dn.txt", parentId = "syn11629934" )
    annotations(f2) <- list( Type="Gene Set", Category="Interferome", Reference="LSP Experiment" )
    synStore(f2)
    
}

