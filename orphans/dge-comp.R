## Orphaned enrichment of DGE differential expression
##   for Nienke's gene sets
##
## by Artem Sokolov

## Check for enrichment of Nienke's sets in our differential expression signatures
diffexp.enrich <- function()
{
    ## Load the relevant data
    load( "edgeR-allDrugs.RData" )
    NN <- synGet( "syn11973633", downloadLocation = "~/data/AMP-AD" )@filePath %>% read.gmt
    MM <- synGet( "syn11801537", downloadLocation = "~/data/AMP-AD" )@filePath %>% read_csv

    ## Rename the gene sets from links to drugs
    MM <- MM %>% filter( link %in% names(NN) )
    stopifnot( all(MM$link == names(NN)) )
    names(NN) <- tolower( MM$name )

    ## Reduce DGE-derived differential expression profiles and
    ##  Nienke's sets to the same set of drugs
    vNames <- intersect( names(RR), names(NN) )
    gRes <- list()
    gtl <- list()
    for( i in vNames )
    {
        cat( "Working with", i, "\n" )
        vz <- setNames( RR[[i]]$logFC, RR[[i]]$Gene )
        gRes[[i]] <- fgsea::fgsea( NN[i], vz, 10000 )
        gtl[[i]] <- plotGSEATable( NN[i], vz, gRes[[i]] )
    }

    ## Reorder the drugs from most to least significant enrichment
    vOrder <- arrange( bind_rows( gRes ), pval )$pathway
    gtl <- gtl[vOrder]

    ## Concatenate all the enrichment plots into a single table
    ## Keep the header from the first result
    ## Keep the footer from the last result
    ## Use the top 10 hits
    gtt <- lapply( gtl[2:9], function(.) {.[2,]} ) %>%
        c( list(gtl[[1]][1:2,]), ., list(gtl[[10]][2:3,]) )
    gt <- arrangeGrob( grobs=gtt, ncol=1 )
    
    ggsave( "drug-gsea.png", gt, width = 7, height = 7 )
}

