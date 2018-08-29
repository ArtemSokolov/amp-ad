## Orphaned SCDE code
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/DGE" )@filePath }

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read.gmt <- function( fn, iName = 1 )
{
    ## Parse the raw sets
    v <- read_lines(fn) %>% as.list() %>%
        lapply( strsplit, "\\t" ) %>% lapply( "[[", 1 )
  
    ## Use the token iName as set names
    ## Drop the second token (which is a reference)
    names(v) <- lapply( v, "[[", iName ) %>% unlist()
    lapply( v, "[", -1:-2 )
}

## Loads and returns a list of all the relevant gene sets
get.genesets <- function()
{
    hset <- read.gmt( "genesets/h.all.v5.2.symbols.gmt" )
    vISG <- syn( "syn11629935" ) %>% scan( what=character() )
    vReNup <- syn( "syn12031302" ) %>% scan( what=character() )
    vReNdn <- syn( "syn12031304" ) %>% scan( what=character() )
    list( Schoggins_ISG = vISG, ReN_Proteomics = vReNup, ##ReNProt_dn = vReNdn,
##         GSE104190 = gs.GSE104190(),
         IFN_alpha = hset$HALLMARK_INTERFERON_ALPHA_RESPONSE,
         IFN_gamma = hset$HALLMARK_INTERFERON_GAMMA_RESPONSE,
         IL2_STAT5 = hset$HALLMARK_IL2_STAT5_SIGNALING,
         DNA_repair = hset$HALLMARK_DNA_REPAIR,
         Mitotic_Spindle = hset$HALLMARK_MITOTIC_SPINDLE,
         Spermatogenesis = hset$HALLMARK_SPERMATOGENESIS,
         Myogenesis = hset$HALLMARK_MYOGENESIS
         )
}

## Differential gene expression analysis
main.scde <- function()
{
    ## Load the data and relevant gene sets
    X <- syn( "syn11948496" ) %>% read_csv( col_types = cols() )
    Y <- syn( "syn11948497" ) %>% read_csv( col_types = cols() )

    ## Some conflicts with synapseClient()
    ## Load after pulling all the necessary data off Synapse
    library( scde )
    
    ## Compose the dichotomy of interest: dsRNA vs. Drug Control
    ## Remove L09 and O14 from the analysis, as they are outliers
    d1 <- c("naked dsRNA", "Drug control")
    D1 <- Y %>% filter( Drug %in% d1 ) %>%
        select( -Concentration ) %>% filter( !(Well %in% c("L09", "O14") ) )
    X1 <- X %>% select( HUGO, one_of( D1$Well ) ) %>% as.data.frame %>%
        column_to_rownames( "HUGO" ) %>% clean.counts
    y1 <- setNames( D1$Drug, D1$Well ) %>% factor( levels = d1 )

    ## Fit the models to individual samples
    ## Expensive operation
    ## Load results from file if available
    ##    ifm <- scde.error.models( counts = X1, groups = y1, n.cores = 7,
    ##                               save.model.plots=FALSE, verbose = 1 )
    ##    save( ifm, file="ifm-dsRNA-control.RData" )
    load( "qc/scde-dsRNA-control.RData" )

    ## Estimate gene expression prior
    ge.prior <- scde.expression.prior( models = ifm, counts = X1, length.out = 400, show.plot = TRUE )

    ## Run differential expression test
    ediff <- scde.expression.difference( ifm, X1, ge.prior, groups = y1,
                                        n.randomizations=100, n.cores = 1, verbose = 1 ) %>%
        rownames_to_column( "Gene" )
    ediff %>% write_csv( "scde-dsRNA-control.csv" )

    ## Check for enrichment of various gene sets
    ediff <- read_csv( "qc/scde-dsRNA-control.csv" )
    vz <- setNames( ediff$Z, ediff$Gene )
    vv <- get.genesets()
    gRes <- fgsea::fgsea( vv, vz, 10000 )
    gt <- plotGSEATable( vv, vz, gRes )
    ggsave( "scde-GSEA.png", gt, width = 5, height = 5 )
}

## Custom version of fgsea::plotGseaTable
plotGSEATable <- function(pathways, stats, fgseaRes,
                          gseaParam=1,
                          colwidths=c(2.2, 3, 0.8, 1, 1))
{
    library( grid )
    library( gridExtra )
    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))
    })

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- fgseaRes[match(pn, fgseaRes$pathway), ] %>%
            mutate( pp = as.character(round( pval, 3 )) ) %>%
            mutate( pp = ifelse( pp == "0", "< 0.001", pp ) )
        list(
            textGrob(pn, just="right", x=unit(0.95, "npc")),
            ggplot() +
            geom_segment(aes(x=p, xend=p,
                             y=0, yend=statsAdj[p]),
                         size=0.2) +
            scale_x_continuous(limits=c(0, length(statsAdj)),
                               expand=c(0, 0)) +
            scale_y_continuous(limits=c(-1, 1),
                               expand=c(0, 0)) +
            xlab(NULL) + ylab(NULL) +
            theme(panel.background = element_blank(),
                  axis.line=element_blank(),
                  axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  panel.grid = element_blank(),
                  axis.title=element_blank(),
                  plot.margin = rep(unit(0,"null"),4),
                  panel.spacing = rep(unit(0,"null"),4)
                  ),
            textGrob(sprintf("%d", annotation$size)),
            textGrob(sprintf("%.2f", annotation$NES)),
            textGrob(annotation$pp)
        )
    })
    
    
    rankPlot <-
        ggplot() +
        geom_blank() +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              plot.margin = unit(c(0,0,0.5,0), "npc"),
              panel.spacing = unit(c(0,0,0,0), "npc")
              )
    
    grid.arrange(grobs=c(
                     lapply(c("Pathway", "Gene ranks", "Size", "NES", "p-value"), textGrob),
                     unlist(ps, recursive = FALSE),
                     list(nullGrob(),
                          rankPlot
                          )),
                 ncol=5, widths=colwidths)
}
