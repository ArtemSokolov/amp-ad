## Differential expression and gene set enrichment
##   for Metformin drug vs. control comparison
##
## by Artem Sokolov

library( tidyverse )
library( synapser )

synLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD", ifcollision = "overwrite.local" )$path }

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
    read_lines(fn) %>% str_split( "\\t" ) %>%
        set_names( map_chr(., nth, iName) ) %>%
        map( ~.x[-2:-1] )
}

main <- function()
{
    ## Retrieve differential expression signature for Metformin
    RR <- syn( "syn15674107" ) %>% read_csv( col_types=cols() ) %>% filter( Drug == "metformin" )
    RR %>% write_tsv( "metformin/dfexp-metformin.tsv" )
    vsig <- with( RR, set_names( logFC, Gene ) )

    ## Load curated gene sets
    lGS <- syn( "syn16551164" ) %>% read_gmt()

    ## Run FGSEA
    gRes <- fgsea::fgsea( lGS, vsig, 1e5 )

    ## Identify the top 10 "up" and top 10 "down" pathways
    Rup <- gRes %>% filter( NES>0 ) %>% arrange( pval ) %>% head(10)
    Rdn <- gRes %>% filter( NES<0 ) %>% arrange( pval ) %>% head(11) %>% slice( -10 )
    
    ## Plot the top k pathways
    gtup <- plotGSEATable( lGS[Rup$pathway], vsig, Rup, 1, c(5, 1.5, 0.4,0.4,0.4) )
    gtdn <- plotGSEATable( lGS[Rdn$pathway], vsig, Rdn, 1, c(5, 1.5, 0.4,0.4,0.4) )
    gt <- arrangeGrob( gtup, gtdn, ncol=1 )
    ggsave( "metformin-gsea.png", gt, width=13.5, height=7 )
}

## Custom version of fgsea::plotGseaTable
plotGSEATable <- function(pathways, stats, fgseaRes, gseaParam=1,
                          colwidths=c(0.8, 1.5, 0.4, 0.4, 0.4))
{
  library( grid )
  library( gridExtra )

  ## Scale the signature
  statsAdj <- sort( stats, decreasing=TRUE )
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  ## Match gene names to their indices in the scale signature
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  
  ## Generates a plot of segments associated with a position vector p
  fsegs <- function(p) {
    ggplot() + xlab(NULL) + ylab(NULL) +
      geom_segment(aes(x=p, xend=p, y=0, yend=statsAdj[p]), size=1) +
      scale_x_continuous(limits=c(0, length(statsAdj)), expand=c(0, 0)) +
      scale_y_continuous(limits=c(-1, 1), expand=c(0, 0)) +
      theme(panel.background = element_blank(), axis.line=element_blank(),
            axis.text=element_blank(), axis.ticks=element_blank(),
            panel.grid = element_blank(), axis.title=element_blank(),
            plot.margin = rep(unit(0,"null"),4),
            panel.spacing = rep(unit(0,"null"),4))  
  }
    
  ## Reformat the input results matrix
  RF <- as_data_frame(fgseaRes) %>% select( pathway, size, pval, NES ) %>%
    mutate( clr = ifelse(NES > 0, "red", "blue") ) %>%
    mutate_at( "pval", round, 3 ) %>% 
    mutate_at( "NES", ~sprintf("%.2f",.x) ) %>%
    mutate_all( as.character ) %>%
    mutate( pval = ifelse(pval == "0", "< 0.001", pval) ) %>%
    mutate( ii = pathways[pathway] )
  
  ## Generate grobs for each pathway
  GRBS <- RF %>% mutate_at("clr", map, ~gpar(col=.x, fontface="bold")) %>%
    mutate_at( vars(size, pval, NES), map, textGrob, gp=gpar(fontface="bold") ) %>%
    mutate_at( "ii", map, fsegs ) %>%
    mutate( pw = map2(pathway, clr, ~textGrob(.x, just="right", x=unit(0.95,"npc"), gp=.y)) )

  ## Convert grobs table to a list
  grbs <- GRBS %>% mutate( .rn = 1:n() ) %>% split( .$.rn ) %>%
    map( select, pw, ii, size, NES, pval ) %>% map( flatten )
  
  ## Axis for the rank plot  
  axisPlot <- ggplot() + geom_blank() + xlab(NULL) + ylab(NULL) +
    scale_x_continuous(limits=c(0, length(statsAdj)), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-1, 1), expand=c(0, 0)) +
    theme(panel.background = element_blank(), axis.line=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold"),
          panel.grid = element_blank(), axis.title=element_blank(),
          plot.margin = unit(c(0,0,0.5,0), "npc"),
          panel.spacing = unit(c(0,0,0,0), "npc"))
  
  grid.arrange(grobs=c(
    lapply(c("Pathway", "Gene ranks", "Size", "NES", "p-value"),
           textGrob, gp=gpar(fontface="bold")),
    flatten(grbs), list(nullGrob(), axisPlot)),
    ncol=5, widths=colwidths)
}
