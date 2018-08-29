## QC plots for the first DGE experiment
##
## by Artem Sokolov

library( tidyverse )
library( synapser )
suppressMessages(synLogin())

## Retrieves a file from synapse to local disk and returns its local path
syn <- function( id, dlc = "~/data/AMP-AD/DGE" )
{ synGet( id, downloadLocation = dlc )$path }

## Makes a boxplot of the overall counts per well
plot_totalCounts <- function( CC )
{
    ## Generate a summary boxplot
    library( ggrepel )
    ggplot( CC, aes( x = "Master Plate", y = TotalCounts ) ) + theme_bw() +
        geom_boxplot( outlier.shape = NA ) +
        geom_point( data = subset(CC, TotalCounts < 1e5) ) +
        scale_y_continuous( trans="log10" ) + xlab("") + coord_flip() +
        geom_text_repel(data=subset(CC, TotalCounts < 1e4), aes(label=Well), segment.color="black",
                        nudge_x = 1, segment.alpha = 0.5, fontface="bold", color = "red" ) +
        geom_text_repel(data=subset(CC, TotalCounts < 1e4), aes(label=TotalCounts),
                        segment.color="black", segment.alpha = 0.5,
                        nudge_x = -1, fontface="bold", color = "red" ) +
        theme( axis.title = element_text( size=12, face="bold" ),
              axis.text = element_text( size=12, face="bold" ),
              axis.text.y = element_blank(), axis.ticks.y = element_blank() )
}

## Generates a heatmap for a pre-selected set of drugs
## X - genes-by-wells data.frame of counts
## Y - meta-data matrix with three columns: Well - Drug - Concentration
## NOTE: heatmap computed over genes that have non-zero values across all selected wells
hmapDrugs <- function( X, Y )
{
    ## Cluster all the stressors and controls
    vStress <- c( "LPS", "dsRNA +lipo", "naked dsRNA", "Lipo control", "Drug control" )
    
    ## Replaces 0s with NAs
    zero2na <- function( v ) { ifelse(v == 0, NA, v) }
    
    ## Identify the wells of interest
    Y1 <- Y %>% filter( Drug %in% vStress ) %>% select( Well, Drug )
    X1 <- X %>% select( HUGO, Y1$Well ) %>% mutate_all( zero2na ) %>% na.omit

    ## Compute pairwise correlations
    RR <- X1 %>% as.data.frame %>% column_to_rownames("HUGO") %>%
        cor( method="sp", use="pairwise.complete.obs" )
    RY <- Y1 %>% as.data.frame %>% mutate( Drug = as.factor(Drug) ) %>% column_to_rownames("Well")

    ## Set up a palette
    pal <- scales::hue_pal()(length(vStress)) %>% setNames( vStress ) %>% list( Drug = . )

    ## Plot the heatmap
    library( pheatmap )
    gg <- pheatmap( RR, annotation_row = RY, annotation_col = RY, annotation_colors = pal, silent=TRUE )

    ## Arrange all annotations to be in a single column
    gt <- gg$gtable
    i1 <- which( gt$layout$name == "legend" )
    i2 <- which( gt$layout$name == "annotation_legend" )
    gr <- gridExtra::arrangeGrob( gt$grobs[[i1]], gt$grobs[[i2]], heights=c(1,2) )
    gg <- gt[,1:4] %>% gtable::gtable_add_cols( gt$widths[6] ) %>%
        gtable::gtable_add_grob( gr, t = 4, b = 5, l = 5, clip="off" )
    gg
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

## Loads and composes a list of gene sets relevant for this QC vignette
## The result is stored on Synapse in syn15672121
composeGenesets <- function()
{
  ## Load all the relevant gene sets
  hset <- read_gmt( "h.all.v6.2.symbols.gmt" )
  vISG <- syn( "syn11629935" ) %>% scan( what=character() )
  vReNup <- syn( "syn12031302" ) %>% scan( what=character() )
  list( Schoggins_ISG = vISG, ReN_Proteomics = vReNup,
        IFN_alpha = hset$HALLMARK_INTERFERON_ALPHA_RESPONSE,
        IFN_gamma = hset$HALLMARK_INTERFERON_GAMMA_RESPONSE,
        IL2_STAT5 = hset$HALLMARK_IL2_STAT5_SIGNALING,
        DNA_repair = hset$HALLMARK_DNA_REPAIR,
        Mitotic_Spindle = hset$HALLMARK_MITOTIC_SPINDLE,
        Spermatogenesis = hset$HALLMARK_SPERMATOGENESIS,
        Myogenesis = hset$HALLMARK_MYOGENESIS
  )
  
  ## Write to a common .gmt file
  imap( GS, ~c(.y,.y,.x)) %>% map( str_flatten, "\t") %>% write_lines( "DGE-qc.gmt")
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

## Differential expression analysis (via edgeR)
## Pathway enrichment using fGSEA
main_edgeR <- function( X, Y )
{
  ## Compose the dichotomy of interest: (dsRNA & dsRNA+lipo) vs. Drug Control
  ## Remove L09 and O14 from the analysis, as they are outliers
  d1 <- c("naked dsRNA" = "dsRNA", "dsRNA +lipo" = "dsRNA", "Drug control" = "Control")
  Y1 <- Y %>% filter( Drug %in% names(d1), !(Well %in% c("L09", "O14")) ) %>%
    mutate( Treatment = d1[Drug] ) %>% select( -Concentration, -Drug ) %>%
    as.data.frame %>% column_to_rownames("Well")
  X1 <- X %>% select( HUGO, rownames(Y1) ) %>% as.data.frame %>% column_to_rownames( "HUGO" )

  ## Run bulk diff. exp. analysis
  library( edgeR )
  dl <- DGEList( counts = X1, samples = Y1 )
  dl <- calcNormFactors( dl )
  mmx <- model.matrix( ~Treatment, data = dl$samples )
  dl <- estimateDisp( dl, mmx )
  gf <- glmFit( dl, mmx ) %>% glmLRT( coef = 2 )
  RR <- topTags(gf, nrow(X1)) %>% as.data.frame %>% rownames_to_column( "Gene" )

  ## Check for enrichment of relevant gene sets
  library( fgsea )
  vv <- syn("syn15672202") %>% read_gmt()
  vz <- setNames( RR$logFC, RR$Gene )
  gRes <- fgsea( vv, vz, 10000 )
  plotGSEATable( vv, vz, gRes )
}

