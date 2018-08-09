## Plots for the 2018-05-25 meeting
##
## by Artem Sokolov

source( "../../R/freezemine.R" )

## Prefixes a filename with the directory corresponding to the analyses in this file
pfn <- function( fn ) {str_c( "2018-05-25/", fn )}

## Bold-face text element of requested size
## Used extensively in ggplot themes throughout the file
etxt <- function(s, ...) {element_text(size = s, face="bold", ...)}

## Plots performance of background sets on ROSMAP across various task definitions
bkMain <- function()
{
    ## Download all relevant background files
    BB <- bkSummary() %>% filter( grepl("ROSMAP", Dataset) ) %>% extract2( "SynID" ) %>%
        map( ~synBkAUC(.x, TRUE) ) %>% bind_rows

    ## Define common theme elements
    thm <- theme( axis.text = etxt(12), axis.title = etxt(14), strip.text = etxt(14),
                 legend.text = etxt(12), legend.title = etxt(14) )

    ## Plot all results, showing how performance of background sets varies from
    ##   task to task
    gg <- ggplot( BB, aes( x=Size, y=AUC, color=Task ) ) + theme_bw() + thm +
        geom_smooth( lwd = 1.5 ) +
        guides( color=guide_legend(override.aes=list(fill="white")) )
    ggsave( "bkROSMAP.png", gg, width=8, height=5 )
}

## Summary plot of performance for a set of drugs
perfSummary <- function()
{
    ## Identify drugs for which in vitro data is available
    XIV <- synStats( "syn12081315" )

    ## Match up drugs against the nominal target names
    TGT <- synGet( "syn11801537", downloadLocation = "~/data/AMP-AD" )@filePath %>%
        read_csv %>% select( Name = name, Target = target_name, URL = link )
    
    ## Load performance on Nienke's set
    XX <- statsSummary() %>% filter( grepl("Nienke", Filename), Dataset=="ROSMAP" ) %$%
        setNames( SynID, Filename ) %>% map( ~synStats(.x, LINCSmap(), TRUE) ) %>%
        bind_rows( .id = "Filename" ) %>% select( -Filename ) %>%
        mutate( Profiled = str_to_lower(Name) %in% XIV$Name )

    ## For each profiled drug, identify its best-performing task
    X <- XX %>% filter(Profiled) %>% arrange(pval) %>% group_by(Name) %>%
        filter(row_number() == 1 ) %>% ungroup %>% filter( pval < 0.05 ) %>%
        inner_join( TGT )

    ## Download all relevant background files
    ## Match up the set sizes for each drug in XX against background
    BB <- bkSummary() %>% filter( grepl("ROSMAP", Dataset) ) %>% extract2( "SynID" ) %>%
        map( ~synBkAUC(.x, TRUE) ) %>% bind_rows
    BK <- X %>% select( Name, Size, Task ) %>% mutate( Size = round(Size/10)*10) %>%
      inner_join( BB ) %>% rename( Drug = Name )
        
    ## Compute p-values
    ## f.pval <- function( x, drug, task ) 
    ## {
    ##   v <- filter(BK, Name==drug, Task==task)$AUC
    ##   sum( v > x ) / length(v)
    ## }
    ## YY <- XX %>% rowwise %>% mutate( pp = f.pval(AUC, Name, Task) ) %>% ungroup

    ## Hand-fix missing Target entries
    N2T <- set_names( X$Target, X$Name )
    N2T["Alisertib"] <- "AURKA"
    N2T["SB202190"] <- "p38 MAPK"
    X <- X %>% mutate( Target = N2T[Name] )
    
    ## Reorder drugs based on AUC performance
    X <- X %>% arrange(AUC) %>% mutate( Name = factor(Name,Name) ) %>%
      mutate( Lbl = str_c(" ", Target) )
    BK <- BK %>% mutate( Drug = factor(Drug, levels(X$Name)))

    library( ggridges )
    gg <- ggplot(BK, aes(x=AUC, y=Drug, fill=Task)) + 
      theme_ridges(center_axis_labels=TRUE) +
      geom_density_ridges2(scale=1.2, size=1, alpha=0.5) +
      geom_segment( data=X, aes(x=AUC, xend=AUC, y=as.numeric(Name), yend=as.numeric(Name)+0.9 ), color="red", lwd=2 ) +
      geom_text( data=X, aes(x=AUC, y=as.numeric(Name) + 0.5, label=Lbl), hjust=0, fontface="bold", size=5 ) +
      scale_fill_manual( values=c("tomato","darkolivegreen","steelblue")) +
      scale_y_discrete(expand=c(0,0.95)) +
      theme( axis.text = etxt(12), axis.title = etxt(14),
             legend.position=c(0,0), legend.justification =c(1,0.5),
             legend.title=element_text(face="bold"), legend.text=element_text(face="bold"),
             legend.box.background = element_rect(color="gray40") )
    ggsave( "candidates.png", gg, width=7, height=7 )
}

## Wrangles viability tables into a single results data frame
wrangleViability <- function()
{
    ## Load individual files
    fns <- list.files(pattern="csv")
    X <- set_names( fns, str_sub(fns,1,-5) ) %>% map( read_csv, col_type=cols() )
    
    ## Combine everything into a common data frame
    XX <- imap( X, ~mutate(.x, Drug = .y)) %>% bind_rows %>% 
      gather( Concentration, Viability, -Drug ) %>% na.omit()
    
    write_csv( XX, "all-viability.csv" )
}

## Plots viability data
plotViability <- function()
{
    ## Define the order of concentration values to plot along the x axis
    ordConc <- c( "Control", "0", "0.01", "0.1", "1", "10" )
  
    ## Load viability data
    X <- read_csv( "all-viability.csv", col_types=cols() ) %>%
      mutate( Concentration = factor(Concentration, ordConc) )
    levels(X$Concentration)[1] <- "DMSO"

    ## Define a theme common to multiple plots
    thm <- theme( axis.text = etxt(12), axis.title = etxt(14),
                  strip.text = etxt(14), strip.background = element_blank() )
        
    ## Plot everything
    gg <- ggplot( X, aes(x = Concentration, y = Viability) ) + theme_bw() +
      stat_summary( fun.y = "mean", geom="bar", alpha=0.75, fill="steelblue" ) +
      geom_jitter( width=0.05 ) + facet_wrap( ~Drug ) + thm
    ggsave( "viability.png", gg, width = 10, height=8)
    
    ## Plot a "zoom in" onto a single drug
    X1 <- X %>% filter( Drug == "Metformin")
    gg1 <- ggplot( X1, aes(x = Concentration, y = Viability) ) + theme_bw() +
      stat_summary( fun.y = "mean", geom="bar", alpha=0.75, fill="steelblue" ) +
      geom_jitter( width=0.05 ) + facet_wrap( ~Drug ) + thm
    ggsave( "viability-Metformin.png", gg1, width=5, height=5 )
}
