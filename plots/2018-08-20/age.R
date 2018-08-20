## Distributions across Braak staging
##
## by Artem Sokolov

source( "api.R" )
library( ggforce )
library( ggthemes )

synapseLogin()

ageByBraak <- function()
{
    ## Load ROSMAP data
    X <- syn( "syn14306482" ) %>% read_tsv(col_types = cols(ID = col_character()))

    ## Braak to Stage mapping associated with the AC task
    b2s <- c( "A", "A", "A", NA, NA, "C", "C" )
    
    ## Isolate AOD and Braak
    Y <- X %>% select( AOD, Braak ) %>%
        mutate( Age = map_if(AOD, ~(. == "90+"), ~"90") ) %>% unnest %>%
        mutate_at( "Age", as.double ) %>%
        mutate( Stage = b2s[Braak+1] ) %>%
        mutate_at( "Braak", factor )

    ## Make plots
    gg1 <- ggplot( Y, aes(x=Braak, y=Age, color=Stage) ) + theme_bw() + bold_theme() + geom_sina() +
        scale_color_manual( values=c("A"="steelblue", "C"="tomato"), na.value="gray40", guide=FALSE )

    gg2 <- ggplot( na.omit(Y), aes(x=Stage, y=Age, fill=Stage) ) + theme_bw() + bold_theme() +
        geom_violin( alpha=0.8 ) + geom_sina( size = 0.5 ) +
        scale_fill_manual( values=c("A"="steelblue", "C"="tomato"), guide=FALSE )

    gridExtra::arrangeGrob( gg1, gg2, nrow=1, widths=c(3,1) ) %>%
        ggsave( "plots/age-by-braak.png", ., width=8, height=4 )
}

AODconf <- function()
{
    ## A couple of mappings for renaming
    m1 <- c( "All-pairs" = "All Pairs", "AOD-one" = "Min Age", "ID-reuse" = "ID reuse" )
    m2 <- c( "Braak" = "Braak", "AOD" = "Age" )
    
    ## Load cached AUC estimates
    AA <- syn( "syn14716566" ) %>% read_csv( col_types=cols() ) %>% select(-Label) %>%
        filter( Estimate %in% c( "All-pairs", "AOD-one", "ID-reuse" ) ) %>%
        mutate( Estimate = m1[Estimate], Labels = m2[Labels] )

    ## Define the palette and the main plotting function
    pal <- few_pal()(6)
    fplot <- function( A )
    {
        ggplot( A, aes( x=Size, y=AUC, color=Estimate ) ) +
            theme_bw() + geom_smooth( se=FALSE ) + bold_theme() +
            facet_wrap( ~Labels, nrow=1 ) +
            scale_color_manual( values=pal[c(1,3,4,5)] )
    }
    
    ## Create a series of summary plots
    AA %>% filter(Estimate=="All Pairs") %>% fplot %>% ggsave("plots/AOD1.png", ., width=8, height=3.5)
    AA %>% filter(Estimate != "Min Age") %>% fplot %>% ggsave("plots/AOD2.png", ., width=8, height=3.5)
    AA %>% fplot %>% ggsave( "plots/AOD3.png", ., width=8, height=3.5 )
}
