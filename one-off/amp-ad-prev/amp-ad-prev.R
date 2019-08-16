## Gene sets found by previously-published AMP-AD studies
##
## by Artem Sokolov

library( tidyverse )

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
  read_lines(fn) %>% str_split( "\\t" ) %>%
    set_names( map_chr(., nth, iName) ) %>%
    map( ~.x[-2:-1] )
}

## Genes present in AMP-AD datasets
ampad <- function()
{
    synapser::synLogin()
    synExtra::synDownloader(".", ifcollision="overwrite.local")("syn16204450") %>%
                                                                  read_gmt()
}

## PMID: 30297968
## Title: Integrative transcriptome analyses of the aging brain
##        implicate altered splicing in Alzheimer’s disease susceptibility
## Table: List of differentially spliced introns associated with clinical AD
##        status in ROS/MAP that replicate in the MSBB dataset
## Index: S8
## Keyphrase: intron splicing
gs1 <- function()
{
    rmt <- file.path( "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0238-1",
                     "MediaObjects/41588_2018_238_MOESM7_ESM.xlsx")
    openxlsx::read.xlsx( rmt ) %>% pull( gene_id ) %>% unique()
}

## PMID: 29598827
## Title: Integrated biology approach reveals molecular and pathological interactions
##        among Alzheimer’s Aβ42, Tau, TREM2, and TYROBP in Drosophila models
## Table: Overlap between HBTRC human AD co-expression network modules and DEGs
##        identified in Tau/TREM2WT/TYROBP and Tau/TREM2R47H/TYROBP files.
## Index: S8a
## Keyphrase: Tau/TREM2-R47H/TYROBP model
gs2 <- function()
{
    rmt <- file.path( "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-018-0530-9",
                     "MediaObjects/13073_2018_530_MOESM9_ESM.xlsx")
    openxlsx::read.xlsx( rmt, sheet=2, startRow=3 ) %>% pluck( "Overlap.HumanOrtholog", 1 ) %>%
        str_split( ";" ) %>% pluck(1)
}

## PMID: 30136084
## Title: Divergent brain gene expression patterns associate with distinct
##        cell-specific tau neuropathology traits in progressive supranuclear palsy
## Index: S1
## Keyphrase: tau neuropathology
gs3 <- function()
{
    rmt <- file.path( "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6208732",
                     "bin/401_2018_1900_MOESM1_ESM.xlsx")
    openxlsx::read.xlsx( rmt, sheet=2 ) %>% filter( !is.na(Gene) ) %>% pull(Gene)
}

## PMID: 29358399
## Title: Tau induces blood vessel abnormalities and angiogenesis-related
##        gene expression in P301L transgenic mice and human Alzheimer’s disease
## Table: Differential expression between subjects in B3 stages (Braak V/VI)
##        and B1 stages (Braak 0/I/II)
## Index: S02
## Keyphrase: angiogenesis
gs4 <- function()
{
    rmt <- file.path( "https://www.pnas.org/highwire/filestream/790703",
                     "field_highwire_adjunct_files/1/pnas.1710329115.sd02.xlsx")
    openxlsx::read.xlsx( rmt, sheet=2 ) %>% pull( Gene ) %>% unique
}

## PMID: 29107053
## Title: Conserved brain myelination networks are altered in Alzheimer's
##        and other neurodegenerative diseases
## Table: Combined differential gene expression (DGE) and co-expression
##        network analyses results for AD vs. non-AD subjects, temporal cortex,
##        "Comprehensive Model"
## Index: S3
## Keyphrase: differential gene expression
gs5 <- function()
{
    rmt <- "https://ars.els-cdn.com/content/image/1-s2.0-S1552526017337664-mmc2.xlsx"
    openxlsx::read.xlsx( rmt, sheet=4 ) %>% filter( Dx.Bonferroni.PValue < 0.05 ) %>%
        pull( gene_symbol ) %>% unique()
}

## PMID: 29110684
## Title: Multiscale network modeling of oligodendrocytes reveals molecular
##        components of myelin dysregulation in Alzheimer’s disease
## Table: Differentially expressed genes found the mouse knockout RNAseq
##        analyses of Ugt8 in the FC (Data 2)
## Index: S5/D2
## Keyphrase: Ugt8 knock-out mouse
gs6 <- function()
{
    rmt <- file.path( "https://static-content.springer.com/esm/art%3A10.1186%2Fs13024-017-0219-3",
                     "MediaObjects/13024_2017_219_MOESM5_ESM.zip")
    lcl <- "PMID-29110684.zip"
    download.file( rmt, lcl )
    zip::unzip( lcl, "Data 2.tsv" )
    read_tsv( "Data 2.tsv", col_types=cols() ) %>%
        filter( Q.Value < 0.05, !is.na(HGCN_Symbol) ) %>%
        pull(HGCN_Symbol) %>% unique()
}

## PMID: 27989508
## Title: A Multi-network Approach Identifies Protein-Specific Co-expression
##        in Asymptomatic and Symptomatic Alzheimer’s Disease
## Table: RNA and Co-expression Analysis and Module Membership
## Index: S8/M14
## Keyphrase: tubulin-binding proteins
gs7 <- function()
{
    synapser::synLogin()
    lcl <- synExtra::synDownloader(".", ifcollision="overwrite.local")("syn19001714")
    read_csv( lcl, col_types=cols() ) %>%
        filter( `net$colors` == "T-M14" ) %>%
        pull(GI) %>% unique
}

## PMID: 27799057
## Title: Integrative network analysis of nineteen brain regions identifies
##        molecular signatures and networks underlying selective regional
##        vulnerability to Alzheimer’s disease
## Table: Genes correlated with cognitive/neuropathological traits.
## Index: S7
## Keyphrase: cognitive traits
gs8 <- function()
{
    rmt <- file.path( "https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-016-0355-3",
                     "MediaObjects/13073_2016_355_MOESM1_ESM.xlsx")
    openxlsx::read.xlsx( rmt, startRow=3, sheet=9 ) %>%
        filter( FDR == 0, !is.na(GeneSymbol) ) %>%
        pull(GeneSymbol) %>% unique() %>%
        str_split( " /// " ) %>% unlist
}

## PMID: 29802388
## Title: A molecular network of the aging human brain provides insights
##        into the pathology and cognitive decline of Alzheimer’s disease
## Table: Table contains the gene to module assignments. 47 modules with
##        at least 20 genes were identified using the SpeakEasy algorithm.
## Index: S3/m109
## Keyphrase: Cognitive Decline
gs9 <- function()
{
    rmt <- file.path("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6599633",
                     "bin/NIHMS962042-supplement-Sup_Table_3.xlsx")
    openxlsx::read.xlsx( rmt, startRow=5 ) %>%
        filter( Module.ID == "m109" ) %>% pull( Gene.Symbol )
}

main <- function()
{
    ## Download raw files from the corresponding papers
    ## Annotate each filename with PMID and Suppl. Table index
    X <- tribble(
        ~PMID, ~Table, ~Parser, ~Keyphrase,
        "30297968", "S8",      gs1, "Intron Splicing",
        "29598827", "S8a",     gs2, "Tau/TREM2-R47H/TYROBP",
        "30136084", "S1",      gs3, "Tau Neuropathology",
        "29358399", "S02",     gs4, "Angiogenesis",
        "29107053", "S3",      gs5, "Differential Gene Expression",
        "29110684", "S5/D2",   gs6, "Ugt8 Knock-out Mouse",
        "27989508", "S8/M14",  gs7, "Tubulin-binding Proteins",
        "27799057", "S7",      gs8, "Cognitive traits",
        "29802388", "S3/m109", gs9, "Cognitive decline"
    )

    ## Retrieve the gene sets
    GS <- X %>% mutate( GeneSet = map(Parser, exec) )

    ## Compose a .gmt file
    GMT <- GS %>% select( -Parser ) %>% mutate_at( "GeneSet", map_chr, str_flatten, "\t" ) %>%
        mutate( Content = glue::glue("PMID:{PMID}-{Table}\t{Keyphrase}\t{GeneSet}") )
    cat( GMT$Content, sep="\n", file="prev-ampad.gmt" )
}
