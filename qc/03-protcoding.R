## Functionality used by 03-protcoding.Rmd
##
## by Artem Sokolov

## Downloads a synapse ID to local disk and returns its filename
syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/QC" )@filePath }

## Custom ggplot theme that boldifies text elements
bold_theme <- function()
{
    etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
    theme_bw() + theme( axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          strip.text = etxt(12) )
}

## Loads all the background AUC score files
loadBKAUCs <- function()
{
    ## The IDs data frame is generated with the following code:
    ##
    ## source( "../R/resmine.R" )
    ## S <- allSettings() %>% filter( Task == "AC", Strategy != "strategy2" ) %>%
    ##     mutate( files = map( settings_md5, synBySettingsMd5 ) ) %>% select( -settings_md5 ) %>%
    ##     mutate( bkfile=map(files, filter, parentName=="score", name=="background_auc.csv") ) %>%
    ##     mutate( bkid = map(bkfile, select, id) ) %>% unnest(bkid) %>% select(-files, -bkfile)
    ## dput(S)
    IDs <- structure(list(Dataset = c("ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", 
                                      "ROSMAPpc", "ROSMAPpc", "ROSMAPpc", "ROSMAPpc"),
                          Linearity = c("Nonlinear", "Linear", "Linear", "Nonlinear",
                                        "Linear", "Nonlinear", "Nonlinear", "Linear"),
                          Strategy = c("strategy3", "strategy3", "strategy1", "strategy1",
                                       "strategy1", "strategy1", "strategy3", "strategy3"),
                          Task = c("AC", "AC", "AC", "AC", "AC", "AC", "AC", "AC"), 
                          id = c("syn12664727", "syn12664699", "syn12664695", "syn12664723", 
                                 "syn14712048", "syn14712492", "syn14715304", "syn14715369" )),
                     row.names = c(NA, -8L), class = c("tbl_df", "tbl", "data.frame"),
                     .Names = c("Dataset", "Linearity", "Strategy", "Task", "id"))

    ## Load individual files and aggregate their content
    IDs %>% mutate( fn = map_chr(id, syn) ) %>%
        mutate( X = map( fn, read_csv, col_types=cols() ) ) %>% select( -id, -fn ) %>%
        mutate( X2 = map( X, gather, Size, AUC ) ) %>% select( -X ) %>% unnest( X2 ) %>%
        mutate_at( "Size", as.integer )
}

## Loads all the relevant score files for Nienke's gene sets
loadScores <- function()
{
    ## The IDs data frame is generated with the following code:
    ##
    ## source( "../R/resmine.R" )
    ## S <- allSettings() %>% filter( Task == "AC", Strategy == "strategy3" ) %>%
    ##     mutate( files = map( settings_md5, synBySettingsMd5 ) ) %>% select( -settings_md5 ) %>%
    ##     mutate( gsfile = map(files, filter, parentName == "score", grepl("Nienke",name)) ) %>%
    ##     mutate( gsid = map(gsfile, select, id) ) %>% unnest( gsid ) %>%
    ##     select(-files, -gsfile, -Strategy, -Task)
    ## dput(S)
    IDs <- structure(list(Dataset = c("ROSMAP", "ROSMAP", "ROSMAPpc", "ROSMAPpc"),
                          Linearity = c("Nonlinear", "Linear", "Nonlinear", "Linear"),
                          id = c("syn12664731", "syn12664703", "syn14715307", "syn14715371")),
                     row.names = c(NA, -4L), class = c("tbl_df", "tbl", "data.frame"),
                     .Names = c("Dataset", "Linearity", "id"))

    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( URL = link, Name = name, Target = target_name )
    
    ## Load individual files and aggregate their content
    IDs %>% mutate( fn = map_chr(id,syn) ) %>%
        mutate( X = map(fn, read_csv, col_types=cols()) ) %>% select( -id, -fn ) %>%
        mutate( X2 = map(X, gather, URL, AUC) ) %>% select( -X ) %>% unnest(X2) %>%
        inner_join( MLINCS, by="URL" ) %>% select( -URL ) %>% spread( Dataset, AUC )
}

## Retrieves the vector of column names from ROSMAP
featsROSMAP <- function()
{
    syn("syn11738012") %>% read_tsv( col_types=cols(), n_max=5 ) %>% colnames
}

## Retrives the vector of column names from ROSMAPpc
featsROSMAPpc <- function()
{
    syn("syn14306482") %>% read_tsv( col_types=cols(), n_max=5 ) %>% colnames
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

## Computes the intersection of Nienke's gene sets against feature sets from both
##   ROSMAP and ROSMAPpc
getSetSizes <- function()
{
    ## Download all the relevant information
    vROS <- featsROSMAP()
    vRpc <- featsROSMAPpc()
    lGS <- syn("syn11973633") %>% read_gmt()
    
    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( URL = link, Name = name, Target = target_name )
    
    ## Compute intersection of gene sets with features from each dataset
    X <- data_frame( URL = names(lGS), Genes = unname(lGS) ) %>%
        inner_join( MLINCS, ., by="URL" ) %>% select( -URL ) %>%
        mutate( iROS = map(Genes, intersect, vROS),
               iRpc = map(Genes, intersect, vRpc) ) %>%
        mutate( nROS = map_int(iROS, length), nRpc = map_int(iRpc, length) )

    ## Identify the discrepancies
    Y <- X %>% filter( nROS != nRpc ) %>% select(-Genes) %>%
        mutate( eROS = map2(iROS, iRpc, setdiff), eRpc = map2(iRpc, iROS, setdiff) ) %>%
        mutate( ExtraGene = ifelse( nROS > nRpc, eROS, eRpc ),
               Dataset = ifelse( nROS > nRpc, "ROSMAP", "ROSMAPpc" ) ) %>%
        unnest(ExtraGene) %>% select( Name, Target, ExtraGene, Dataset )

    ## Compose a summary table
    Z <- Y %>% select( -Target ) %>% group_by( ExtraGene ) %>% arrange( Name ) %>%
        summarize( `Affected Drugs` = str_flatten( Name, ", " ), Dataset = unique(Dataset) ) %>%
        select( `Extra Gene` = ExtraGene, Dataset, everything() )
    Z
}

## A cache of getSetSizes() output created via dput()
getSetSizes_cache <- function()
{
structure(list(`Extra Gene` = c("AIM1", "BCRP1", "MAP3K20"), 
    Dataset = c("ROSMAP", "ROSMAP", "ROSMAPpc"), `Affected Drugs` = c("Alisertib, AT9283, Axitinib, Barasertib, BAY61-3606, BMS345541, BMS 777607, Cediranib, CHIR-99021, Crizotinib, Dactolisib, Erlotinib, Go 6976, GSK1070916, GW843682X, KIN001-111, KIN001-135, Linifanib, Linsitinib, LY2090314, MK 1775, MLN8054, Nintedanib, NVP-TAE226, NVP-TAE684, Palbociclib, PF04217903, PF562271, PHA-665752, PHA-767491, PKC412, PP1, R406, Rucaparib, Ruxolitinib, SB202190, SB 239063, Sotrastaurin, S-Ruxolitinib, Staurosporine, Sunitinib, TAK-632, TAK-715, URMC-099, Y-27632, ZM-447439", 
    "Gefitinib", "Bosutinib, Doramapimod, Imatinib, Ki20227, Masitinib, Motesanib, Neratinib, Nilotinib, RAF 265, SB 203580, Vandetanib"
    )), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, 
-3L), .Names = c("Extra Gene", "Dataset", "Affected Drugs"))
}

## A closer check of AZD 5438
check1 <- function()
{
    ## Load the mapping of LINCS URLs to drug names and their nominal targets
    MLINCS <- syn( "syn11801537" ) %>% read_csv( col_types=cols() ) %>%
        select( URL = link, Name = name, Target = target_name )
    lGS <- syn("syn11973633") %>% read_gmt()

    ## Identify genes associated with AZD 5438
    iURL <- MLINCS %>% filter( Name == "AZD 5438" ) %>% .$URL
    vGenes <- lGS[[iURL]]

    ## Download ROSMAP and ROSMAPpc
    ROS <- syn("syn11738012") %>% read_tsv( col_types=cols() )
    Rpc <- syn("syn14306482") %>% read_tsv( col_types=cols() )

    ## Compare the two datasets in the gene space defined by AZD 5438
    ROS1 <- ROS %>% arrange( ID ) %>% select( ID, AOD, Braak, one_of(vGenes) )
    Rpc1 <- Rpc %>% arrange( ID ) %>% select( ID, AOD, Braak, one_of(vGenes) )

    v1 <- as.numeric( filter( ROS1, AOD != "90+" )$AOD )
    v2 <- as.numeric( filter( Rpc1, AOD != "90+" )$AOD )
}

## A check of runs files for AZD 5438
check2 <- function()
{
    rROS <- syn("syn12663890") %>% read_csv( col_types = cols() )
    rRpc <- syn("syn14712042") %>% read_csv( col_types = cols() )

    rROSBK <- syn("syn12663896") %>% read_csv( col_types = cols() )
    rRpcBK <- syn("syn14714693") %>% read_csv( col_types = cols() )
}
