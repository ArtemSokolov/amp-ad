#filter TAS table based on binding drug target for visualization function inputs
TAS_drug_binder<-function(TAS_tibble,drug_target,meta_tibble,condition){

  #set Drugs to lowercase
  TAS_tibble<-TAS_tibble %>% mutate_at( "name", str_to_lower )
  meta_tibble<-meta_tibble %>%  mutate_at( "Drug", str_to_lower )

  #select TAS data for drug
  cat('Selecting for Drugs that target:',drug_target,'\n')
  binders<-TAS_tibble %>% filter(symbol==drug_target & tas %in% c(1,2,3)) %>% pull(name) %>% unique()
  cat('Number of Drugs that DO bind:',length(binders),'\n')
  nonbinders<-TAS_tibble %>% filter(symbol==drug_target & tas == 10) %>% pull(name) %>% unique()
  cat('Number of Drugs that do NOT bind:',length(nonbinders),'\n')

  print('Filtering meta table for drugs')
  well_binders<- meta_tibble %>% filter(Drug %in% binders) %>% select(Well) %>% mutate( Binder = "Yes" )
  well_nonbinders<-meta_tibble %>% filter(Drug %in% nonbinders) %>% select(Well) %>% mutate( Binder = "No" )
  cat('Found',nrow(well_binders),'wells that match a binding Drug \n')
  cat('Found',nrow(well_nonbinders),'wells that match a nonbinding Drug \n')
  output<-bind_rows(well_binders,well_nonbinders)
  return(as.data.frame(output))
}

#plots TF targets from edge list in GSEA analysis
#plots TF targets from edge list in GSEA analysis using synapse for grabbing files
TF_enriched_targets_Violin_Plot<-function(count_file_syn_id='syn20494174',
                                          meta_file_syn_id='syn20495336',
                                          TAS_file_syn_id='syn18268627',
                                          folder='syn20555538',
                                          drug_target,
                                          TF_name){

  #load in counts, meta, TAS, DGE, and GSEA data for the particular drug_target and TF name
  print('Loading Data')
  counts<-syn_csv(count_file_syn_id)
  meta<-syn_csv(meta_file_syn_id)
  TAS<-syn_csv(TAS_file_syn_id)
  reference<-syn_csv(folder) #folder containing all DGE and GSEA synapse ids to load
  reference$genetarget <- str_to_lower(reference$genetarget) #make same case
  drug_target <- str_to_lower(drug_target) #make same case
  DGE<-reference %>% filter(file_type=='DGE' & genetarget == drug_target & fTox==FALSE & fConc==FALSE & mpi==TRUE & dge_method == 'EdgeR') %>%
    select(synid) %>% syn_tsv()
  GSEA<-reference %>% filter(file_type=='GSEA' & genetarget == drug_target & fTox==FALSE & fConc==FALSE & mpi==TRUE & dge_method == 'edger') %>%
    select(synid) %>% syn_tsv()

  #filter for samples that are binders or nonbinders
  meta<-TAS_drug_binder(TAS_tibble=TAS,drug_target=str_to_upper(drug_target),meta_tibble=meta) #filter for samples that are binders or nonbinders for drug target
  #filter for TF target
  TARGETS = unlist(strsplit(unlist(c(GSEA %>% filter(pathway==TF_name) %>% select(leadingEdge))),split = " "))#TF targets

  #visualization data  setup
  counts <- counts %>% filter(HUGO %in% TARGETS) #filter to targets
  tmp <- counts %>% filter(HUGO %in% TARGETS) #filter to targets
  genes<-counts$HUGO
  counts<-counts %>% select(meta$Well)#select samples
  samples <- colnames(counts)
  counts <- data.table::data.table(t(counts))
  colnames(counts)<-genes
  counts$Well <- samples
  meta <- meta[c("Well","Binder")]
  final <- merge(counts,meta,by="Well")
  final <- subset(final, select = -c(Well) )
  name<-names(final)
  final <- final %>% gather(Gene,Expression,name[1:length(name)-1])

  #plot
  p<-ggplot(final, aes_string(x="Gene", y="Expression", fill="Binder")) +
    geom_violin(trim=FALSE, color = NA, scale = "width") +
    ggtitle(TF_name) +
    xlab('TF Targets Indicating Enrichment') +
    theme(axis.text.x = element_text(angle = 90))

  #display the results
  #print(p)

  #DGE Results For Targets
  print('Target Genes Significance')
  print(DGE %>% filter(Gene %in% TARGETS))

  return(p)
}
