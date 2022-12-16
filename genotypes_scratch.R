library(tidyverse)
library(markerGeneProfile)
library(ggpubr)
library(modelr)
library(MASS)
library(cowplot)
library(broom)

# get CMC clinical meta
cmc_clinical_path = '/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/Clinical/'
cmc_clinical_csv = paste0(cmc_clinical_path, 'CMC_Human_clinical_metadata.csv')

# get CMC rnaseq meta
cmc_rnaseq_path = '/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/Metadata/'
cmc_rnaseq_meta_path = paste0(cmc_rnaseq_path, 'CMC_Human_rnaSeq_metadata.csv')

cmc_rnaseq_meta = read_csv(cmc_rnaseq_meta_path)
colnames(cmc_rnaseq_meta) = make.names(colnames(cmc_rnaseq_meta), unique = T)

cmc_clinical_meta = read_csv(cmc_clinical_csv)
colnames(cmc_clinical_meta) = make.names(colnames(cmc_clinical_meta), unique = T)

# create a new df called cmc_meta that merges the clinical and rnaseq meta - note that this duplicates data from some subjects

cmc_meta = full_join(cmc_clinical_meta, cmc_rnaseq_meta)


cmc_meta$Dx = factor(cmc_meta$Dx, levels = c('Control', 'SCZ', 'BP', 'AFF', 'undetermined'))
cmc_meta$Reported.Gender = factor(cmc_meta$Reported.Gender, levels = c('Male', 'Female'))
cmc_meta$Ethnicity = factor(cmc_meta$Ethnicity, levels = c('Caucasian', 'African-American', 'Asian', 'Hispanic', '(Multiracial)'))

# convert ages to numbers
cmc_meta = cmc_meta %>% mutate(Age_norm = case_when(Age.of.Death == '90+' ~ 90, 
                                                    TRUE ~ as.numeric(Age.of.Death)))
# get numeric subject id - this is required for mapping to etienne's data from Pitt cohort
cmc_meta = cmc_meta %>% mutate(Subject = parse_number(SampleID))

# load in Alex's estimates of PRSs for SZ, BP, etc.
datit=read_csv('/external/rprshnas01/netdata_kcni/stlab/CMC_genotypes/SNPs/Release3/Metadata/CMC_Human_SNP_metadata.csv')

datit= datit %>% select(.,Individual_ID,Genotyping_Sample_ID)

prs=read_tsv("/external/rprshnas01/netdata_kcni/stlab/Alex_PRS/CMC_genotypes/PRSCS_cmcrodrigo")
prs$ID=gsub("0_", "", prs$ID)
datmrg1=merge(datit,prs,by.x = 'Genotyping_Sample_ID', by.y='ID')
# datcmc1=merge(datcmc,datmrg1,by.x="Individual.ID",by.y="Individual_ID")

joined_df = left_join(cmc_meta, datmrg1, by = c("Individual.ID" ="Individual_ID")) 


joined_df %>% ggplot(aes(x = Dx, y = SZ)) + geom_boxplot() + facet_wrap(~Institution)

# SST PRSs - based on Alex's estimates from CMC based on Fernanda's code

sst_prs = read_delim('/external/rprshnas01/kcni/asegura/SST_Fernanda/results/CMC_PRS_SST.txt', delim = '\t')
  
sst_prs$ID = str_sub(sst_prs$ID, 3)

joined_df = left_join(joined_df, sst_prs, by = c("Genotyping_Sample_ID" = "ID"))


# Now load in our data from CMC RNAseq data

### Now read in gene expression data from CommonMind subjects
gene_counts_ob = readRDS('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/geneCountsMerged.RDS')

# this is the expression data from Pitt Penn and MSSM
cmc_expr_1 = read_tsv('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/MSSM.Penn.Pitt_DLPFC.featureCount.tsv.gz')
# 
cmc_expr_2 = read_tsv('/external/rprshnas01/netdata_kcni/stlab/Public/commonmind/ControlledAccess/Data/RNAseq/Release4/QuantitatedExpression/NIMH.HBCC.featureCount.tsv.gz')
# 
cmc_expr = cbind(cmc_expr_1, cmc_expr_2[, 7:387]) 

cmc_expr_counts_only = cmc_expr[, 7:1011] %>% as.data.frame()
rownames(cmc_expr_counts_only) = cmc_expr$Geneid %>% substr(., 1, 15) %>% make.names(., unique = T)

### perform normalization for gene expression data from CMC

cmc_expr_cpm = edgeR::cpm(cmc_expr_counts_only, prior.count = 1)
cmc_expr_cpm = cmc_expr_cpm %>% as.data.frame() %>% tibble::rownames_to_column(var = "ensembl_id")

sst_mrna = cmc_expr_cpm %>% filter(ensembl_id == 'ENSG00000157005') %>% 
  select(-ensembl_id) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'SampleID')

colnames(sst_mrna) = c('SampleID', 'SST_mRNA')

# joined_df = left_join(joined_df, sst_mrna)



### GET MARKERS FOR MGP ANALYSIS
# note that this is the list of markers from micaela's paper - you get similar but diff results if you use the markers from the aging paper
sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

# I find it helpful to map some gene symbols to ensembl ids manually using mappings from hgnc, you can get those from here: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

hgnc_mapping = read_tsv('~/mathys_analysis/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list_ensembl <- lapply(new_cell_types, function(cell_type){
  return(new_markers %>% filter(subclass == cell_type, 
                                ensembl_id %in% rownames(cmc_expr_counts_only),
  ) %>% pull(ensembl_id))
})
names(new_marker_list_ensembl) <- new_cell_types
print(new_cell_types)

# new_marker_list_ensembl reflects the markers that we'll use in the MGP calculation


#### PERFORM MGP ANALYSIS

use_samples = cmc_meta %>% filter(!is.na(SampleID)) %>% pull(SampleID)

use_samples = intersect(use_samples, colnames(cmc_expr_cpm))

# perform MGP calculations, I set seekConsensus and removeMinority to FALSE as we're using Micaela's consensus markers
estimations <-  mgpEstimate(
  exprData= cmc_expr_cpm[, c('ensembl_id', use_samples)], # filters samples for those defined in use_samples
  genes=new_marker_list_ensembl,
  geneColName='ensembl_id',
  outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
  geneTransform = NULL, # this is the default option for geneTransform
  groups=NULL, #if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus=FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority=FALSE) 

# calculate qc metrics as well
pc1_exp = lapply(new_cell_types, function(cell_type){
  s = estimations$trimmedPCAs[[cell_type]] %>% summary()
  return(s$importance[3, 1])
}) %>% unlist()

mgp_qc_df = data.frame(subclass = new_cell_types, removedMarkerRatios = estimations$removedMarkerRatios, pc1_exp = pc1_exp )


mgp_df_big = estimations$estimates %>% as.data.frame() %>% tibble::rownames_to_column(var = "SampleID") 

# merge mgp df with metadata df
cmc_mgps_wide = left_join(mgp_df_big , cmc_meta) # this is the final mgp df 

prs_mgp_df = left_join(joined_df, cmc_mgps_wide %>% dplyr::select(SampleID, SST))

prs_mgp_df = left_join(prs_mgp_df, sst_mrna )


# 
# mgp_results = cmc_mgps_wide # calculated in mgp_scratch
# 
# joined_df = left_join(joined_df, mgp_results %>% select(SampleID:VLMC), by = "SampleID") %>% as.data.frame()

prs_mgp_df %>% ggplot(aes(x = Dx, y = SST)) + geom_boxplot()


### estimate beta coeffs for linear models including basic covariates and SST PRS
model_effects = prs_mgp_df %>% filter(!Institution == 'NIMH-HBCC', 
                     Dx %in% c('Control', 'SCZ','BP'), 
                     Sex %in% c('XY', 'XX'), Age_norm > 20) %>% 
  # group aging_paper_label stacked data by cell_type
  group_by(Institution) %>%
  
  # fit all the cell_type_prop data accorting to the model 
  # using the broom package to tidy the results 
  do(tidy(lm(scale(SST) ~ scale(Age_norm)  +  
               scale(PMI..in.hours.)  + scale(RIN) + 
               Sex + Dx + scale(SST_PRS_0.1),  data = .))) %>%
  
  # unstack the data and adjust for multiple comparisons using the Benjamini-Hochberg method
  ungroup() %>% 
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  
  # clean up the names the the term column
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept",
                       `SexM` = "gender:Male",
                       `DxSCZ` = "SCZ",
                       `DxBP` = "BP",
                       `DxMDD` = "MDD",
                       `EthnicityAfrican-American` = "ethnicity:Black", ## Q in the model ???
                       `scale(Age_norm)` = "age",
                       `scale(PMI..in.hours.)` = "PMI", 
                       `scale(RIN)` = "RIN",
                       `scale(SST_PRS_0.1)` = 'SST_PRS_0.1', 
                       `scale(SST_PRS_1)` = 'SST_PRS_1'))

model_effects %>% filter(term == 'SST_PRS_0.1') %>% 
  ggplot(aes(x = Institution, y = estimate)) + 
  geom_bar(stat = 'identity', width = 0.8, aes(y = estimate), colour = "black") +
  geom_errorbar((aes(ymin = estimate - std.error, ymax = estimate + std.error)))

# ggsave(filename = "sst_mgp_vs_prs_figure.pdf", device = "pdf", width = 8, height = 4, units = "in", dpi = 500, limitsize = F)

# print out the beta coeffs and association stats from linear model
model_effects %>% filter(term == 'SST_PRS_0.1')


# create a new dataframe with residuals after accounting for Age, PMI, RIN, Sex, and Dx per institution
inst_list = c('Penn', 'Pitt', 'MSSM')
prs_mgp_df_w_resids = lapply(inst_list, function(inst){
  
  SST_mgp_model = lm(scale(SST) ~ scale(Age_norm)  +  
                       scale(PMI..in.hours.) + scale(RIN) + 
                       Sex + Dx, data = prs_mgp_df  %>% filter(Institution == inst))
  
  new_df = prs_mgp_df %>% filter(Institution == inst) %>% 
    add_residuals(SST_mgp_model)
  return(new_df)
}) %>% bind_rows()


sst_resid_mgp_vs_prs_plot = prs_mgp_df_w_resids %>% 
  filter(!Institution == 'NIMH-HBCC') %>% 
  ggplot(aes(x = SST_PRS_0.1, y = resid)) + 
  geom_point(alpha = .5) + geom_smooth(method = 'lm', se = F) + 
  stat_cor(method = "pearson", size = 5) + 
  theme_cowplot() + 
  xlab('SST polygenic score') + 
  ylab('SST rel. cell proportions (AU)') + 
  facet_wrap(~Institution)

sst_resid_mgp_vs_prs_plot

ggsave(filename = "sst_mgp_vs_prs_figure.pdf", device = "pdf", width = 8, height = 4, units = "in", dpi = 500, limitsize = F)

sst_resid_mgp_vs_prs_all_insts_plot = prs_mgp_df_w_resids %>% 
  ggplot(aes(x = SST_PRS_0.1, y = resid)) + 
  geom_point(alpha = .5) + geom_smooth(method = 'lm', se = F) + 
  stat_cor(method = "pearson", size = 5) + 
  theme_cowplot() + 
  xlab('SST polygenic score') + 
  ylab('SST rel. cell proportions (AU)')

sst_resid_mgp_vs_prs_all_insts_plot

ggsave(filename = "sst_resid_mgp_vs_prs_all_insts_plot.pdf", device = "pdf", width = 4, height = 4, units = "in", dpi = 500, limitsize = F)


# 
# SST_mrna_model = lm(scale(log2(SST_mRNA+1)) ~ scale(Age_norm)  +  
#                       scale(PMI..in.hours.)  + scale(RIN) + 
#                       Sex + Dx + Institution, data = prs_mgp_df  %>% filter(!Institution == 'NIMH-HBCC'))
# 
# 
# sst_resid_mrna_vs_prs_plot =prs_mgp_df %>% filter(!Institution == 'NIMH-HBCC') %>% 
#   add_residuals(SST_mrna_model) %>% ggplot(aes(x = SST_PRS_0.1, y = resid, color = Dx, group = 1)) + 
#   geom_point(alpha = .5) + geom_smooth(method = 'lm') + 
#   stat_cor(method = "pearson", size = 5) + 
#   ylim(c(-2, 3)) + 
#   theme_cowplot() + 
#   xlab('SST PRS (0.1 threshold)') + 
#   ylab('SST cell proportions (demog. adj.)') + 
#   facet_wrap(~Institution)
# 
# prs_mgp_df %>% filter(!Institution == 'NIMH-HBCC') %>% 
#   add_residuals(SST_mgp_model) %>% ggplot(aes(x = Age_norm, y = SST, color = Dx, group = 1)) + 
#   geom_point(alpha = .5) + geom_smooth(method = 'lm') + 
#   stat_cor(method = "pearson", size = 5) + 
#   # ylim(c(-2, 3)) + 
#   theme_cowplot() + 
#   xlab('Age at death') + 
#   ylab('SST rel. cell proportions (AU)') + 
#   facet_wrap(~Institution)
# 
# prs_mgp_df %>% filter(!Institution == 'NIMH-HBCC') %>% 
#   add_residuals(SST_mgp_model) %>% ggplot(aes(x = Age_norm, y = SST, color = Dx, group = 1)) + 
#   geom_point(alpha = .25) + geom_smooth(method = 'lm') + 
#   stat_cor(method = "pearson", size = 5) + 
#   # ylim(c(-2, 3)) + 
#   theme_cowplot() + 
#   xlab('Age at death') + 
#   ylab('SST rel. cell proportions (AU)') + 
#   facet_wrap(~Institution)
# 
# 
# joined_df %>% ggplot(aes(x = SST_PRS_0.1, y = log2(SST+ 1))) + geom_point() + geom_smooth()  + facet_wrap(~Institution)
