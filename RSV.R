library(tidyverse)
library(readxl)
library(sqldf)
library(data.table)
library(HGNChelper)
library(STRINGdb)
library(RandomWalkRestartMH) 
library(gprofiler2)
library(RColorBrewer)
library(xlsx)

### Load data
data <- read_delim("./processed_data/vep_slivar_annovar.vcf.hg38_multianno.txt",
                   col_names = TRUE, na = ".")
control <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/control_AF.txt", 
                      col_names = c("CHR", "POS", "REF", "ALT", "AF"))
#revel <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/revel_with_transcript_ids", delim = ",", na = ".") %>%
#  mutate(chr=str_c("chr", chr, sep = ""))
#gnomad_constraint <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/gnomad.v2.1.1.lof_metrics.by_gene.txt") %>% 
#  select(gene, pLI, oe_lof_upper) 
gnomad_constraint <- read_delim("~/Documents/R_projects/datasets/gnomAD_constraint/20201106_gnomAD.LOEUF.pLI.annotations.tsv.gz")
dbnsfp <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/dbNSFP4.1_gene", na = ".")
IUI_and_PanelApp_PID <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/IUI_and_PanelApp_PID.txt", delim = "\n",
                                   col_names = "hgnc_symbol")
RSV_genes_disgenet <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/RSV_disgenet.tsv")
RSV_genes_open_target_v6 <- read_delim("~/Documents/R_projects/PRI_downstream_vep_gnomad3/RSV_genes_open_targets_v6.tsv")

### Fill all NA CADD scores by uploading to https://cadd.gs.washington.edu/score
# Rename columns of data
data <- data %>% 
  rename(CHR=Otherinfo4, POS=Otherinfo5, ID=Otherinfo6, REF=Otherinfo7, ALT=Otherinfo8, gnomAD_AF=AF)

cadd_na <- data  %>%          #Subset file to thise with NA patho score
  filter(is.na(REVEL_score) & is.na(CADD_phred)) %>%
  select(CHR:ALT) %>%
  mutate(ID = ".") %>%
  rename("#CHR" = CHR)
write_tsv(cadd_na, file = "./processed_data/cadd_na.tsv", col_names = T) #Submit this file to the CADD website
cadd_na_results <- read_delim("./processed_data/cadd_na_results.tsv", skip = 2,
                                    col_names = c("CHR",	"POS",	"REF",	"ALT",	"RawScore",	"PHRED")) %>% 
  mutate(CHR = str_c("chr", CHR))
#replace NA CADDs with downloaded values
data <-  data %>% 
  left_join(cadd_na_results) %>% 
  mutate(PHRED = coalesce(PHRED, CADD_phred)) #Fill NA values in PHRED with CADD_phred values

### Filter on MAF
data_gnomad <- data %>% 
  filter(gnomAD_AF<=0.01 | is.na(gnomAD_AF)) %>%
  filter(ExAC_ALL<=0.01 | is.na(ExAC_ALL)) %>% 
  filter(ALL.sites.2015_08<=0.01 | is.na(ALL.sites.2015_08)) %>% 
  filter(esp6500siv2_all<=0.01 | is.na(esp6500siv2_all))

### compare AF of case vs control
data_gnomad_control <- data_gnomad %>% 
  rowwise() %>% 
  mutate(case_af = str_split(Otherinfo11, ";")[[1]][2]) %>%         #Extract case AF
  mutate(case_af = as.numeric(str_sub(case_af, start = 4))) %>%     
  left_join(control, by = c("CHR", "POS", "REF", "ALT")) %>%        #left_join with control
  rename(control_af = AF) %>%                                     #Rename to control_af
  replace_na(list(control_af = 0))                                  #replace NAs with 0 in control_af

data_gnomad_control <- data_gnomad_control %>% 
  filter(case_af >= control_af) %>% #Compare AF case vs control
  filter((case_af >= gnomAD_AF) | is.na(gnomAD_AF)) %>% 
  filter((case_af >= ExAC_ALL) | is.na(ExAC_ALL)) %>% 
  filter((case_af >= ALL.sites.2015_08) | is.na(ALL.sites.2015_08)) %>% 
  filter((case_af >= esp6500siv2_all) | is.na(esp6500siv2_all))

### Remove MT variant
data_gnomad_control_MT <- data_gnomad_control %>% filter(Chr != "chrM")

### Filter on REVEL and CADD
data_gnomad_control_MT_patho <- data_gnomad_control_MT %>% 
  filter(REVEL_score>=0.5 | (is.na(REVEL_score) & PHRED>=20))

### At least 1 homo OR haploinsufficient
# count number of hetero and homo
data_gnomad_control_MT_patho$n_hom <- 0
data_gnomad_control_MT_patho$n_het <- 0
for (i in 1:nrow(data_gnomad_control_MT_patho)) {
  for(j in 13:139){
    column_idx <- str_c("Otherinfo", j)
    genotype <- str_split(data_gnomad_control_MT_patho[i, column_idx], ":")[[1]][1]
    n_dots <- str_count(genotype, "\\.")
    n_ones <- str_count(genotype, "1")
    n_zeros <- str_count(genotype, "0")
    if(n_dots >= 1) {next}
    if(n_ones==1 && n_zeros==1) {
      data_gnomad_control_MT_patho$n_het[i] <- data_gnomad_control_MT_patho$n_het[i]+1
    }
    if(n_ones==2 && n_zeros==0){
      data_gnomad_control_MT_patho$n_hom[i] <- data_gnomad_control_MT_patho$n_hom[i]+1
    }
  }
}
# add hgnc symbols for left_join
new_hgnc_table <- getCurrentHumanMap()
hgnc_symbols <- checkGeneSymbols(data_gnomad_control_MT_patho$Gene.ensGene,
                                 map=new_hgnc_table, unmapped.as.na=T)
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="PSEN2;RP5-1087E8.6"] <- "PSEN2"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="LY75;LY75-CD302"] <- "LY75"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="BBS5;RP11-724O16.1"] <- "BBS5"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="PP7080"] <- "SLC9A3"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="AC008937.3"] <- "SETD9"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="PLPBP;ADGRA2"] <- "PRSS1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="OXR1;OXR1"] <- "OXR1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="MUC2;MUC5AC"] <- "MUC2"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="HSPB2;HSPB2-C11orf52"] <- "HSPB2"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="AP002884.2;SDHD"] <- "SDHD"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="ESR2;SYNE2"] <- "SYNE2"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="CORO7;CORO7-PAM16"] <- "CORO7"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="ACSM1;THUMPD1"] <- "ACSM1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="MC1R;RP11-566K11.2"] <- "TUBB3"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="P2RX5;P2RX5-TAX1BP3"] <- "P2RX5-TAX1BP3"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP11-192H23.4;RSKR"] <- "KIAA0100"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="P2RY11;PPAN-P2RY11"] <- "P2RY11"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="L3MBTL1;RP1-138B7.8"] <- "L3MBTL1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="CTA-246H3.13;CTA-246H3.12"] <- "LRP5L"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP5-972B16.2;TSPAN7"] <- "TSPAN7"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="MUTYH;RP4-534D1.6"] <- "MUTYH"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="DCTN1;RP11-287D1.3"] <- "DCTN1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="PRSS50;RP11-316F12.1"] <- "PRSS50"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="HULC;RNU6ATAC21P"] <- "HULC"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="NONE;LINC00839"] <- "LINC00839"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP11-307N16.6;SPATA13"] <- "SPATA13"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="CHMP4A;RP11-468E2.1"] <- "CHMP4A"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP11-343C2.11;VPS4A"] <- "VPS4A"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="DSG1;DSG1-AS1"] <- "DSG1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="ALDH16A1;CTD-3148I10.9"] <- "ALDH16A1"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP11-115F18.1;RNU6-713P"] <- "RNU6-713P"
hgnc_symbols$Suggested.Symbol[hgnc_symbols$x=="RP11-468E2.1;RP11-468E2.2;TM9SF1"] <- "TM9SF1"
# Add official gene symbols
data_gnomad_control_MT_patho$hgnc_symbol <- hgnc_symbols$Suggested.Symbol
# Add information from gnomAD and dbnsfp
data_gnomad_control_MT_patho <- data_gnomad_control_MT_patho %>% 
  left_join(gnomad_constraint, by=c("hgnc_symbol"="genes")) %>% 
  left_join(dbnsfp, by=c("hgnc_symbol"="Gene_name"))
# filter on pli and homozygous samples
genes_to_keep <- data_gnomad_control_MT_patho %>%
  filter(n_hom>=1 | GnomAD_pLI>=0.9 | is.na(GnomAD_pLI)) %>% 
  select(hgnc_symbol)
data_gnomad_control_MT_patho_hompli <- data_gnomad_control_MT_patho %>% 
  filter(hgnc_symbol %in% genes_to_keep$hgnc_symbol)

### Random Walk
#Retrieve StringDB
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=400, input_directory="")
string_id_mapped <- string_db$map(as.data.frame(data_gnomad_control_MT_patho_hompli),
                                  "hgnc_symbol", removeUnmappedRows = F ) %>% 
  select(Gene.ensGene, Gene.refGene, hgnc_symbol, STRING_id) %>% 
  distinct(hgnc_symbol, .keep_all = T)
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="CFAP276"] <- "9606.ENSP00000358965"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="GARIN4"] <- "9606.ENSP00000294829"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="GBA3"] <- "9606.ENSP00000314508"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="ZBED10P"] <- "9606.ENSP00000343242"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="C8ORF34-AS1"] <- "9606.ENSP00000427820"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="P2RX5-TAX1BP3"] <- "9606.ENSP00000225328"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="IGLV11-55"] <- "9606.ENSP00000431254"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="IGLV9-49"] <- "9606.ENSP00000431254"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="IGLV3-19"] <- "9606.ENSP00000431254"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="IGLV3-12"] <- "9606.ENSP00000431254"
string_id_mapped$STRING_id[string_id_mapped$hgnc_symbol=="PPP4R3C"] <- "9606.ENSP00000483228"
undefined <- string_id_mapped %>% filter(is.na(STRING_id))
string_id_mapped <- string_id_mapped %>% filter(!is.na(STRING_id))
#Read networks and clean them.
string_db_graph <- string_db$get_graph()
string_db_graph <- igraph::simplify(string_db_graph, remove.multiple = TRUE, remove.loops = TRUE) 
STRING_Multiplex <- create.multiplex(list(PPI=string_db_graph))
AdjMatrix_STRING_PATH <- compute.adjacency.matrix(STRING_Multiplex)
AdjMatrixNorm_STRING_PATH <- normalize.multiplex.adjacency(AdjMatrix_STRING_PATH)
#Prepare seeds
RSV_hgnc_open_target_v6 <- checkGeneSymbols(RSV_genes_open_target_v6$symbol, map = new_hgnc_table) %>% filter(!is.na(Suggested.Symbol))
RSV_hgnc_disgenet <- checkGeneSymbols(RSV_genes_disgenet$Gene, map = new_hgnc_table) %>% filter(!is.na(Suggested.Symbol))
RSV_seed_gene_union <- data.frame(hgnc_symbol=union(RSV_hgnc_disgenet$Suggested.Symbol, RSV_hgnc_open_target_v6$Suggested.Symbol) )
RSV_seed_gene_string_id_union <- string_db$map(RSV_seed_gene_union, "hgnc_symbol", removeUnmappedRows = T)$STRING_id
#Random walk on RSV seeds
RWR_STRING_Results_union <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_STRING_PATH,
                                                            STRING_Multiplex,RSV_seed_gene_string_id_union)
RSV_RWR_STRING_Results_Scores_union <- left_join(string_id_mapped,
                                                   RWR_STRING_Results_union$RWRM_Results, by=c("STRING_id"="NodeNames")) %>% 
  rename("Score_RSV_union"="Score")
RSV_RWR_STRING_Results_NAs_union <- RSV_RWR_STRING_Results_Scores_union %>% filter(is.na(Score_RSV_union ))
RSV_RWR_STRING_Results_notNAs_union  <- RSV_RWR_STRING_Results_Scores_union  %>% filter(!is.na(Score_RSV_union ))
counts_RSV_union <- data.frame(STRING_id=RSV_RWR_STRING_Results_notNAs_union$STRING_id,
                                 Score_RSV_union=RSV_RWR_STRING_Results_notNAs_union$Score_RSV_union,
                                 rank_RSV_union=rank(-RSV_RWR_STRING_Results_notNAs_union$Score_RSV_union),
                                 count_score=rep(0, nrow(RSV_RWR_STRING_Results_notNAs_union)), 
                                 count_rank=rep(0, nrow(RSV_RWR_STRING_Results_notNAs_union)))
#Random walk on random seeds
set.seed(12345)
n_permute <- 1000
genes_index_on_STRING <- which(igraph::V(string_db_graph)$name %in% RSV_RWR_STRING_Results_notNAs_union$STRING_id)
whole_ix_exclude_candidates <- seq(1, length(igraph::V(string_db_graph)$name))[-genes_index_on_STRING]
for (i in 1:n_permute) {
  #generate a random set with the size of original seeds (here 392 genes) 
  random_number <- sample(whole_ix_exclude_candidates, length(RSV_seed_gene_string_id_union), replace = F)
  temp_seed_node <- igraph::V(string_db_graph)[random_number]$name
  RWR_STRING_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_STRING_PATH,
                                                      STRING_Multiplex,temp_seed_node)
  temp_RWR_STRING_Results <- left_join(counts_RSV_union,
                                       RWR_STRING_Results$RWRM_Results, by=c("STRING_id"="NodeNames"))
  #If the rank is less than RSV_related rank, increase count
  must_increase_ids_rank <- temp_RWR_STRING_Results %>%
    filter(!is.na(Score)) %>%
    mutate(rank_random=rank(-Score)) %>% 
    filter(rank_random<=rank_RSV_union) %>%
    select(STRING_id)
  counts_RSV_union[counts_RSV_union$STRING_id %in% must_increase_ids_rank$STRING_id, "count_rank"] = counts_RSV_union[counts_RSV_union$STRING_id %in% must_increase_ids_rank$STRING_id, "count_rank"] +1
  #If the Score is more than RSV_related Score, increase count
  must_increase_ids_score <- temp_RWR_STRING_Results %>%
    filter(!is.na(Score)) %>%
    filter(Score>=Score_RSV_union) %>%
    select(STRING_id)
  counts_RSV_union[counts_RSV_union$STRING_id %in% must_increase_ids_score$STRING_id, "count_score"] = counts_RSV_union[counts_RSV_union$STRING_id %in% must_increase_ids_score$STRING_id, "count_score"] +1
}
#Calculate p-values
pvalues_RSV_union <- counts_RSV_union %>% mutate(pvalue_rank=count_rank/n_permute) %>% mutate(pvalue_score=count_score/n_permute)
pvalues_RSV_union$hgnc_symbol <- RSV_RWR_STRING_Results_notNAs_union$hgnc_symbol
pvalues_RSV_union$pvalue_adjusted_rank <- p.adjust(pvalues_RSV_union$pvalue_rank, method = "fdr")
pvalues_RSV_union$pvalue_adjusted_score <- p.adjust(pvalues_RSV_union$pvalue_score, method = "fdr")
data_gnomad_control_MT_patho_hompli_rw <- data_gnomad_control_MT_patho_hompli %>%
  filter(hgnc_symbol %in% pvalues_RSV_union$hgnc_symbol[pvalues_RSV_union$pvalue_adjusted_score <= 0.05]) %>% 
  left_join(pvalues_RSV_union %>% select(hgnc_symbol, Score_RSV_union))
#Extract useful data for adjacent genes
adjacent_genes_info <- data_gnomad_control_MT_patho_hompli_rw %>%
  select(Chr:Alt, hgnc_symbol, case_af, control_af, gnomAD_AF, n_het, n_hom, GnomAD_pLI, Otherinfo11)
#Extract useful data for seed genes
seed_genes_info <- data_gnomad_control_MT_patho_hompli %>%
  filter(Gene.ensGene %in% RSV_RWR_STRING_Results_NAs_union$Gene.ensGene) %>% 
  select(Chr:Alt, hgnc_symbol, case_af, control_af, gnomAD_AF, n_het, n_hom, GnomAD_pLI, Otherinfo11)

### LoF
HC_indexes <- c()
LC_indexes <- c()
for(i in 1:nrow(data_gnomad_control_MT_patho_hompli)){
  if(!is.na(str_match(string=data_gnomad_control_MT_patho_hompli$Otherinfo11[i], "\\|HC\\|"))[1,1]){
    HC_indexes <- c(i, HC_indexes)
  }
  if(!is.na(str_match(string=data_gnomad_control_MT_patho_hompli$Otherinfo11[i], "\\|LC\\|"))[1,1]){
    LC_indexes <- c(i, LC_indexes)
  }
}
lof <- c(HC_indexes, LC_indexes)
data_gnomad_control_MT_patho_hompli_lof <- data_gnomad_control_MT_patho_hompli[lof, ]

### Manually look at variants with NA LOEUF
data_gnomad_control_MT_patho_cadd_homloeuf_lof_loeuf <- data_gnomad_control_MT_patho_cadd_homloeuf_lof %>% 
  filter(!is.na(oe_lof_upper) | (hgnc_symbol %in% c("GARIN4", "IFNA5", "OR14I1", "P2RX5-TAX1BP3", "PKD1L3")))

### PID
data_gnomad_control_MT_patho_cadd_homloeuf_pid <- data_gnomad_control_MT_patho_cadd_homloeuf %>% 
  filter(hgnc_symbol %in% IUI_and_PanelApp_PID$hgnc_symbol)


