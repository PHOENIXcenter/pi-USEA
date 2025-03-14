###USEA(Ubiquitin ligases/Deubiquitylases_Subtrates Enrichment Analysis)
###20250302-04
###From:LiuNing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(c("dplyr", "tidyr", "data.table"))
BiocManager::install(c("AnnotationDbi", "clusterProfiler"))
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

######################################################################################Part1
#interaction
##############################################E3REAL
setwd("file/data/part1_interaction")#Change your file address
E3Net <- read.csv("E3Net.csv", header = TRUE, sep = ",")
UbiNet <- read.csv("UbiNet.csv", header = TRUE, sep = ",")
E3S <- read.csv("E3_S.csv", header = TRUE, sep = ",")

common_columns <- Reduce(intersect, list(names(E3Net), names(UbiNet), names(E3S)))
E3Net_selected <- E3Net %>% select(common_columns)
UbiNet_selected <- UbiNet %>% select(common_columns)
E3S_selected <- E3S %>% select(common_columns)
combined_df <- bind_rows(E3Net_selected, UbiNet_selected, E3S_selected)
combined_df<-distinct(combined_df)#1108,189

combined_df_real <- combined_df %>% select(-contains("ORGANISM", ignore.case  = TRUE))
combined_df_real <- combined_df_real %>% 
  mutate(across(everything(), trimws)) %>%  
  distinct(ENZY_ACC_ID, SUB_ACC_ID, .keep_all = TRUE)#562,189

num_unique <- length(unique(combined_df_real$ENZY_ACC_ID))
print(num_unique)
write.csv(combined_df_real, "E3real.csv", row.names = FALSE)

##############################################E3PREDICT
#####REAL+predict
M.musculus_E3 <- read.csv("M.musculus_E3.csv", header = TRUE, sep = ",")
combinedpre <- bind_rows(M.musculus_E3[, c(1, 2)], combined_df)
combined_df_pre<-distinct(combinedpre)#1049683,248

combined_df_pre <- combined_df_pre %>% select(-contains("ORGANISM", ignore.case  = TRUE))
combined_df_pre <- combined_df_pre %>% 
  mutate(across(everything(), trimws)) %>%  
  distinct(ENZY_ACC_ID, SUB_ACC_ID, .keep_all = TRUE)#1049068,248

num_unique <- length(unique(combined_df_pre$ENZY_ACC_ID))
print(num_unique)

write.csv(combined_df_pre, "E3pre.csv", row.names = FALSE)
##############################################DUBREAL
DUBreal <- read.csv("DUB_S.csv", header = TRUE, sep = ",")
DUBreal <- DUBreal %>% 
  mutate(across(everything(), trimws)) %>%  
  distinct(ENZY_ACC_ID, SUB_ACC_ID, .keep_all = TRUE)#79,39

num_unique_DUBreal <- length(unique(DUBreal$ENZY_ACC_ID))
cat("The number of DUBreal:", num_unique_DUBreal, "\n")

write.csv(DUBreal, "DUBreal.csv", row.names = FALSE)
##############################################DUBPREDICT
#####REAL+predict
DUBpre <- read.csv("M.musculus_DUB.csv", header = TRUE, sep = ",")
transdsi_dub<- read.csv("TransDSI.csv", header = TRUE, sep = ",")#19461
transdsi_dub_exploded <- transdsi_dub %>%
  separate_rows(`Contributing.Residues..SUB.`, sep = ",")

transdsi <- transdsi_dub_exploded %>%
  rename(ENZY_ACC_ID = SwissProt.ID..DUB.,
         SUB_ACC_ID = SwissProt.ID..SUB.,
  )
transdsi <- transdsi_dub_exploded %>%
  rename(ENZY_ACC_ID = SwissProt.ID..DUB.,
         SUB_ACC_ID = SwissProt.ID..SUB.,
         Ubiquitination_site = Contributing.Residues..SUB.)

transdsi$Ubiquitination_site <- gsub("\\[|\\]", "", transdsi$Ubiquitination_site)

transdsi_dub_selected <- transdsi %>%
  select(ENZY_ACC_ID, SUB_ACC_ID)
DUBpre_selected <- DUBpre %>%
  select(ENZY_ACC_ID, SUB_ACC_ID)

DUBpre<- bind_rows(transdsi_dub_selected, DUBpre_selected)

combinedpre <- bind_rows(DUBpre[, c(1, 2)], DUBreal)
combined_df_pre<-distinct(combinedpre)

combined_df_pre <- combined_df_pre %>% select(-contains("ORGANISM", ignore.case  = TRUE))
combined_df_pre <- combined_df_pre %>% 
  mutate(across(everything(), trimws)) %>%  
  distinct(ENZY_ACC_ID, SUB_ACC_ID, .keep_all = TRUE)#1081570,169

num_unique <- length(unique(combined_df_pre$ENZY_ACC_ID))
print(num_unique)

combined_df_pre <- combined_df_pre %>% select(where(~ !all(is.na(.))))
write.csv(combined_df_pre, "DUBpre.csv", row.names = FALSE)


######################################################################################Part2
#Site
####################################################################PSP,IUUCD,PLMD，transdsi
setwd("file/data/part2_site")#Change your file address
PSP_site <- read.csv("PSP_site.csv", header = TRUE, sep = ",")#（from PSP）26341
IUUCD<- read.csv("IUUCD.csv", header = TRUE, sep = ",")#IUUCD 538
Ubisite <- read.csv("Ubisite.csv", header = TRUE, sep = ",") #5004

df_processed <- IUUCD %>%
  separate_rows(Ubiquitination_site, sep = ";") %>%
  mutate(Ubiquitination_site = paste0("K", Ubiquitination_site)) %>%
  distinct()

PSP_site_extracted <- PSP_site %>%
  select(ENZYME, SUB_MOD_RSD) %>%
  rename(Gene_name = ENZYME, Ubiquitination_site = SUB_MOD_RSD)

df_combined <- bind_rows(df_processed, PSP_site_extracted)

SUB_GENE <- bitr(Ubisite$Uniprot.ID, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Mm.eg.db") %>%
  distinct()
#SUB_GENE <- bitr(Ubisite$Uniprot.ID, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") %>%
#  distinct()
Ubisite <- Ubisite %>%
  rename(UNIPROT = Uniprot.ID) %>%
  distinct()

merged_Ubisite <- left_join(Ubisite, SUB_GENE, by = "UNIPROT")

merged_Ubisite_selected <- merged_Ubisite %>%
  select(SYMBOL, Sites) %>%
  rename(Gene_name = SYMBOL, Ubiquitination_site = Sites) %>%
  mutate(Ubiquitination_site = paste0("K", Ubiquitination_site),
         Gene_name = paste(toupper(substring(Gene_name, 1, 1)), tolower(substring(Gene_name, 2)), sep = ""))

df_final <- bind_rows(df_combined, merged_Ubisite_selected) %>%
  mutate(across(everything(), trimws)) %>%  
  distinct(Gene_name, Ubiquitination_site, .keep_all = TRUE)  
transdsi_selected <- transdsi %>%
  select(Gene_name = SUB_ACC_ID, Ubiquitination_site)
combined_data <- bind_rows(df_final , transdsi_selected)

Site_unique <- combined_data %>%
  mutate(across(everything(), trimws)) %>%  
  distinct(Gene_name, Ubiquitination_site, .keep_all = TRUE)

unique_combinations_count <- Site_unique %>%
  distinct(Gene_name, Ubiquitination_site) %>%
  nrow()
cat("The number of combinations:", unique_combinations_count, "\n")#134793
write.csv(Site_unique, "Site_unique.csv", row.names = FALSE)


######################################################################################Part3
#calculation
#####################################################################################
setwd("file/data/part3_calculation")#Change your file address
exp_data<- read.csv("test_data.csv", header = TRUE, sep = ",")
#Change exp_data to your own experimental data

combined_df_pre <- read.csv("E3pre.csv", header = TRUE, sep = ",")
combined_df_real <- read.csv("E3real.csv", header = TRUE, sep = ",")
DUBpre <- read.csv("DUBpre.csv", header = TRUE, sep = ",")
DUBreal <- read.csv("DUBreal.csv", header = TRUE, sep = ",")
site <- read.csv("Site_unique.csv", header = TRUE, sep = ",")

###########################################################################
num_unique_combined_df_pre <- length(unique(combined_df_pre$ENZY_ACC_ID))
num_unique_combined_df_real <- length(unique(combined_df_real$ENZY_ACC_ID))
num_unique_DUBpre <- length(unique(DUBpre$ENZY_ACC_ID))
num_unique_DUBreal <- length(unique(DUBreal$ENZY_ACC_ID))

 
cat("The number of E3pre 中不同 ENZY_ACC_ID 的种类数:", num_unique_combined_df_pre, "\n")
cat("The number of E3real:", num_unique_combined_df_real, "\n")
cat("The number of DUBpre:", num_unique_DUBpre, "\n")
cat("The number of DUBreal:", num_unique_DUBreal, "\n")

unique_combinations_count <- site %>%
  distinct(Gene_name, Ubiquitination_site) %>%
  nrow()
cat("The number of combinations:", unique_combinations_count, "\n")

common_DUB <- intersect(DUBreal$ENZY_ACC_ID, DUBpre$ENZY_ACC_ID)
num_common_DUB <- length(common_DUB)
cat("The number of common:", num_common_DUB, "\n")

common_combined <- intersect(combined_df_pre$ENZY_ACC_ID, combined_df_real$ENZY_ACC_ID)
num_common_combined <- length(common_combined)
cat("The number of common:", num_common_combined, "\n")




######################################################################
pre_subgene <- bitr(combined_df_pre$SUB_ACC_ID, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
#pre_subgene <- bitr(combined_df_pre$SUB_ACC_ID, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
merged <- merge(x = combined_df_pre, y = pre_subgene ,by.x = "SUB_ACC_ID", by.y = "UNIPROT", all = TRUE)	##和上面效果一致
merged$SUB_GENE[is.na(merged$SUB_GENE)] <- merged$SYMBOL[is.na(merged$SUB_GENE)]
######################################################################E3pre+site
dfpre_gene<- merge(x =merged, y = site, by.x = "SYMBOL", by.y = "Gene_name", all = TRUE)
result_pre_gene <- merge(dfpre_gene, exp_data, by.x = c("SUB_ACC_ID", "Ubiquitination_site"), by.y = c("SUB_ACC_ID", "SUB_MOD_RSD"))

ENZYME_stats <- result_pre_gene %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

print(ENZYME_stats)
df <- left_join(result_pre_gene, ENZYME_stats, by = "ENZY_ACC_ID")
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))
df_summary <- df %>%
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE)
  )

print(df_summary)

m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd

df$z_score <- (df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2 * (1 - pnorm(abs(df$z_score)))  

df$fdr <- p.adjust(df$p_value, method = "fdr") 


df <- df %>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

print(head(df))
write.csv(df, "E3pre_score.csv", row.names = TRUE)
#next############################################################E3withoutsite 
result<- merge(x =merged, y = exp_data, by.x = "SUB_ACC_ID", by.y = "SUB_ACC_ID")

ENZYME_stats <- result %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

print(ENZYME_stats)
df <- left_join(result, ENZYME_stats, by = "ENZY_ACC_ID")

df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

print(df_summary)

m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd

df$z_score <- (df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  # 计算 FDR
df<- df %>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

write.csv(df, "E3pre_nosite.csv", row.names = FALSE)##without site


################################################################################E3real
dfreal_gene<- merge(x =combined_df_real, y = site, by.x = "SUB_GENE", by.y = "Gene_name", all = TRUE)
#32276-71646
result_real_gene <- merge(dfreal_gene, exp_data, by.x = c("SUB_ACC_ID", "Ubiquitination_site"), by.y = c("SUB_ACC_ID", "SUB_MOD_RSD"))
#8————3049
#####################
ENZYME_stats <- result_real_gene %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result_real_gene, ENZYME_stats, by = "ENZY_ACC_ID")

df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- (df$logFC_mean - m) * sqrt(df$logFC_count) / s
df$p_value <- 2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  # 计算 FDR
df <- distinct(df)
df<- df %>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)
write.csv(df, "E3real_score.csv", row.names = TRUE)

#next############################################################withoutsite 
result<- merge(x =combined_df_real, y = exp_data, by.x = "SUB_ACC_ID", by.y = "SUB_ACC_ID")
#### 
ENZYME_stats <- result %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result, ENZYME_stats, by = "ENZY_ACC_ID")
 
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- (df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  # 计算 FDR
df<- df %>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)
write.csv(df, "E3real_nosite.csv", row.names = FALSE)

#####################################################################################DUBpre
pre_subgene <- bitr(DUBpre$SUB_ACC_ID, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Mm.eg.db)
#pre_subgene <- bitr(DUBpre$SUB_ACC_ID, fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
merged <- merge(x = DUBpre, y = pre_subgene ,by.x = "SUB_ACC_ID", by.y = "UNIPROT", all = TRUE)	
merged$SUB_GENE[is.na(merged$SUB_GENE)] <- merged$SYMBOL[is.na(merged$SUB_GENE)]
##########################################################################DUBpre+site
dfpre_gene<- merge(x =merged, y = site, by.x = "SYMBOL", by.y = "Gene_name", all = TRUE)
result_pre_gene <- merge(dfpre_gene, exp_data, by.x = c("SUB_ACC_ID", "Ubiquitination_site"), by.y = c("SUB_ACC_ID", "SUB_MOD_RSD"))

ENZYME_stats <- result_pre_gene %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result_pre_gene, ENZYME_stats, by = "ENZY_ACC_ID")
 
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- -(df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2*(1 - pnorm(abs(df$z_score)))

df$fdr <- p.adjust(df$p_value, method = "fdr")  # 计算 FDR
df <- distinct(df)

num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)
write.csv(df, "DUBpre_score.csv", row.names = TRUE)



#next############################################################withoutsite 
result<- merge(x =DUBpre, y = exp_data, by.x = "SUB_ACC_ID", by.y = "SUB_ACC_ID")

#### 
ENZYME_stats <- result %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result, ENZYME_stats, by = "ENZY_ACC_ID")
 
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- -(df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  #

df<- df%>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

write.csv(df, "DUBpre_nosite.csv", row.names = FALSE)##without site


##################################################################################DUBreal
dfreal_gene<- merge(x =DUBreal, y = site, by.x = "SUB_GENE", by.y = "Gene_name", all = TRUE)
#32276-59974 
result_real_gene <- merge(dfreal_gene, exp_data, by.x = c("SUB_ACC_ID", "Ubiquitination_site"), by.y = c("SUB_ACC_ID", "SUB_MOD_RSD"))
#8————232
#####################
ENZYME_stats <- result_real_gene %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result_real_gene, ENZYME_stats, by = "ENZY_ACC_ID")
 
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- -(df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <-2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  
df <- distinct(df)

df<- df%>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

write.csv(df, "DUBreal_score.csv", row.names = TRUE)

#next############################################################withoutsite DUBreal
result<- merge(x =DUBreal, y = exp_data, by.x = "SUB_ACC_ID", by.y = "SUB_ACC_ID")

#### 
ENZYME_stats <- result %>%
  group_by(ENZY_ACC_ID) %>%
  summarise(logFC_mean = mean(FC), logFC_count = n())

 
print(ENZYME_stats)
df <- left_join(result, ENZYME_stats, by = "ENZY_ACC_ID")
 
df_summary <- df %>% 
  summarise(
    logFC_mean_avg = mean(logFC_mean, na.rm = TRUE),
    logFC_mean_sd = sd(logFC_mean, na.rm = TRUE))

 
print(df_summary)
 
m <- df_summary$logFC_mean_avg
s <- df_summary$logFC_mean_sd
 
df$z_score <- -(df$logFC_mean - m) * sqrt(df$logFC_count) / s

df$p_value <- 2*(1 - pnorm(abs(df$z_score)))
df$fdr <- p.adjust(df$p_value, method = "fdr")  
df<- df%>% select(where(~ !all(is.na(.))))
num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

num_unique <- length(unique(df$ENZY_ACC_ID))
print(num_unique)

write.csv(df, "DUBreal_nosite.csv", row.names = FALSE)##without site
