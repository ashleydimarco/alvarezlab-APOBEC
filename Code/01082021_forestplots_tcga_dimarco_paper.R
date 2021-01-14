#################################################
# Sarah C. Van Alsten                           #
# Date: January 8, 2021                         #
# Institution: UNC Chapel Hill                  #
# Purpose: Make Forest Plots for RFDs in TCGA   #
# For Ashley Di Marco APOBEC Paper              #
# Packages: tidyverse, janitor                  #
# Last Updated: January 8, 2021                 #
#################################################


# Open Libraries ----------------------------------------------------------

library(tidyverse)




# Data Management ---------------------------------------------------------

#inputting of the data

tcga <- readxl::read_excel("data/11032020_TCGA_Repair_Score_sample_data_PCA.xlsx", sheet = 2) %>%
  janitor::clean_names()

#for tcga, not all clinical characteristics yet accounted for read in those files to merge
y <- read.table("data/brca_merged_with_abc.txt", header = 1)

tcga <- tcga %>%
  full_join(y, by = c("sample" = "bcr_patient_barcode.1"))




#recoding variables before running models
tcga$breast_carcinoma_estrogen_receptor_status <- factor(tcga$breast_carcinoma_estrogen_receptor_status, levels=c("Positive","Negative"))
tcga$er_pos <- tcga$breast_carcinoma_estrogen_receptor_status
tcga$p53_RNA <- factor(tcga$p53_RNA, levels=c("WT-like","Mut-like"))
tcga$pam50 <- factor(tcga$PAM50, levels=c("LumA","LumB", "Normal", "Her2", "Basal"))
tcga$pathologic_stage_AW <- factor(tcga$pathologic_stage_AW, levels=c("Stage I","Stage II", "Stage III", "Stage IV", "Stage X"))
tcga$race2 <- factor(ifelse(tcga$race == "BLACK OR AFRICAN AMERICAN", "Black", 
                            ifelse(tcga$race == "WHITE", "White", NA)))
tcga$race2 <- relevel(tcga$race2, ref = "White")
tcga$race3 <- ifelse(tcga$race == "BLACK OR AFRICAN AMERICAN", "AA", 
                     ifelse(tcga$race == "WHITE" | tcga$race == "ASIAN" | tcga$race == "AMERICAN INDIAN OR ALASKA NATIVE", "Non-AA", NA))
tcga$menopause <- factor(ifelse(tcga$menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre",
                                ifelse(tcga$menopause_status == "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", "Post", NA)))
tcga$node <- factor(ifelse(tcga$number_of_lymphnodes_positive_by_he == 0, "Negative",
                           ifelse(tcga$number_of_lymphnodes_positive_by_he != "[Not Available]", "Positive", NA)))
tcga$her2_status <- factor(ifelse(tcga$her2_immunohistochemistry_level_result %in% c( "2+", "3+"), "Positive", "Negative"))
tcga$dna_group_cat <- factor(tcga$dna_group, levels = c("Repair High", "APOBEC", "Heterogeneous Repair", "HR/FA"))

#apobec group (high vs low)
di <- readxl::read_excel("data/dimarco/DiMarco-APOBEC-Immune-Signature-Data_original.xlsm")
di <- janitor::clean_names(di)
#Per Andrea's email: "The enrichment score of 2 is used as a cutoff for "APOBEC high" and "APOBEC low"
di$apobec_high <- ifelse(di$apobec > 2, 1, 0)

tcga <- di %>%
  select(sample_id, apobec_high) %>%
  mutate(sample = str_replace_all(substr(sample_id,1,12), "-", "\\.")) %>%
  inner_join(tcga, by = "sample")

tcga$apobec_high <- factor(ifelse(tcga$apobec_high == 1, "High", "Low"), levels = c("Low", "High"))


#ROR-P
tcga$High_RORP <- ifelse(tcga$ROR.P.Group..Subtype...Proliferation. == "high", 1, 0)
tcga$Med_RORP <- ifelse(tcga$ROR.P.Group..Subtype...Proliferation. == "med", 1, 0)
tcga$Low_RORP <- ifelse(tcga$ROR.P.Group..Subtype...Proliferation. == "low", 1, 0)

#APOBEC high vs low (low is referent) and look at subtype, 
#node positivity, stage, grade, p53 mutant status (both RNA and IHC if you have for TCGA), ER, PR, HER2, race, age
tcga$p53_IHC <- tcga$tp53



# Run Models to Get RFDs --------------------------------------------------


#function that will add the point estimate (95% CI), panel (184 or 51), outcome, and level to a dataframe
addRes <- function(var, df){
  
  fmla <- paste0(var, "~ apobec_high")
  
  mod <-  glm(formula = fmla,
              family=binomial(link="identity"),
              data=df)
  
  #confidence interval
  ci <- confint(mod)
  
  #coefficients
  est <- summary(mod)$coefficients
  
  #bind everything together
  res <- cbind.data.frame(est, ci, rep(var, nrow(est)))
  
  
  #give consistent names
  names(res) <- c("estimate", "std.error", "Z", "p.value", "conf.low", "conf.high", "outcome")
  
  #return it
  return(res)
}


tcga <- tcga %>%
  mutate(p53 = ifelse(p53_RNA == "WT-like", 0, 1),
         p53_IHC = ifelse(p53_IHC == "NO", 0, 1),
         race = ifelse(race3 == "Non-AA", 0, 1),
         er_pos = ifelse(er_pos == "Positive", 0, 1),
         menopausalstatus = ifelse(menopause == "Post", 0, 1),
         node_status = ifelse(node == "Negative", 0, 1),
         her2_positive = ifelse(her2_status != "Positive", 0, 1),
         age_cat = ifelse(Age_Cat == "<=50", 0, 1))

tcga <- tcga %>%
  mutate(lumA = ifelse(pam50 == "LumA", 1, 0),
         lumB = ifelse(pam50 == "LumB", 1, 0),
         norm = ifelse(pam50 == "Normal", 1, 0),
         her2 = ifelse(pam50 == "Her2", 1, 0),
         basal = ifelse(pam50 == "Basal", 1, 0))

tcga <- tcga %>%
  mutate(stageI = ifelse(stage == "Stage I", 1, 0),
         stageII = ifelse(stage == "Stage II", 1, 0),
         stageIII = ifelse(stage == "Stage III", 1, 0),
         stageIV = ifelse(stage == "Stage IV", 1, 0))

tcga_df_un <- map_dfr(.x = list("p53", "age_cat", "race", "er_pos", "menopausalstatus", "node_status", "her2_positive",
                                "lumA", "lumB", "norm", "her2", "basal", "stageI", "stageII", "stageIII", "stageIV", "High_RORP",
                                "Med_RORP", "Low_RORP", "p53_IHC"),
                      .f = ~addRes(var = .x, df = tcga))


tcga_df_un <- tcga_df_un %>%
  mutate(predictor = rownames(.)) %>%
  filter(!str_detect(predictor, "ntercept")) %>%
  mutate(predictor = str_remove_all(predictor,  "\\.(.*)")) %>%
  mutate(predictor = str_remove_all(predictor, "apobec_high"))

#add observations that represent the reference group of apobec low group
#add points for the reference group,

refgroups <- cbind.data.frame(rep(0,20),rep(0,20),rep(0,20),
                 rep(0,20),rep(0,20),rep(0,20),
                 tcga_df_un$outcome, rep("Low", 20))
names(refgroups) <- names(tcga_df_un)

tcga_df_un <- tcga_df_un %>%
  bind_rows(refgroups)

#add a group label for what the outcome is part of
tcga_df_un <- tcga_df_un %>%
  mutate(group = case_when(outcome == "p53" ~ "P53",
                           outcome == "race" ~ "Race",
                           outcome == "er_pos" ~ "ER Positive",
                           outcome == "menopausalstatus" ~ "Menopausal",
                           outcome == "node_status" ~ "Node Status",
                           outcome == "her2_positive" ~ "HER2",
                           outcome == "age_cat" ~ "Age",
                           str_detect(outcome, "stage") ~ "Stage",
                           str_detect(outcome, "RORP") ~ "RORP Group",
                           T ~ "PAM50")) 
tcga_df_un <- tcga_df_un %>%
  mutate(estimate = as.numeric(estimate),
         lowerCI = as.numeric(conf.low),
         upperCI = as.numeric(conf.high)) %>%
  #want some sort of ordering indicator sof reference group always ends up on the bottom
  mutate(ordering = case_when(outcome == "age_cat" ~ 1,
                              outcome == "p53" ~ 2,
                              outcome == "er_pos" ~ 3,
                              outcome == "lumA" ~ 4,
                              outcome == "lumB" ~ 5,
                              outcome == "her2" ~ 6,
                              outcome == "basal" ~ 7,
                              outcome == "norm" ~ 8,
                              outcome == "her2_positive" ~ 12,
                              outcome == "node_status" ~ 13,
                              outcome == "menopausalstatus" ~ 14, 
                              outcome == "race" ~ 15,
                              outcome == "stageIV" ~ 19,
                              outcome == "stageIII" ~ 18,
                              outcome == "stageII" ~ 17,
                              outcome == "stageI" ~ 16,
                              outcome == "Low_RORP" ~ 20,
                              outcome == "Med_RORP" ~ 21,
                              outcome == "High_RORP" ~ 22))%>%
  mutate(outcome = case_when(outcome == "race" ~ "African American",
                             outcome == "menopausalstatus" ~ "Premenopausal",
                             outcome == "her2_positive" ~ "HER2 Positive",
                             outcome == "norm" ~ "Normal",
                             outcome == "er_pos" ~ "ER Negative",
                             outcome == "age_cat" ~ "50 + Years",
                             outcome == "node_status" ~ "Node Positive",
                             outcome == "basal" ~ "Basal",
                             outcome == "her2" ~ "Her2",
                             outcome == "lumB" ~ "LumB",
                             outcome == "lumA" ~ "LumA",
                             str_detect(outcome, "stage") ~ str_replace(outcome, "stage", "Stage "),
                             outcome == "High_RORP" ~ "High",
                             outcome == "Med_RORP" ~ "Intermediate",
                             outcome == "Low_RORP" ~ "Low",
                             outcome == "p53" ~ "Mutant-Like",
                             T ~ outcome))


# Make Graph -------------------------------------------------------------


pdf("data/dimarco/tcga_rfd_apobec_plot.pdf")

tcga_df_un %>%
 # filter(outcome %in% c("African American", "p53", "LumB", "LumA", "Basal", "Her2"))%>%
  ggplot(aes(y = estimate, 
             ymin = lowerCI,
             ymax = upperCI, 
             x = reorder(outcome, ordering), 
             color = predictor, 
             group = predictor)) +
  geom_pointrange(position =position_dodge(width = .4), size = 1, shape =20) + 
  coord_flip() + 
  labs(y = "RFD (95% CI)", x = "", color = "DNA Group") + 
  theme_bw() +
  theme(plot.caption = element_text(hjust = .5),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank()
        )+
  geom_hline(yintercept = 0, linetype = "dotdash", color = "grey") +
  facet_grid(group~., scales = "free_y", space = "free") +
  scale_color_manual(values = c("purple","black"))

dev.off()



# Focus on Immune Clusters ------------------------------------------------


#the immune clusters from Ashley
imm <- readxl::read_excel("data/dimarco/Supplementary File 2_v1.xls")

#per Ashley's email "we're interested in seeing why 1/2 the APOBEC-high HER2 tumors are in cluster 1 and half are in cluster 2.
#probably makes sense then to just look at the APOBEC-high HER2 tumors and do RFDs btwn cluster 1 and cluster 2 then
tcga2 <- imm %>%
  mutate(sample = str_replace_all(substr(Sample_ID, 1, 12), "-", "\\.")) %>%
  right_join(tcga, by = "sample")


tcga2$her2_immune_cluster <- ifelse(tcga2$Cluster_Number == "APOBEC-high HER2-1", 1,
                                    ifelse(tcga2$Cluster_Number == "APOBEC-high HER2-2", 0, 
                                           ifelse(tcga2$Cluster_Number == "APOBEC-high Basal-1",3,
                                                  ifelse(tcga2$Cluster_Number == "APOBEC-high Basal-2", 4, NA))))

tcga2$her2_basal <- ifelse(tcga2$her2_immune_cluster %in% c(1,0), "HER2", 
                           ifelse(tcga2$her2_immune_cluster %in% c(3,4), "Basal", NA))

tcga2$cluster_only <- ifelse(tcga2$her2_immune_cluster %in% c(1,3), "Cluster 1",
                             ifelse(tcga2$her2_immune_cluster %in% c(0,4), "Cluster 2", NA))


#function that will add the point estimate (95% CI), panel (184 or 51), outcome, and level to a dataframe
addRes2 <- function(var, df){
  
  fmla <- paste0("her2_immune_cluster ~", var)
  
  mod <-  glm(formula = fmla,
              family=binomial(link="identity"),
              data=df)
  
  #confidence interval
  ci <- confint(mod)
  
  #coefficients
  est <- summary(mod)$coefficients
  
  #bind everything together
  res <- cbind.data.frame(est, ci, rep(var, nrow(est)))
  
  
  #give consistent names
  names(res) <- c("estimate", "std.error", "Z", "p.value", "conf.low", "conf.high", "outcome")
  
  #return it
  return(res)
}


tcga_her2 <- map_dfr(.x = list("p53", "age_cat", "race", "er_pos", "menopausalstatus", "node_status", "her2_positive",
                                "lumA", "lumB", "norm", "her2", "basal", "stageI", "stageII", "stageIII", "stageIV", "High_RORP",
                                "Med_RORP", "Low_RORP", "p53_IHC"),
                      .f = ~addRes2(var = .x, df = tcga2))

names(tcga2)

table(her2 = tcga2$her2_immune_cluster, other = tcga2$p53) #%>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$race2) #%>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$er_pos)# %>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$menopausalstatus)# %>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$node_status)# %>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$her2_positive)# %>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$stage) #%>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$p53_IHC) #%>% fisher.test()
table(her2 = tcga2$her2_immune_cluster, other = tcga2$ROR.P.Group..Subtype...Proliferation.)
table(her2 = tcga2$her2_immune_cluster, other = tcga2$age_cat)
table(her2 = tcga2$her2_immune_cluster, other = tcga2$node_status)
table(her2 = tcga2$her2_immune_cluster, other = tcga2$er_pos)
table(her2 = tcga2$her2_immune_cluster, other = tcga2$age_cat)

#maybe make a bar chart showing the counts for each in immune classes 1 and 2?

mydf <- rbind.data.frame(c("P53 (IHC)", "Mutant-Like", "Class 2", 10, "HER2"),
                         c("P53 (IHC)", "Mutant-Like", "Class 1", 14, "HER2"),
                         c("Age"      , "<= 50 years", "Class 2", 16, "HER2"),
                         c("Age"      , "<= 50 years", "Class 1", 15, "HER2"),
                         c("Stage"    , "Stage I"    , "Class 2", 2, "HER2"),
                         c("Stage"    , "Stage I"    , "Class 1", 3, "HER2"),
                         c("Stage"    , "Stage II"   , "Class 2", 15, "HER2"),
                         c("Stage"    , "Stage II"   , "Class 1", 12, "HER2"),
                         c("Stage"    , "Stage III"  , "Class 2", 3, "HER2"),
                         c("Stage"    , "Stage III"  , "Class 1", 6, "HER2"),
                         c("ROR-P"    ,"Intermediate", "Class 2", 15, "HER2"),
                         c("ROR-P"    ,"Intermediate", "Class 1", 9, "HER2"),
                         c("ROR-P"    , "High"      , "Class 2", 7, "HER2"),
                         c("ROR-P"    , "High"      , "Class 1", 12, "HER2"),
                         c("Node"     , "Postitive" , "Class 2", 6, "HER2"),
                         c("Node"     , "Postitive" , "Class 1", 7, "HER2"),
                         c("Menopause", "Post-Menopausal", "Class 2", 6, "HER2"),
                         c("Menopause", "Post-Menopausal", "Class 1", 6, "HER2"),
                         c("ER"       , "Positive"  , "Class 2", 11, "HER2"),
                         c("ER"       , "Positive"  , "Class 1", 13, "HER2"),
                         c("Race"     , "White"    , "Class 2", 12, "HER2"),
                         c("Race"     , "White"    , "Class 1", 7, "HER2"),
                         c("P53 (IHC)", "Mutant-Like", "Class 2", 1, "Basal"),
                         c("P53 (IHC)", "Mutant-Like", "Class 1", 31, "Basal"),
                         c("Age"      , "<= 50 years", "Class 2", 0, "Basal"),
                         c("Age"      , "<= 50 years", "Class 1", 19, "Basal"),
                         c("Stage"    , "Stage I"    , "Class 2", 0, "Basal"),
                         c("Stage"    , "Stage I"    , "Class 1", 5, "Basal"),
                         c("Stage"    , "Stage II"   , "Class 2", 1, "Basal"),
                         c("Stage"    , "Stage II"   , "Class 1", 24, "Basal"),
                         c("Stage"    , "Stage III"  , "Class 2", 0, "Basal"),
                         c("Stage"    , "Stage III"  , "Class 1", 2, "Basal"),
                         c("ROR-P"    ,"Intermediate", "Class 2", 1, "Basal"),
                         c("ROR-P"    ,"Intermediate", "Class 1", 8, "Basal"),
                         c("ROR-P"    , "High"      , "Class 2", 0, "Basal"),
                         c("ROR-P"    , "High"      , "Class 1", 26, "Basal"),
                         c("Node"     , "Postitive" , "Class 2", 0, "Basal"),
                         c("Node"     , "Postitive" , "Class 1", 11, "Basal"),
                         c("Menopause", "Post-Menopausal", "Class 2", 1, "Basal"),
                         c("Menopause", "Post-Menopausal", "Class 1", 8, "Basal"),
                         c("ER"       , "Positive"  , "Class 2", 1, "Basal"),
                         c("ER"       , "Positive"  , "Class 1", 13, "Basal"),
                         c("Race"     , "White"    , "Class 2", 0, "Basal"),
                         c("Race"     , "White"    , "Class 1", 16, "Basal"))

names(mydf) <- c("predictor", "level", "immune_class", "count", "subtype")


pdf("data/dimarco/counts_for_her2_apobec_high.pdf")
mydf %>%
  mutate(count = as.numeric(count))%>%
  ggplot(aes(y = count, x = level, fill = immune_class)) +
  geom_bar(position = "dodge", stat = "identity", width = .4) +
  coord_flip() +
  facet_grid(predictor~subtype, space = "free", drop = T, scales = "free_y") +
  theme_bw() + theme(strip.text = element_text(size = 7))+
  labs(fill = "Immune Cluster", x = "", y = "N")
dev.off()

pdf("data/dimarco/counts_for_her2_apobec_high.pdf")
mydf %>%
  mutate(count = as.numeric(count),
         immune_class = paste0(subtype, ": ", immune_class))%>%
  ggplot(aes(y = count, x = level, fill = immune_class)) +
  geom_bar(position = "dodge", stat = "identity", width = .4) +
  coord_flip() +
  facet_grid(predictor~., space = "free", drop = T, scales = "free_y") +
  theme_bw() + theme(strip.text = element_text(size = 7))+
  labs(fill = "Immune Cluster", x = "", y = "N")
dev.off()


#in addition to just focusing on the HER2 tumors, it could be helpful to at least look at ALL




