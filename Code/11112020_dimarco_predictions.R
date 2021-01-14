#################################################
# Sarah C. Van Alsten                           #
# Date: November 11, 2020                       #
# Institution: UNC Chapel Hill                  #
# Purpose: Predict Signature Calls for DiMarco  #
#          Immune Paper                         #
# Packages: tidyverse                           #
#################################################

library(tidyverse)

#She said she could send us her TCGA IDs with the APOBEC enrichment scores (just in case they are different from ours) 
#and I figured we could see if our predictor pulled out those calls


#read in the data
di <- readxl::read_excel("data/dimarco/DiMarco-APOBEC-Immune-Signature-Data_original.xlsm")
di <- janitor::clean_names(di)
#Per Andrea's email: "The enrichment score of 2 is used as a cutoff for "APOBEC high" and "APOBEC low"
di$apobec_high <- ifelse(di$apobec > 2, 1, 0)




y <- read.table("data/tcga_brca_all_signature_versions.txt", header = 1, sep = "\t")

#limit the observations in di to those that are also in y where we have a prediction for APOBEC signature
table(di$sample_id %in% substr(y$Barcode.x, 1, 15))
di$sample_id[which(!di$sample_id %in% substr(y$Barcode.x, 1, 15))]
y$Barcode.x[which(!substr(y$Barcode.x, 1, 15) %in% di$sample_id)]

which(!substr(y$Barcode.x, 1, 15) %in% di$sample_id)

#use the FULL set of available genes, not just the ones in the 184 panel
x2 <- read.table("data/TCGA_BRCA_1201_Tumors_Log.txt", header = T, sep = "\t")
rownames(x2) <- x2[,1]
x2 <- x2[,-1]


#need to be median centered and then imputed
source("C:/Users/Owner/OneDrive/Documents/UNC/Research/dna_repair/code/clanc.R")
source("C:/Users/Owner/OneDrive/Documents/UNC/Research/IMMUNE/immune_signatures/code/arrayTools_6_22.R")

#x2 <- x2 %>% medianCtr()

#write it out, then impute with readarray
#write.table(x2, "data/alltcgalogmed.txt", sep = "\t", col.names = NA)

#x3 <- readarray("data/alltcgalogmed.txt")
write.table(x3$xd, "data/alltcgalogmedimpute.txt")

#find the most variable genes
geneList <- rownames(x3$xd)

#if we want UNADJUSTED scores output instead, use this function instead
fullUnadjFun <- function(baseString, geneList1 = geneList, data1 = y){

  A2 <- read.table("data/alltcgalogmedimpute.txt")
  
  #extract apobec calls
  apobec.sig <- data1$apobec_high
  
  #can't have missing values for the outcome variables in SAM.
  #make subsets that don't include missings
  apobec.A <- A2[,!is.na(apobec.sig)]

  ##########################################################################

  data.apobec.a <- list(x=as.matrix(apobec.A), y=apobec.sig[!is.na(apobec.sig)], genenames=as.character(rownames(A2)), logged2=TRUE)
  samr.apobec.a <- samr(data.apobec.a, resp.type="Two class unpaired", nperms=500, random.seed = 240957)
  
  
  # Get Delta Tables + Plots ------------------------------------------------
  
  delta.tables <- purrr::map(list(samr.apobec.a),
                             .f = samr.compute.delta.table)
  
  #select the appropriate delta values
  del1 <- delta.tables[[1]] %>% 
    as.data.frame() %>% 
    arrange(`# med false pos`, delta) %>%
    filter(delta !=0) %>%
    slice(1) %>%
    dplyr::select(delta)%>%
    as.numeric()
  
  #Plots of differentially expressed genes
  plots<- map2(.x = list(samr.apobec.a),
               .y = list(del1),
               .f = samr.plot)
  
  #get significant genes
  siggenes.table <- pmap(.l = list(samr.obj = list(samr.apobec.a),
                                   del = list(del1),
                                   data = list(data.apobec.a),
                                   delta.table = delta.tables),
                         .f = samr.compute.siggenes.table, all.genes = T)
  
  #get p-values for all of these
  newSig <- purrr::map(siggenes.table, getP)
  
  outNames <- c(paste0("data/dimarco/12152020_apobec_", baseString, ".txt"))
  
  
  outNames <- as.list(outNames)
  
  #write out the sig genes
  map2(.x = newSig, .y =outNames,
       .f = writeSigGene)
  
  #read in and bind as volcano plot format
  volc <- map_dfr(.x = outNames,
                  .f = readSigGene)
  

  
  names(volc) <- c("x", "row", "gene_name", "gene_id", "t_stat", "numerator", "denominator",
                   "fold_change", "q_value", "p_value", "expression", "signature")
  
  
  volc$name <- sapply(str_split(volc$gene_name, "\\|"), "[", 1)
  volc$name <- toupper(volc$name)
  return(volc)
  
  
}


#thus, need to compute a p-value based on the t-statistic we get and compare to the q-value
getP <- function(sigg){
  sigg$genes.up<- cbind(sigg$genes.up, pt(q = abs(as.numeric(sigg$genes.up[,4])), df = 992, lower.tail = F))
  sigg$genes.lo<- cbind(sigg$genes.lo, pt(q = abs(as.numeric(sigg$genes.lo[,4])), df = 992, lower.tail = F))
  return(sigg)
}

#read in the data that was just written out and bind it into one large dataframe
readSigGene <- function(filename){
  
  lo <- read.delim(paste0(filename, "_low.txt")) %>% 
    mutate(express = "Down") %>%
    mutate(signature = substr(filename, 6, nchar(filename)))
  
  hi <- read.delim(paste0(filename, "_high.txt")) %>%
    mutate(express = "Up") %>%
    mutate(signature = substr(filename, 6, nchar(filename)))
  
  return(rbind.data.frame(lo, hi))
  
}


#write out each of the significant gene tables: one for up regulate, one for down regulated
writeSigGene <- function(sigobj, filename){
  write.table(sigobj$genes.up, paste0(filename, "_high.txt"), sep='\t', col.names=NA)
  write.table(sigobj$genes.lo, paste0(filename, "_low.txt"), sep='\t', col.names=NA)
  
}



m1 <- fullUnadjFun(baseString = "test")

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .1)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]
write.table(geneSub, "data/dimarco/top10percentVariable.txt", sep = "\t", col.names = NA)

y$barcode_short2 <- str_replace_all(substr(y$barcode_short, 1, 16), "-", "\\.")

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]


# di <- di %>% filter(sample_id %in% substr(y$Barcode.x, 1, 15))
# y <- y %>% filter(substr(Barcode.x, 1, 15) %in% di$sample_id)
# y$Barcode.x <- substr(y$Barcode.x, 1, 15)
# 
# y <- di %>%
#   dplyr::select(sample_id, apobec_high) %>%
#   full_join(y, by = c("sample_id" = "Barcode.x"))
# 
# table(y$Apo_targeted, y$apobec_high)
# table(y$APOBEC_wes, y$apobec_high)
# table(y$apo_expr, y$apobec_high)



# Helper Functions -----------------------------------------------

percCorr <- function(tab){
  
  total <- sum(tab[1:2,1:2])
  perCorr <-(sum(tab[1,1], tab[2,2])/total)*100
  return(perCorr)
  
  
}

sens <- function(tab){
  
  rowTot <- sum(tab[2,])
  perCorr <- (tab[2,2]/rowTot)*100
  return(perCorr)
}

spec <- function(tab){
  
  rowTot <- sum(tab[1,])
  perCorr <- (tab[1,1]/rowTot)*100
  return(perCorr)
}

ppv <- function(tab){
  
  rowTot <- sum(tab[,2])
  perCorr <- (tab[2,2]/rowTot)*100
  return(perCorr)
}

npv <- function(tab){
  
  rowTot <- sum(tab[,1])
  perCorr <- (tab[1,1]/rowTot)*100
  return(perCorr)
}

allStats <- function(tab){
  if (ncol(tab)==2 & nrow(tab) == 2){
    return(list(perc_correct = percCorr(tab),
                sensitivity = sens(tab),
                specificity = spec(tab),
                ppv = ppv(tab),
                npv = npv(tab)))
  } else {
    if (ncol(tab) ==1 & nrow(tab)==2){
      if (colnames(tab) == "1"){
        perc_correct = tab[1,1]/sum(tab[,1])
        sensitivity = NA
        specificity = tab[1,1]/sum(tab[,1])
        ppv = NA
        npv =  tab[2,1]/sum(tab[,1])
        return(list(perc_correct = perc_correct,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    ppv = ppv,
                    npv = npv))
      
      } else {
        perc_correct = tab[1,1]/sum(tab[,1])
        specificity = NA
        sensitivity = tab[1,1]/sum(tab[,1])
        npv = NA
        ppv =  tab[2,1]/sum(tab[,1])
        return(list(perc_correct = perc_correct,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    ppv = ppv,
                    npv = npv))
      }
    } else if (ncol(tab) ==1 & nrow(tab)==1){
      if (colnames(tab) == "1"){
        perc_correct = 100
        sensitivity = NA
        specificity = 100
        ppv = NA
        npv =  100
        return(list(perc_correct = perc_correct,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    ppv = ppv,
                    npv = npv))
        
      } else {
        perc_correct = 100
        specificity = NA
        sensitivity = 100
        npv = NA
        ppv =  100
        return(list(perc_correct = perc_correct,
                    sensitivity = sensitivity,
                    specificity = specificity,
                    ppv = ppv,
                    npv = npv))
      }
    }
    
  }

}



#i is the number of genes per groups
#x1 is the complete set of genes to evaluate
#signature is which of the mutation signatures to evaluate
runTune <- function(i, x1 = geneSub, signature = "apobec_high"){
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  mysig <- y[, names(y) == signature]
  
  Pre <-clancPredict(x1[,!is.na(mysig)],mysig[!is.na(mysig)],
                     x1[,!is.na(mysig)], i) 

  
  t1 <- table(mysig[!is.na(mysig)],Pre$predictions[,2])
  
  t1Stats <- allStats(t1)
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i))
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes")
  
  return(df)
  
}

#dat= data frame returned by runTune or runTuneTrain
#signature = the mutational signature to label on the graph
rocPlot <- function(dat, signature){
  #plot an ROC curve for using any genes, including those with lots of missing data
  p1 <- dat %>%
    filter(mod_type == "full")%>%
    mutate(spec1 = 1-as.numeric(specificity)/100) %>%
    ggplot(aes(x = as.numeric(specificity)/100, color = as.numeric(ngenes), y = as.numeric(sensitivity)/100))+
    geom_point() + 
    theme_bw() + 
    labs(color = "Number Genes", x = "Specificity", y = "Sensitivity") +
    scale_color_viridis_c() +
    ggtitle(paste0("All Genes: ", signature))
  
  return(p1)
}

#extra is whether to return a particular model (for a certain model)
#i is the number of genes per groups
#x1 is the complete set of genes to evaluate
#x3 is the set of genes with >50% non missing
#signature is which of the mutation signatures to evaluate
runTuneTrain <- function(i, x1 = geneSub, extra = NULL, signature = "apobec_high"){
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  if (signature == "Aging"){
    y$samp <- ifelse(rownames(y) %in% which(y$Aging == 2), 1, 0)
  } else if (signature == "UV"){
    y$samp <- ifelse(rownames(y) %in% which(y$UV == 2), 1, 0)
  } else if (signature == "HR"){
    y$samp <-ifelse(rownames(y) %in% which(y$HR == 2), 1, 0)
  } else if (signature == "APOBEC"){
    y$samp <- ifelse(rownames(y) %in% which(y$APOBEC == 2), 1, 0)
  } else if (signature == "APOBEC3A3B"){
    y$samp <- ifelse(rownames(y) %in% which(y$APOBEC3A3B == 2), 1, 0)
  } else if (signature == "apobec_high"){
    y$samp <- ifelse(rownames(y) %in% which(y$apobec_high == 2), 1, 0)
  }
  
  train <- y %>%
    mutate(rown = row_number())%>%
    group_by(samp)%>%
    slice_sample(prop = .7)
  test <- y %>%
    mutate(rown = row_number())%>%
    filter(!rown %in% train$rown)
  
  #now subset x1 and x3 into training and test sets
  x1train <- x1[, train$rown]
  
  x1test <- x1[, test$rown]

  sigCol <- train[,names(train)== signature]
  sigCol2 <- test[,names(test)== signature]
  
  if (!is.null(extra)){
    Pre <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                       x1train[,!is.na(sigCol)], extra)
    return(Pre)
  }
  
  
  #in sample prediction
  Pre <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                     x1train[,!is.na(sigCol)], i) 
  
  t1 <- table(sigCol[!is.na(sigCol)],Pre$predictions[,2])
  
  t1Stats <- allStats(t1)
  
  #out of sample prediction
  Pre3 <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                      x1test[,!is.na(sigCol2)], i) 

  
  t1 <- table(sigCol[!is.na(sigCol)],Pre$predictions[,2])
  t3 <- table(sigCol2[!is.na(sigCol2)],Pre3$predictions[,2])
  
  t1Stats <- allStats(t1)
  t3Stats <- allStats(t3)
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i, "in_sample"),
                         c(unlist(allStats(t3)), "full", i, "out_sample"))
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes", "in_out")
  
  return(df)
  
}

rocPlotTrain <- function(dat, signature){
  #plot an ROC curve for using any genes, including those with lots of missing data
  p1 <- dat %>%
    filter(mod_type == "full")%>%
    mutate(spec1 = 1-as.numeric(specificity)/100) %>%
    ggplot(aes(x = as.numeric(specificity)/100, color = as.numeric(ngenes), y = as.numeric(sensitivity)/100))+
    geom_point() + 
    theme_bw() + 
    facet_grid(.~in_out)+
    labs(color = "Number Genes", x = "Specificty", y = "Sensitivity") +
    scale_color_viridis_c() +
    ggtitle(paste0("All Genes: ", signature))
  

  return(p1)
}

makeConcPlot <- function(dat, sens){
  
  if (sens == "sensitivity"){
    return(
      dat %>%
        ggplot(aes(x = ngenes, y = as.numeric(sensitivity), color = in_out)) +
        geom_point() +
        theme_bw() +
        #ggtitle("All Genes") +
        labs(color = "Prediction Type", x = "Number Genes Per Group", y = "Sensitivity")
    )
  } else {
    return(
      dat %>%
        ggplot(aes(x = ngenes, y = as.numeric(specificity), color = in_out)) +
        geom_point() +
        theme_bw() +
        #ggtitle("All Genes") +
        labs(color = "Prediction Type", x = "Number Genes Per Group", y = "Specificity")
    )
  }
  
}

concordancePlot <- function(dat, signature2){
  
  #all
  p1 <- dat %>%
    filter(mod_type == "full")%>%
    mutate(ngenes = as.numeric(ngenes))%>%
    makeConcPlot("sensitivity") +
    ggtitle(paste0("All Genes: Sensitivity, ", signature2))
  
  
  #all
  p3 <-  dat %>%
    filter(mod_type == "full")%>%
    mutate(ngenes = as.numeric(ngenes))%>%
    makeConcPlot("specificity") +
    ggtitle(paste0("All Genes: Specificity, ", signature2))
  
  
  return(list(p1,p3))
  
}


#function to do everything: run the algorithm on full sample, make plots
#then run algorithm on test/training set and make plots; return both data frames
runClanc <- function(signature1, genes = NULL){
  
  if(is.null(genes)){
  #for all genes
  myTune <- purrr::map_dfr(.x = 1:100, .f = ~runTune(i = .x, signature = signature1))
  plot1 <- rocPlot(myTune,  signature1)
  

  #for the full set of genes, retune
  myTuneTT <- purrr::map_dfr(1:100, .f = ~runTuneTrain(i = .x, signature = signature1))
  plot3 <- rocPlotTrain(myTuneTT,  signature1)
  cplot3 <- concordancePlot(myTuneTT, signature1)
  
  
  return(list(all.genes.df = myTune, roc1 = plot1,
              all.genes.df.TT = myTuneTT, roc3 = plot3, concordance3 = cplot3))
  } else {
    #for all genes
    myTune <- purrr::map_dfr(.x = 1:100, .f = ~runTune(i = .x, signature = signature1,x1 = genes))
    plot1 <- rocPlot(myTune,  signature1)
    
    
    #for the full set of genes, retune
    myTuneTT <- purrr::map_dfr(1:100, .f = ~runTuneTrain(i = .x, signature = signature1, x1 = genes))
    plot3 <- rocPlotTrain(myTuneTT,  signature1)
    cplot3 <- concordancePlot(myTuneTT, signature1)
    
    
    return(list(all.genes.df = myTune, roc1 = plot1,
                all.genes.df.TT = myTuneTT, roc3 = plot3, concordance3 = cplot3))
    
  }
}




#y$apobec_high <- factor(ifelse(y$apobec_high == 1, 2, 1))
apo.pred <- runClanc("apobec_high")
apo.pred.df <- apo.pred[[1]] %>%
  mutate(youden = as.numeric(sensitivity)+ as.numeric(specificity)-1)


apo.pred.df %>%
  #filter(as.numeric(ngenes)<1) %>%
  arrange(desc(youden)) %>%
  slice(1:10)


apo.pred[[2]]

apo.tune.df <- apo.pred[[3]]

apo.tune.df %>%
  #filter(as.numeric(ngenes)<15) %>%
  mutate(youden = as.numeric(sensitivity)+ as.numeric(specificity)-1)%>%
  arrange(desc(youden)) %>%
  slice(1:20) # 10 genes performs well on both in sample and out of sample

# 10 Fold Cross Validation to Get Bounds for the Sens/Spec ----------------

#function to do 10 fold CV of the tune/train in order to generate sort of confidence intervals
#for each of the points along the number of genes per group
runClancCV <- function(signature1){
  
  for (i in 1:10){
    if (i==1){
      full <-  purrr::map_dfr(.x = 1:50, .f = ~runTuneTrain(i = .x, x1 = geneSub, signature = signature1))
    } else{
      f <- purrr::map_dfr(.x = 1:50, .f = ~runTuneTrain(i = .x, x1 = geneSub, signature = signature1))
      full <- bind_rows(full, f)
    }
    
  }
  return(full)
}

runClancCV2 <- function(signature1){
  
  for (i in 1:10){
    if (i==1){
      full <-  purrr::map_dfr(.x = 1:50, .f = ~runTuneTrain(i = .x, x1 = x3$xd, signature = signature1))
    } else{
      f <- purrr::map_dfr(.x = 1:50, .f = ~runTuneTrain(i = .x, x1 = x3$xd, signature = signature1))
      full <- bind_rows(full, f)
    }
    
  }
  return(full)
}

plotCV <- function(d, signature){
  d2a <-  d %>%
    filter(mod_type=="full")%>%
    mutate(sensitivity = as.numeric(sensitivity),
           specificity = as.numeric(specificity),
           ngenes = as.numeric(ngenes))%>%
    mutate(in_out = case_when(in_out == "in_sample"~ "Training Set",
                              T ~ "Test Set")) %>%
    group_by(ngenes) %>%
    summarise(avg_sens = mean(sensitivity),
              avg_spec = mean(specificity),
              min_sens = min(sensitivity),
              max_sens = max(sensitivity),
              min_spec = min(specificity),
              max_spec = max(specificity),
              sd_sens = sd(sensitivity),
              sd_spec = sd(specificity),
              .groups = "keep") 
  
  d2b <-  d %>%
    filter(mod_type=="full" & in_out != "in_sample")%>%
    mutate(sensitivity = as.numeric(sensitivity),
           specificity = as.numeric(specificity),
           ngenes = as.numeric(ngenes))%>%
    mutate(in_out = case_when(in_out == "in_sample"~ "Training Set",
                              T ~ "Test Set")) %>%
    group_by(ngenes) %>%
    summarise(avg_sens = mean(sensitivity),
              avg_spec = mean(specificity),
              min_sens = min(sensitivity),
              max_sens = max(sensitivity),
              min_spec = min(specificity),
              max_spec = max(specificity),
              sd_sens = sd(sensitivity),
              sd_spec = sd(specificity),
              .groups = "keep")
  
  d2 <- bind_rows(d2a, d2b)
  d2$in_out <- c(rep("Training Set", nrow(d2a)),
                 rep("Test Set", nrow(d2b)))
  
  
  sensPlot <- d2 %>%
    ggplot(aes(x = ngenes, y = avg_sens, color = in_out)) +
    geom_point(position = position_dodge(width = .5))+
    geom_ribbon(aes(ymin = (avg_sens - sd_sens), ymax = (avg_sens + sd_sens), fill = in_out), position = position_dodge(width = .5),
                alpha = .1) +
    theme_bw() +
    #ggtitle(paste0("10-Fold Cross Validation of Sensitivity for ", signature)) +
    labs(color = "Prediction Type", x = "Genes Per Group", y = "Sensitivity (%)") +
    geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed") +
    guides(fill = F) + ylim(c(45,80))
  
  specPlot <- d2 %>%
    ggplot(aes(x = ngenes, y = avg_spec, color = in_out)) +
    geom_point(position = position_dodge(width = .5))+
    geom_ribbon(aes(ymin = (avg_spec - sd_spec), ymax = (avg_spec + sd_spec), fill = in_out), position = position_dodge(width = .5),
                alpha = .1) +
    theme_bw() +
    #ggtitle(paste0("10-Fold Cross Validation of Specificity for ", signature)) +
    labs(color = "Prediction Type", x = "Genes Per Group", y = "Specificity (%)") +
    geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed") +
    guides(fill = F)+ ylim(c(45,80))
  
  return(list(sensPlot, specPlot))
}

apo.cv <- runClancCV("apobec_high")
apo.cv3 <- runClancCV("apobec_high")
apo.cv2 <- runClancCV2("apobec_high")
apoCVplot <- plotCV(apo.cv, "apobec")
plotCV(apo.cv2, "apobec")



summarized <- apo.cv %>%
  group_by(ngenes, in_out)  %>%
  mutate(sensitivity = as.numeric(sensitivity),
         specificity = as.numeric(specificity))%>%
  summarise(avg_sens = mean(sensitivity, na.rm =T),
            avg_spec = mean(specificity, na.rm =T),
            min_sens = min(sensitivity),
            max_sens = max(sensitivity),
            min_spec = min(specificity),
            max_spec = max(specificity),
            sd_sens = sd(sensitivity),
            sd_spec = sd(specificity),
            .groups = "keep")

summarized2 <-apo.cv %>%
  mutate(youden = as.numeric(sensitivity)/100+as.numeric(specificity)/100-1)%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(youden),
            sd_y = sd(youden),
            .groups = "keep")

#4 actually looks pretty good for sensitivity here


#okay, so picking the most variably expressed genes didn't seem to do much. What if instead we just
#give it the FULL gene list and let it do it's thing?
apo.pred2 <- runClanc("apobec_high", genes = x3$xd[,which(colnames(x3$xd) %in% y$barcode_short2)])
apo.pred2[[2]]


#see how a five gene predictor does

apo5 <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], classes = y$apobec_high[!is.na(y$apobec_high)], y = geneSub, ngenes = 5)
write.table(apo5$centroids, "data/dimarco/apobec_enrichment_gene_predictor_5.txt", sep = "\t")

rownames(apo5$centroids)
table(true = y$apobec_high, pred = apo5$predictions[,2])
table(true = y$apobec_high, pred = apo5$predictions[,2]) %>% allStats()

apo19 <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], classes = y$apobec_high[!is.na(y$apobec_high)], y = geneSub, ngenes = 19)
write.table(apo5$centroids, "data/dimarco/apobec_enrichment_gene_predictor_5.txt", sep = "\t")

rownames(apo19$centroids)
table(true = y$apobec_high, pred = apo19$predictions[,2])
table(true = y$apobec_high, pred = apo19$predictions[,2]) %>% allStats()

#make a plot of the youden's index
#also plot the Youden's index, which could be a good addition to the fourth part of the figure
youden.plot <- apo.cv %>%
  mutate(youden = as.numeric(sensitivity)/100+as.numeric(specificity)/100-1)%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(youden),
            sd_y = sd(youden),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(y = avg_y, x = as.numeric(ngenes), color = in_out)) + 
  geom_point(position = position_dodge(width = .5)) + 
  geom_smooth(se = F) + 
  geom_ribbon(aes(ymin = avg_y - sd_y, ymax = avg_y + sd_y,fill = in_out), position = position_dodge(width = .5), alpha = .1) +
  theme_bw() + labs(x = "Genes Per Group", y = "Youden's Index") + guides(color = F, fill = F)

#maybe also try the F1 score: 2*((PPV*sens)/(PPV+sens))
f1.plot <- apo.cv %>%
  mutate(f1score = 2*((as.numeric(ppv)/100)*((as.numeric(sensitivity)/100))/
                        ((as.numeric(ppv)/100)+((as.numeric(sensitivity)/100)))))%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(f1score),
            sd_y = sd(f1score),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(y = avg_y, x = as.numeric(ngenes), color = in_out)) + 
  geom_point(position = position_dodge(width = .5)) + 
  geom_smooth(se = F) + 
  geom_errorbar(aes(ymin = avg_y - sd_y, ymax = avg_y + sd_y), position = position_dodge(width = .5)) +
  theme_bw() + labs(x = "Genes Per Group", y = "F1 Score") + guides(color = F)


#15 genes looks like it does pretty well. Make a figure showing sensitivity, specificity, f1 score and then a confusion matrix


apo21 <- clancPredict(x = x2[,!is.na(y$apobec_high)], classes = y$apobec_high[!is.na(y$apobec_high)], y = x2, ngenes = 21)
write.table(apo11$centroids, "data/dimarco/apobec_enrichment_gene_predictor_21.txt", sep = "\t")

rownames(apo21$centroids)
table(true = y$apobec_high, pred = apo21$predictions[,2])
table(true = y$apobec_high, pred = apo21$predictions[,2]) %>% allStats()

apo15 <- clancPredict(x = x2[,!is.na(y$apobec_high)], classes = y$apobec_high[!is.na(y$apobec_high)], y = x2, ngenes = 15)
write.table(apo15$centroids, "data/dimarco/apobec_enrichment_gene_predictor_15.txt", sep = "\t")

rownames(apo15$centroids)
table(true = y$apobec_high, pred = apo15$predictions[,2])


#make a plottable confusion matrix
conf <- data.frame(table(Expression = apo5$predictions[,2], WES = y$apobec_high))

plotTable <- conf %>%
  mutate(goodbad = ifelse(Expression == WES, "correct", "incorrect")) %>%
  group_by(WES) %>%
  mutate(prop = Freq/sum(Freq))

#plot matrix
confMatrix <- plotTable %>%
  mutate(WES = ifelse(WES == 2, "Yes", "No"),
         Expression = ifelse(Expression == 2, "Yes", "No"))%>%
  ggplot(mapping = aes(x = WES, y = Expression,
                       fill = goodbad
  )) +
  geom_tile(alpha = .3, color = "black") +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1, size = 6) +
  scale_fill_manual(values = c(correct = "#2ca25f", incorrect = "red")) +
  labs(x = "APOBEC High")+
  theme_bw() + guides(alpha = FALSE, fill = FALSE) +
  theme(plot.background = element_rect(color = "white"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())

library(patchwork)
#finally, make the overall figure
pdf("data/dimarco/tcga_classifier_results_v2.pdf")
(apoCVplot[[1]] + apoCVplot[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + confMatrix) + plot_annotation(tag_levels = "a")
dev.off()


#maybe also try to add the mouse validation to the same figure?

new.df <- data.frame(WES = c(1,2,1,2),
                     Expression = c(1,1,2,2),
                     Freq = c(7,5,5,7))

plotTable2 <- new.df%>%
  mutate(goodbad = ifelse(Expression == WES, "correct", "incorrect")) %>%
  group_by(WES) %>%
  mutate(prop = Freq/sum(Freq))

#plot matrix
confMatrix2 <- plotTable2 %>%
  mutate(WES = ifelse(WES == 2, "Yes", "No"),
         Expression = ifelse(Expression == 2, "Yes", "No"))%>%
  ggplot(mapping = aes(x = WES, y = Expression,
                       fill = goodbad
  )) +
  geom_tile(alpha = .3, color = "black") +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1, size = 6) +
  scale_fill_manual(values = c(correct = "#2ca25f", incorrect = "red")) +
  labs(x = "APOBEC High")+
  theme_bw() + guides(alpha = FALSE, fill = FALSE) +
  theme(plot.background = element_rect(color = "white"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())


cf1 <- (confMatrix + ggtitle ("TCGA")+ theme(plot.title = element_text(hjust  = .5))) + 
    (confMatrix2 + ggtitle("Mouse") + theme(plot.title = element_text(hjust  = .5)))

#finally, make the overall figure
pdf("data/dimarco/tcga_classifier_results_v3.pdf",width = 8.5, height = 10)
(apoCVplot[[1]] + apoCVplot[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + cf1) + plot_annotation(tag_levels = "a")
dev.off()

tiff(file = "data/dimarco/tcga_classifier_results_v3.tiff", width = 6600, height = 7000, units = "px", res = 800)
(apoCVplot[[1]] + apoCVplot[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + cf1) + plot_annotation(tag_levels = "a")
dev.off()

#per andrea, try removing the genes that were not available in the mouse data from TCGA and seeing how the classifier
#performs

tcga.subset <- x2[!rownames(x2) %in% c("apobec3d_140564_51", "apobec3f_200316_51", "apobec3h_164668_51","apobec3g_60489_51",
                                       "apobec3c_27350_51", "magef1_64110_cta", "znf165_7718_cta", "mtmr15_22909_184" ),]

tcga.subset.pred <- clancPredict(x = x2[!is.na(y$apobec_high)],
                                 classes = y$apobec_high[!is.na(y$apobec_high)],
                                 y = tcga.subset, 21)

t(table(pred = tcga.subset.pred$predictions[,2], true = y$apobec_high)) %>% allStats()

# Updated RNASeq Data -----------------------------------------------------

fvb <- read.table("data/dimarco/FVB_Tumors_DESeq2_Normalizedcounts_Mixture_File3.txt", sep = "\t", header = T)
nsg <- read.table("data/dimarco/NSG_Tumors_DESeq2_Normalizedcounts_Mixture.txt", sep = "\t", header = T)

#these appear to be already normalized per the name, but most likely are not median centered. Double check
fvb.genes <- fvb[,-1]
nsg.gense <- nsg[,-1]

fvb.genes <- fvb.genes %>%
  medianCtr()

nsg.gense <- nsg.gense %>%
  medianCtr()

#fvb$GeneSymbol[(str_detect(fvb$GeneSymbol, "C17"))]

rownames(fvb.genes) <- fvb$GeneSymbol
rownames(nsg.gense) <- nsg$GeneSymbol

#need to rename some of the rows so they have the same names as in x2
fvb.sub <- fvb.genes[rownames(fvb.genes) %in% toupper(sapply(str_split(rownames(apo5$centroids), "\\|"), "[[",1)),]
nsg.sub <- nsg.gense[rownames(nsg.gense) %in% toupper(sapply(str_split(rownames(apo5$centroids), "\\|"), "[[",1)),]

#nags, znf425, loc147727


# #some of the genes are missing. Which ones?
# toupper(sapply(str_split(rownames(apo21$centroids), "_"), "[[",1))[!
#           toupper(sapply(str_split(rownames(apo21$centroids), "_"), "[[",1)) %in% rownames(nsg.sub)]
# 


#RAD51D should replace RAD51L3; replace all the apobec3 with mouse apobec; faap100 = c17orf
# fvb.sub <- fvb.genes[rownames(fvb.genes) %in% c("APOBEC3","FAAP100", "RAD51D",
#                                                 toupper(sapply(str_split(rownames(apo21$centroids), "_"), "[[",1)) ),]
# nsg.sub <- nsg.gense[rownames(nsg.gense) %in% c("APOBEC3", "FAAP100", "RAD51D",
#                                                 toupper(sapply(str_split(rownames(apo21$centroids), "_"), "[[",1))),]
# 
# #duplicate the multiple version of apobec3 using the single gene that mice have (need 6 total, so five more)
# #still missing magef1 and ZNF165 but that can't be helped, mtmr15
# nsg.sub <- rbind(nsg.sub, nsg.sub[1,],nsg.sub[1,],nsg.sub[1,],nsg.sub[1,],nsg.sub[1,])
# fvb.sub <- rbind(fvb.sub, fvb.sub[1,],fvb.sub[1,],fvb.sub[1,],fvb.sub[1,],fvb.sub[1,])

#now rename to have the same names as in the human version
genelist <- sort(rownames(apo5$centroids))

#rename rows to correspond exactly to the genelsit
rownames(nsg.sub) <- rownames(fvb.sub) <-c(genelist[1], genelist[2], genelist[4:8])

"LOC100130148|100130148"
"SIRPG|55423"
"VNN2|8875"

# rownames(nsg.sub) <- rownames(fvb.sub) <- c(genelist[1], genelist[7], genelist[8], genelist[10],
#                                             genelist[11], genelist[12], genelist[13], genelist[14],
#                                             genelist[15], genelist[16], genelist[17], genelist[18],
#                                             genelist[19], genelist[9], genelist[20], genelist[21], genelist[22],
#                                             genelist[23], genelist[25],  genelist[27], genelist[28], genelist[29],
#                                             genelist[30], genelist[31], genelist[32], genelist[33], genelist[34],
#                                             genelist[35], genelist[36], genelist[37], genelist[38], genelist[39],
#                                             genelist[40], genelist[41], genelist[2], genelist[3], genelist[4], genelist[5],
#                                             genelist[6])
# 
# nsg.sub <- nsg.sub[1:34,]
# fvb.sub <- fvb.sub[1:34,]

for(i in 1:length(genelist)){
  print(noquote(genelist[i]))
}


#can now use the classifier to predict the apobec_high nature of the samples
fvb.pred <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)], y = fvb.sub,
                         ngenes = 5)
nsg.pred <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)], y = nsg.sub,
                         ngenes = 5)
fvb.pred$predictions
nsg.pred$predictions


# Can also try to train predictor on the mouse data alone -------------------

#couple of approaches to try: 1. Train on one mouse dataset, test on the other
runTuneB <- function(i, x1 = fvb.genes){
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  mysig <- c(rep(1,6), rep(2,6))
  
  Pre <-clancPredict(x1[,!is.na(mysig)],mysig[!is.na(mysig)],
                     x1[,!is.na(mysig)], i) 
  
  
  t1 <- table(mysig[!is.na(mysig)],Pre$predictions[,2])
  
  t1Stats <- allStats(t1)
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i))
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes")
  
  return(df)
  
}

#function to do everything: run the algorithm on full sample, make plots
#then run algorithm on test/training set and make plots; return both data frames
runClancB <- function(df){
  
  #for all genes
  myTune <- purrr::map_dfr(.x = 1:100, .f = ~runTuneB(i = .x,df))
  plot1 <- rocPlot(myTune,  "")
  
  
  return(list(all.genes.df = myTune, roc1 = plot1))
}



fvb.clanc <- runClancB(fvb.genes)
fvb.clanc[[1]] #basically completely perfect prediction

nsg.clanc <- runClancB(nsg.gense)
nsg.clanc[[2]] #again, basically completely perfect prediction

clancPredict(x = fvb.genes, classes = c(rep(1,6), rep(2,6)), y = nsg.gense,1)
clancPredict(y = fvb.genes, classes = c(rep(1,6), rep(2,6)), x = nsg.gense,10)

#2. 10 fold CV within a set; probably best way to do it is to
# make a 9/3 split; train on the 9 and test on the 3. Repeat 100 times to
# get idea of how well predictions are made

newTuneTrain <- function(df, i){
  
  #split into test and training set
  test <- sample(1:12, 3, F)
  train <- c(1:12)[-test]
  
  #get "ground truth" of that the apobec/control arms are
  test.status <- ifelse(test <= 6, 1, 2)
  train.status <- ifelse(train <= 6, 1, 2)
  
  
  #subset the expression data into test/train
  test.ex <- df[,test]
  train.ex <- df[,train]
  
  #in sample prediction
  Pre <-clancPredict(train.ex, train.status,
                     train.ex, i) 
  t1 <- table(train.status,Pre$predictions[,2])
  
  #out of sample prediction
  Pre2 <-clancPredict(train.ex, train.status,
                     test.ex, i) 
  t2 <- table(test.status, Pre2$predictions[,2])
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i, "in_sample"),
                         c(unlist(allStats(t2)), "full", i,  "out_sample"))
  
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes", "in_out")
  
  return(df)
  
}

newClancCV <- function(df1){
  for (i in 1:20){
    if (i == 1){
      full <-  purrr::map_dfr(.x = 1:23, .f = ~newTuneTrain(i = .x, df = df1))
      
    } else {
      f <- purrr::map_dfr(.x = 1:23, .f = ~newTuneTrain(i = .x, df = df1))
      full <- bind_rows(full, f)
    }
  }
  
  return(full)
}

fvb.tt <- newClancCV(fvb.genes) #still essentially getting perfect predictions



# Make List of Genes Consistent Between Mice and Humans --------------------

#duplicate the multiple version of apobec3 using the single gene that mice have (need 6 total, so five more)
#still missing magef1 and ZNF165 but that can't be helped, mtmr15

#read in cyrus' gene list
cyrus.gene <- readxl::read_excel("data/Copy of 2020-07-24 DNA Repair Gene List for Andrea_AW.xlsx", sheet = 4)
cyrus.gene <- c(cyrus.gene$`184 Genes`, cyrus.gene$`51 Genes`, cyrus.gene$`CTA Genes`)
cyrus.gene <- cyrus.gene[!is.na(cyrus.gene)]

#select these genes from the mouse models
fvb.cyrus <- fvb.genes[rownames(fvb.genes) %in% cyrus.gene,]
nsg.cyrus <- nsg.gense[rownames(nsg.gense) %in% cyrus.gene,]


#now, remove the various apobec3 genes, magef1, znf165, mtmr15 from consideration in x2
x2.new <- x2[!str_detect(rownames(x2), "apobec3") &
               !str_detect(rownames(x2), "magef1") & 
               !str_detect(rownames(x2), "znf165") &
               !str_detect(rownames(x2), "mtmr15"),]



#rework the same function from above but now make x2.new the default matrix
runTuneTrain <- function(i, x1 = x2.new, extra = NULL, signature = "Aging"){
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  if (signature == "Aging"){
    y$samp <- ifelse(rownames(y) %in% which(y$Aging == 2), 1, 0)
  } else if (signature == "UV"){
    y$samp <- ifelse(rownames(y) %in% which(y$UV == 2), 1, 0)
  } else if (signature == "HR"){
    y$samp <-ifelse(rownames(y) %in% which(y$HR == 2), 1, 0)
  } else if (signature == "APOBEC"){
    y$samp <- ifelse(rownames(y) %in% which(y$APOBEC == 2), 1, 0)
  } else if (signature == "APOBEC3A3B"){
    y$samp <- ifelse(rownames(y) %in% which(y$APOBEC3A3B == 2), 1, 0)
  } else if (signature == "apobec_high"){
    y$samp <- ifelse(rownames(y) %in% which(y$apobec_high == 2), 1, 0)
  }
  
  train <- y %>%
    mutate(rown = row_number())%>%
    group_by(samp)%>%
    slice_sample(prop = .7)
  test <- y %>%
    mutate(rown = row_number())%>%
    filter(!rown %in% train$rown)
  
  #now subset x1 and x3 into training and test sets
  x1train <- x1[, train$rown]
  
  x1test <- x1[, test$rown]
  
  sigCol <- train[,names(train)== signature]
  sigCol2 <- test[,names(test)== signature]
  
  if (!is.null(extra)){
    Pre <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                       x1train[,!is.na(sigCol)], extra)
    return(Pre)
  }
  
  
  #in sample prediction
  Pre <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                     x1train[,!is.na(sigCol)], i) 
  
  t1 <- table(sigCol[!is.na(sigCol)],Pre$predictions[,2])
  
  t1Stats <- allStats(t1)
  
  #out of sample prediction
  Pre3 <-clancPredict(x1train[,!is.na(sigCol)], sigCol[!is.na(sigCol)],
                      x1test[,!is.na(sigCol2)], i) 
  
  
  t1 <- table(sigCol[!is.na(sigCol)],Pre$predictions[,2])
  t3 <- table(sigCol2[!is.na(sigCol2)],Pre3$predictions[,2])
  
  t1Stats <- allStats(t1)
  t3Stats <- allStats(t3)
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i, "in_sample"),
                         c(unlist(allStats(t3)), "full", i, "out_sample"))
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes", "in_out")
  
  return(df)
  
}


runTune <- function(i, x1 = x2.new, signature = "Aging"){
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  mysig <- y[, names(y) == signature]
  
  Pre <-clancPredict(x1[,!is.na(mysig)],mysig[!is.na(mysig)],
                     x1[,!is.na(mysig)], i) 
  
  
  t1 <- table(mysig[!is.na(mysig)],Pre$predictions[,2])
  
  t1Stats <- allStats(t1)
  
  df <- rbind.data.frame(c(unlist(allStats(t1)), "full", i))
  names(df) <- c("percent_correct", "sensitivity", "specificity", "ppv", "npv", "mod_type", "ngenes")
  
  return(df)
  
}

runClancCV <- function(signature1){
  
  for (i in 1:10){
    if (i==1){
      full <-  purrr::map_dfr(.x = 1:23, .f = ~runTuneTrain(i = .x, x1 = x2.new, signature = signature1))
    } else{
      f <- purrr::map_dfr(.x = 1:23, .f = ~runTuneTrain(i = .x, x1 = x2.new, signature = signature1))
      full <- bind_rows(full, f)
    }
    
  }
  return(full)
}

#train the classifier on this subset of genes data
apo.redo <- runClanc("apobec_high")
apo.redo.cv <- runClancCV("apobec_high")

#22 genes is the only one that comes up on both in sample and out of sample. Maybe plot
apo.redo.cv %>%
  group_by(in_out)%>%
  arrange(desc(sensitivity)) %>%
  slice(1:10)

plotCV(apo.redo.cv)


#make a plot of the youden's index
#also plot the Youden's index, which could be a good addition to the fourth part of the figure
(youden.plot.new <- apo.redo.cv %>%
  mutate(youden = as.numeric(sensitivity)/100+as.numeric(specificity)/100-1)%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(youden),
            sd_y = sd(youden),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(y = avg_y, x = as.numeric(ngenes), color = in_out)) + 
  geom_point(position = position_dodge(width = .5)) + 
  geom_smooth(se = F) + 
  geom_errorbar(aes(ymin = avg_y - sd_y, ymax = avg_y + sd_y), position = position_dodge(width = .5)) +
  theme_bw() + labs(x = "Genes Per Group", y = "Youden's Index") + guides(color = F))

#maybe also try the F1 score: 2*((PPV*sens)/(PPV+sens))
(f1.plot.new <- apo.redo.cv %>%
  mutate(f1score = 2*((as.numeric(ppv)/100)*((as.numeric(sensitivity)/100))/
                        ((as.numeric(ppv)/100)+((as.numeric(sensitivity)/100)))))%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(f1score),
            sd_y = sd(f1score),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(y = avg_y, x = as.numeric(ngenes), color = in_out)) + 
  geom_point(position = position_dodge(width = .5)) + 
  geom_smooth(se = F) + 
  geom_errorbar(aes(ymin = avg_y - sd_y, ymax = avg_y + sd_y), position = position_dodge(width = .5)) +
  theme_bw() + labs(x = "Genes Per Group", y = "F1 Score") + guides(color = F)
)


#17 genes
apo.new.pred <- clancPredict(x = x2.new[,!is.na(y$apobec_high)], 
                            classes = y$apobec_high[!is.na(y$apobec_high)],
                             y = x2.new[,!is.na(y$apobec_high)], 17)
apo.new.pred$centroids

#how well does it predict on the given data?
table( true = y$apobec_high[!is.na(y$apobec_high)],pred = apo.new.pred$predictions[,2]) %>% allStats()
#66% correct, 73.8% sensitivity, 63.2% specificity, 42.3% ppv, 86.8% npv

#now, pick the same genes in the mouse models and try again
new.fvg.sub <- fvb.genes[rownames(fvb.genes) %in% toupper(sapply(str_split(rownames(apo.new.pred$centroids), "_"), "[[", 1)),]
new.nsg.sub <- nsg.gense[rownames(nsg.gense) %in% toupper(sapply(str_split(rownames(apo.new.pred$centroids), "_"), "[[", 1)),]

#missing 4 of the genes, which ones?

newgenelist <- toupper(sapply(str_split(rownames(apo.new.pred$centroids), "_"), "[[", 1))

newgenelist[!newgenelist %in% rownames(new.fvg.sub)]
#"C17ORF70" "RAD51L3"  "SPINLW1"  "POTEE"
#C17ORF70 should be "FAAP100", RAD51L3 should be "RAD51D", SPINLW1 should be EPPIN, POTEE should be "ANKRD26"
new.fvg.sub <- fvb.genes[rownames(fvb.genes) %in% c(toupper(sapply(str_split(rownames(apo.new.pred$centroids), "_"), "[[", 1)),
                                                    "FAAP100", "RAD51D", "EPPIN", "ANKRD26"),]
new.nsg.sub <- nsg.gense[rownames(nsg.gense) %in% c(toupper(sapply(str_split(rownames(apo.new.pred$centroids), "_"), "[[", 1)),
                                                    "FAAP100", "RAD51D", "EPPIN", "ANKRD26"),]

#now rename everything so it has the same names as what is in x2.new
g <- sort(rownames(apo.new.pred$centroids))

rownames(new.fvg.sub) <- rownames(new.nsg.sub) <- c(g[25], g[2], g[2], g[3], g[5], g[6], g[7], g[8:14], g[32], g[4],
                                                    g[15:24], g[26:31], g[33], g[34])

nsg.pred <- clancPredict(x = x2.new[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)],
                         y = new.nsg.sub,
                         ngenes = 17)
nsg.pred$predictions


fvb.pred <- clancPredict(x = x2.new[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)],
                         y = new.fvg.sub,
                         ngenes = 17)
fvb.pred$predictions





