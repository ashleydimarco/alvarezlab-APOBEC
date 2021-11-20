#################################################
# Sarah C. Van Alsten                           #
# Date: January 30, 2020                        #
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

x2 <- x2 %>% medianCtr()

#write it out, then impute with readarray
write.table(x2, "data/alltcgalogmed.txt", sep = "\t", col.names = NA)

x3 <- readarray("data/alltcgalogmed.txt")
write.table(x3$xd, "data/alltcgalogmedimpute.txt")

#find the most variable genes
geneList <- rownames(x3$xd)

A2b <<- read.table("data/alltcgalogmedimpute.txt")

#if we want UNADJUSTED scores output instead, use this function instead
#rows argument denotes the row numbers of the training set in which to select the most variable genes
fullUnadjFun <- function(baseString, rown, geneList1 = geneList, data1 = y){
  
  A2 <- A2b[,rown]
  
  #extract apobec calls
  apobec.sig <- data1$apobec_high
  apobec.sig <- apobec.sig[rown]
  
  #can't have missing values for the outcome variables in SAM.
  #make subsets that don't include missings
  apobec.A <- A2[,!is.na(apobec.sig)]
  
  print(dim(apobec.A))
  print(length(apobec.sig))
  print(length(apobec.sig[!is.na(apobec.sig)]))
  
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
#signature is which of the mutation signatures to evaluate
runTuneTrain <- function(i, rown, x1 = geneSub, extra = NULL, signature = "apobec_high"){
  
  
  y.new <- y[rown,]
  
  #depending on what signature is given, subset x1/x3 based into stratified training and test set
  if (signature == "Aging"){
    y.new$samp <- ifelse(rownames(y.new) %in% which(y.new$Aging == 2), 1, 0)
  } else if (signature == "UV"){
    y.new$samp <- ifelse(rownames(y.new) %in% which(y.new$UV == 2), 1, 0)
  } else if (signature == "HR"){
    y.new$samp <-ifelse(rownames(y.new) %in% which(y.new$HR == 2), 1, 0)
  } else if (signature == "APOBEC"){
    y.new$samp <- ifelse(rownames(y.new) %in% which(y.new$APOBEC == 2), 1, 0)
  } else if (signature == "APOBEC3A3B"){
    y.new$samp <- ifelse(rownames(y.new) %in% which(y.new$APOBEC3A3B == 2), 1, 0)
  } else if (signature == "apobec_high"){
    y.new$samp <- ifelse(rownames(y.new) %in% which(y.new$apobec_high == 2), 1, 0)
  }
  
  train <- y.new %>%
    mutate(rown = row_number())%>%
    group_by(samp)%>%
    slice_sample(prop = .7)
  test <- y.new %>%
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



di <- di %>% filter(sample_id %in% substr(y$Barcode.x, 1, 15))
y <- y %>% filter(substr(Barcode.x, 1, 15) %in% di$sample_id)
y$Barcode.x <- substr(y$Barcode.x, 1, 15)

y <- di %>%
  dplyr::select(sample_id, apobec_high) %>%
  full_join(y, by = c("sample_id" = "Barcode.x"))

table(y$Apo_targeted, y$apobec_high)
table(y$APOBEC_wes, y$apobec_high)
table(y$apo_expr, y$apobec_high)


#needs to be in 1/2 format for samr selection
y$apobec_high <- ifelse(y$apobec_high == 1, 2, 1)

#define the 10 folds for cross validation
set.seed(29047)
y$sample_fold <- sample(1:nrow(y),size = nrow(y),replace = F)
y$sample_fold <- factor(as.numeric(gtools::quantcut(y$sample_fold, q = 10)))



# Iteration 1 of CV -------------------------------------------------------

#leave fold 1 out and perform full adj fun
m1 <- fullUnadjFun("fold1", which(!y$sample_fold==1))

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)


#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

y$barcode_short2 <- str_replace_all(substr(y$barcode_short, 1, 16), "-", "\\.")
geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==1)]


full1 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==1),
                                                        x1 = geneSub, signature = "apobec_high"))



# Iteration 2 of CV -------------------------------------------------------

#leave fold 2 out and perform full adj fun
m1 <- fullUnadjFun("fold2", which(!y$sample_fold==2))

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==2)]


full2 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==2),
                                                        x1 = geneSub, signature = "apobec_high"))



# Iteration 3 of CV -------------------------------------------------------

#leave fold 3 out and perform full adj fun
m1 <- fullUnadjFun("fold1", which(!y$sample_fold==3))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==3)]


full3 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==3),
                                                        x1 = geneSub, signature = "apobec_high"))


# Iteration 4 of CV -------------------------------------------------------

#leave fold 4 out and perform full adj fun
m1 <- fullUnadjFun("fold4", which(!y$sample_fold==4))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==4)]


full4 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==4),
                                                        x1 = geneSub, signature = "apobec_high"))




# Iteration 5 of CV -------------------------------------------------------

#leave fold 5 out and perform full adj fun
m1 <- fullUnadjFun("fold5", which(!y$sample_fold==5))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==5)]


full5 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==5),
                                                        x1 = geneSub, signature = "apobec_high"))



# Iteration 6 of CV -------------------------------------------------------

#leave fold 6 out and perform full adj fun
m1 <- fullUnadjFun("fold6", which(!y$sample_fold==6))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==6)]


full6 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==6),
                                                        x1 = geneSub, signature = "apobec_high"))




# Iteration 7 of CV -------------------------------------------------------

#leave fold 7 out and perform full adj fun
m1 <- fullUnadjFun("fold7", which(!y$sample_fold==7))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==7)]


full7 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==7),
                                                        x1 = geneSub, signature = "apobec_high"))




# Iteration 8 of CV -------------------------------------------------------

#leave fold 8 out and perform full adj fun
m1 <- fullUnadjFun("fold8", which(!y$sample_fold==8))

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==8)]


full8 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==8),
                                                        x1 = geneSub, signature = "apobec_high"))




# Iteration 9 of CV -------------------------------------------------------

#leave fold 9 out and perform full adj fun
m1 <- fullUnadjFun("fold9", which(!y$sample_fold==9))


#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==9)]


full9 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==9),
                                                        x1 = geneSub, signature = "apobec_high"))




# Iteration 10 of CV ------------------------------------------------------

#leave fold 10 out and perform full adj fun
m1 <- fullUnadjFun("fold10", which(!y$sample_fold==10))

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]

geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]
geneSub <- geneSub[,which(!y$sample_fold==10)]


full10 <-  purrr::map_dfr(.x = 1:513, .f = ~runTuneTrain(i = .x, rown = which(!y$sample_fold==10),
                                                         x1 = geneSub, signature = "apobec_high"))



#attach all 10 of the data sets into one
full <- full1 %>%
  bind_rows(full2) %>%
  bind_rows(full3) %>%
  bind_rows(full4) %>%
  bind_rows(full5) %>%
  bind_rows(full6) %>%
  bind_rows(full7) %>%
  bind_rows(full8) %>%
  bind_rows(full9) %>%
  bind_rows(full10)



full <- full %>%
  mutate(youden = as.numeric(sensitivity)+ as.numeric(specificity)-1)


# Plotting ----------------

plotCV <- function(d, signature){
  d2a <-  d %>%
    filter(mod_type=="full" & in_out == "in_sample")%>%
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
    geom_point(position = position_dodge(width = .5), size = .5, alpha = .8)+
    geom_errorbar(aes(ymin = (avg_sens - sd_sens), ymax = (avg_sens + sd_sens)), position = position_dodge(width = .5),
                  alpha = .1) +
    theme_bw() +
    #ggtitle(paste0("10-Fold Cross Validation of Sensitivity for ", signature)) +
    labs(color = "Prediction Type", x = "Genes Per Group", y = "Sensitivity (%)") +
    geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed", span = 100) +
    guides(fill = F) + ylim(c(50,75))
  
  specPlot <- d2 %>%
    ggplot(aes(x = ngenes, y = avg_spec, color = in_out)) +
    geom_point(position = position_dodge(width = .5), size = .5, alpha = .8)+
    geom_errorbar(aes(ymin = (avg_spec - sd_spec), ymax = (avg_spec + sd_spec)), position = position_dodge(width = .5),
                  alpha = .1) +
    theme_bw() +
    #ggtitle(paste0("10-Fold Cross Validation of Specificity for ", signature)) +
    labs(color = "Prediction Type", x = "Genes Per Group", y = "Specificity (%)") +
    geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed", span = 100) +
    guides(fill = F)+ ylim(c(50,75))
  
  return(list(sensPlot, specPlot))
}



plots <- plotCV(full)

summarized <- full %>%
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

summarized2 <- full%>%
  mutate(youden = as.numeric(sensitivity)/100+as.numeric(specificity)/100-1)%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(youden),
            sd_y = sd(youden),
            .groups = "keep")

summarized2 %>%
  group_by(in_out) %>%
  arrange(desc(avg_y)) %>%
  slice(1:10)

#13 comes up for both, and is better than 5 actually looks pretty good for sensitivity here



# Final Evaluations -------------------------------------------------------

#most variable genes in entire dataset:

#leave fold 10 out and perform full adj fun
m1 <- fullUnadjFun("foldALL", which(y$sample_fold %in% 1:10))

#now, want to get the most variably expressed genes
m1.sub <- m1 %>%
  arrange(desc(fold_change)) %>%
  #arrange(p_value) %>%
  slice_head(prop = .05)

#based on the genes in this subsetted list, get a frame of TCGA genes which I can then feed into Clanc
geneSub <- x3$xd[rownames(x3$xd) %in% m1.sub$gene_name,]
geneSub <- geneSub[,which(colnames(geneSub) %in% y$barcode_short2)]

#see how  the 13 gene predictor does

apo13 <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], classes = y$apobec_high[!is.na(y$apobec_high)], y = geneSub, ngenes = 13)
write.table(apo13$centroids, "data/dimarco/apobec_enrichment_gene_predictor_13.txt", sep = "\t")

rownames(apo13$centroids)
table(true = y$apobec_high, pred = apo13$predictions[,2])
table(true = y$apobec_high, pred = apo13$predictions[,2]) %>% allStats()

#make a plot of the youden's index
#also plot the Youden's index, which could be a good addition to the fourth part of the figure
youden.plot <- full %>%
  mutate(youden = as.numeric(sensitivity)/100+as.numeric(specificity)/100-1)%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(youden),
            sd_y = sd(youden),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(x = as.numeric(ngenes), y = avg_y, color = in_out)) +
  geom_point(position = position_dodge(width = .5), size = .5, alpha = .8)+
  geom_errorbar(aes(ymin = (avg_y - sd_y), ymax = (avg_y + sd_y)), position = position_dodge(width = .5),
                alpha = .1) +
  theme_bw() +
  #ggtitle(paste0("10-Fold Cross Validation of Sensitivity for ", signature)) +
  labs(x = "Genes Per Group", y = "Youden's Index") +
  geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed", span = 100) +
  guides(fill = F, color = F) 


#maybe also try the F1 score: 2*((PPV*sens)/(PPV+sens))
f1.plot <- full %>%
  mutate(f1score = 2*((as.numeric(ppv)/100)*((as.numeric(sensitivity)/100))/
                        ((as.numeric(ppv)/100)+((as.numeric(sensitivity)/100)))))%>%
  group_by(ngenes, in_out)%>%
  summarise(avg_y = mean(f1score),
            sd_y = sd(f1score),
            .groups = "keep")%>%
  mutate(in_out = case_when(in_out == "in_sample" ~ "Training Set",
                            T ~ "Test Set"))%>%
  ggplot(aes(x = as.numeric(ngenes), y = avg_y, color = in_out)) +
  geom_point(position = position_dodge(width = .5), size = .5, alpha = .8)+
  geom_errorbar(aes(ymin = (avg_y - sd_y), ymax = (avg_y + sd_y)), position = position_dodge(width = .5),
                alpha = .1) +
  theme_bw() +
  #ggtitle(paste0("10-Fold Cross Validation of Sensitivity for ", signature)) +
  labs(x = "Genes Per Group", y = "Youden's Index") +
  geom_smooth(alpha = .1, aes(color = in_out, fill = in_out), se = F, linetype = "dashed", span = 100) +
  guides(fill = F, color = F) 

#13 genes looks like it does pretty well. Make a figure showing sensitivity, specificity, f1 score and then a confusion matrix



#make a plottable confusion matrix
conf <- data.frame(table(Expression = apo13$predictions[,2], WES = y$apobec_high))

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
(plots[[1]] + plots[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + confMatrix) + plot_annotation(tag_levels = "a")
dev.off()



# Updated RNASeq Data for Mice -----------------------------------------------------

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
fvb.sub <- fvb.genes[rownames(fvb.genes) %in% toupper(sapply(str_split(rownames(apo13$centroids), "\\|"), "[[",1)),]
nsg.sub <- nsg.gense[rownames(nsg.gense) %in% toupper(sapply(str_split(rownames(apo13$centroids), "\\|"), "[[",1)),]

#some of the genes are missing. Which ones?
toupper(sapply(str_split(rownames(apo13$centroids), "\\|"), "[[",1))[
  !toupper(sapply(str_split(rownames(apo13$centroids), "\\|"), "[[",1)) %in% rownames(nsg.sub)]

#"LOC100131691" "ZNF703"       "ZNF425"       "VNN2"


#now rename to have the same names as in the human version
genelist <- sort(rownames(apo13$centroids))

#rename rows to correspond exactly to the genelsit
rownames(nsg.sub) <- rownames(fvb.sub) <-c(genelist[1:9], genelist[11:23])


for(i in 1:length(genelist)){
  print(noquote(genelist[i]))
}


#can now use the classifier to predict the apobec_high nature of the samples
fvb.pred <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)], y = fvb.sub,
                         ngenes = 13)
nsg.pred <- clancPredict(x = geneSub[,!is.na(y$apobec_high)], 
                         classes = y$apobec_high[!is.na(y$apobec_high)], y = nsg.sub,
                         ngenes = 13)
fvb.pred$predictions
nsg.pred$predictions


#finally, add the mouse prediction data to the figure
#maybe also try to add the mouse validation to the same figure?

#just do for FVB

new.df <- data.frame(WES = c(1,2,1,2),
                     Expression = c(1,1,2,2),
                     Freq = c(5,1,3,3))

#only in the NSG mice

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
pdf("data/dimarco/tcga_classifier_results_v4.pdf",width = 8.5, height = 10)
(plots[[1]] + plots[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + cf1) + plot_annotation(tag_levels = "a")
dev.off()

tiff(file = "data/dimarco/tcga_classifier_results_v4.tiff", width = 6600, height = 7000, units = "px", res = 800)
(plots[[1]] + plots[[2]] + plot_layout(guides = "collect"))/
  (youden.plot + cf1) + plot_annotation(tag_levels = "a")
dev.off()




# Forest Plots ------------------------------------------------------------


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



