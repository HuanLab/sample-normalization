#Closing the knowledge gap of post-acquisition sample normalization for untargeted metabolomics

#This R script that simulates feature tables which can then be used to evaluate post-acquisition
#sample normalization methods

#Brian Low
#June 14, 2024

################################################################################

#Load packages

library("MAFFIN")
library("doParallel")
library("foreach")
library("doRNG")

message("Finished loading packages.")

#Load data; we need to use some starting values to get experimental intensity and
#missing value distributions to be used in the simulations

setwd("C:/Users/User/Desktop/Brian/Normalization/20221202_TPR_Spider_Iterated/20240614_manuscript")
input = read.csv("initial_values.csv")

missing_dist = input$missing_dist
missing_dist = missing_dist[is.finite(missing_dist)]

#Set parameters

unchanged = 2500 #Number of non-dysregulated features
dys = 5000 - unchanged #Number of dysregulated features
up_regulated = c(0.5,0.75) #Percentage of dysregulated features that are upregulated
total_missing = 0.3 #Percentage of injected missing values
sample_imbalance = T #TRUE will introduce random sample variation (i.e., a percentage
#of features between each sample are dysregulated)
replicates = 10 #Number of samples in each experimental group
iterations = 100 #Number of iterations for simulation and benchmarking
output = F #True to save results as a .csv file

################################################################################

#Set up backend

options(warn = 2)

numCores = detectCores()
cl = makeCluster(numCores - 1)
registerDoParallel(cl)

if(getDoParWorkers() != numCores - 1){
  stop("Parallel processing not registered!")
} else {
  message(paste(numCores, "cores detected. Using", numCores - 1, "cores for simulations..."))
}

clusterEvalQ(cl,library("MAFFIN"))

################################################################################

#Functions used to normalize and imputate the simulated feature tables

#Formatting to allow compatability with the MAFFIN package 
MAFFIN_compatible = function(data){
  
  WT_names = rep("WT", replicates)
  KO_names = rep("KO", replicates)
  data = rbind(c(WT_names, KO_names), data)
  identifier = seq(0, (nrow(data)-1),1)
  data = cbind(identifier, data)
  
  return(as.data.frame(data))
  
}

#Sum method
sumNorm_df = function(data){
  
  ref_index = which.min(colSums(data == 0))
  
  working_df = MAFFIN_compatible(data)
  
  #Normalize
  
  estimated_k = suppressMessages(SumNorm(working_df, RunEvaluation = F)$NormFactor)
  estimated_k = estimated_k/estimated_k[ref_index]
  
  for(i in 2:ncol(working_df)){
    working_df[2:nrow(working_df),i] = as.numeric(working_df[2:nrow(working_df),i])/
      estimated_k[i-1]
  }
  
  #Missing value imputation
  
  for(i in 2:nrow(working_df)){
    b = as.numeric(working_df[i,2:ncol(working_df)])
    b[b == 0] = (min(b[b > 0]))/5
    working_df[i,2:ncol(working_df)] = b
  }
  
  results = list(working_df, estimated_k)
  
  return(results)
  
}

#Median method
medianNorm_df = function(data){
  
  ref_index = which.min(colSums(data == 0))
  
  working_df = MAFFIN_compatible(data)
  
  #Normalize
  
  estimated_k = suppressMessages(MedianNorm(working_df, RunEvaluation = F)$NormFactor)
  estimated_k = estimated_k/estimated_k[ref_index]
  
  for(i in 2:ncol(working_df)){
    working_df[2:nrow(working_df),i] = as.numeric(working_df[2:nrow(working_df),i])/
      estimated_k[i-1]
  }
  
  #Missing value imputatation
  
  for(i in 2:nrow(working_df)){
    b = as.numeric(working_df[i,2:ncol(working_df)])
    b[b == 0] = (min(b[b > 0]))/5
    working_df[i,2:ncol(working_df)] = b
  }
  
  results = list(working_df, estimated_k)
  
  return(results)
  
}

#PQN method
PQNNorm_df = function(data){
  
  ref_index = which.min(colSums(data == 0))
  
  working_df = MAFFIN_compatible(data)
  estimated_k = c()
  ref = as.numeric(working_df[2:nrow(working_df), ref_index + 1])
  
  #Get normalization factors
  
  for(i in 2:ncol(working_df)){
    fc = as.numeric(working_df[2:nrow(working_df),i])/ref
    fc = fc[is.finite(log2(fc))]
    estimated_k[i-1] = median(fc) 
  }
  
  #Normalize
  
  for(i in 2:ncol(working_df)){
    working_df[2:nrow(working_df),i] = as.numeric(working_df[2:nrow(working_df),i])/
      estimated_k[i-1]
  }
  
  #Missing value imputation
  
  for(i in 2:nrow(working_df)){
    b = as.numeric(working_df[i,2:ncol(working_df)])
    b[b == 0] = (min(b[b > 0]))/5
    working_df[i,2:ncol(working_df)] = b
  }
  
  results = list(working_df, estimated_k)
  
  return(results)
  
}

#MDFC method
MDFCNorm_df = function(data, bw = 0.3){
  
  ref_index = which.min(colSums(data == 0))
  
  working_df = MAFFIN_compatible(data)
  estimated_k = c()
  ref = as.numeric(working_df[2:nrow(working_df), ref_index + 1])
  
  #Get normalization factors
  
  for(i in 2:ncol(working_df)){
    
    fc = log2(as.numeric(working_df[2:nrow(working_df),i])/ref)
    fc = fc[is.finite(fc)]
    d = density(fc, bw = bw)
    
    if(i == ref_index+1){
      estimated_k[i-1] = 1
    } else {
      estimated_k[i-1] = 2^d$x[which.max(d$y)]
    }
    
  }
  
  #Normalize
  
  for(i in 2:ncol(working_df)){
    working_df[2:nrow(working_df),i] = as.numeric(working_df[2:nrow(working_df),i])/
      estimated_k[i-1]
  }
  
  #Missing value imputation
  
  for(i in 2:nrow(working_df)){
    b = as.numeric(working_df[i,2:ncol(working_df)])
    b[b == 0] = (min(b[b > 0]))/5
    working_df[i,2:ncol(working_df)] = b
  }
  
  results = list(working_df, estimated_k)
  
  return(results)
}

#Quantile method
quantile_df = function(data){
  
  #Normalize
  
  working_df = MAFFIN_compatible(data)
  working_df = suppressMessages(QuantileNorm(working_df, RunEvaluation = F)$NormedTable)
  
  #Missing value imputation
  
  for(i in 2:nrow(working_df)){
    b = as.numeric(working_df[i,2:ncol(working_df)])
    b[b == 0] = (min(b[b > 0]))/5
    working_df[i,2:ncol(working_df)] = b
  }
  
  return(working_df)
  
}

#Class-specific quantile method
cquantile_df = function(data){
  
  working_df = MAFFIN_compatible(data)
  
  group_a = grep("WT", working_df[1,])
  group_b = grep("KO", working_df[1,])
  
  df1 = working_df[,c(1,group_a)]
  df2 = working_df[,c(1,group_b)]
  
  #Normalize each subset independently
  
  quantile_1 = suppressMessages(QuantileNorm(df1, 
                                             RunEvaluation = F)$NormedTable)
  quantile_2 = suppressMessages(QuantileNorm(df2, 
                                             RunEvaluation = F)$NormedTable)
  
  #Merge results
  
  class_quantiles = cbind(quantile_1, quantile_2[,-1])
  
  #Missing value imputation
  
  for(i in 2:nrow(class_quantiles)){
    b = as.numeric(class_quantiles[i,2:ncol(class_quantiles)])
    b[b == 0] = (min(b[b > 0]))/5
    class_quantiles[i,2:ncol(class_quantiles)] = b
  }
  
  return(class_quantiles)
}

################################################################################

#Function used to inject missing values into the feature table
simulate_missing = function(data, mcar, mar, mnar, sim_na){
  
  current_na = mean(is.na(data))
  row_na = sample(missing_dist, size = nrow(data), replace = T)
  
  #Remove missing values from MCAR mechanism
  
  for(i in 1:ncol(data)){
    mcar_na = runif(nrow(data), 0, 1)
    temp_sample = as.numeric(data[,i])
    mcar_rm = which(mcar_na < mcar)
    temp_sample[mcar_rm] = NA
    data[,i] = temp_sample
  }
  
  #Update current missing percentage
  
  current_na = mean(is.na(data))
  
  #Remove missing values from MAR mechanism
  
  while(current_na < sum(c(mar, mcar))){
    
    #Choose two different features
    
    feature1 = sample(1:nrow(data), 1)
    feature2 = sample(setdiff(1:nrow(data), feature1), 1)
    
    #Cutoff percent using missing distribution
    
    cutoff_percent = row_na[feature1]
    cutoff_index = round(ncol(data)*cutoff_percent)
    
    #Remove the highest cutoff_index values in feature 2 based on feature 1
    
    sorted_feature1 = order(as.numeric(data[feature1,]), na.last = F)
    
    if(cutoff_index >= 1){
      remove_index = sorted_feature1[(length(sorted_feature1) - 
                                        cutoff_index + 1):length(sorted_feature1)]
      data[feature2, remove_index] = NA
    }
    
    #Update current missing percentage
    
    current_na = sum(colSums(is.na(data)))/(ncol(data)*nrow(data))
    
  }
  
  #Update current missing percentage
  
  current_na = mean(is.na(data))
  
  #Resort the and feature table and missing distribution so that 
  #features with low intensities will more likely be removed

  data$mean_intensity = apply(data, 1, mean, na.rm = T)
  data = data[order(data$mean_intensity, decreasing = F),]
  data$mean_intensity = NULL
  row_na = sort(row_na, decreasing = T)
  
  #Remove features with the MNAR mechanism
  
  while(current_na < sim_na){
    
    #Randomly select a feature and use its corresponding missing percentage
    
    feature1 = sample(1:nrow(data), 1)
    temp_feature = as.numeric(data[feature1, ])
    
    remove_index = sample(c(0,1), size = replicates*2, prob = 
                            c(1-row_na[feature1], row_na[feature1]), replace = T)
    remove_number = sum(remove_index == 1)
    order_index = order(temp_feature)
    
    #Remove lowest intensities
    
    if(remove_number >= 1){
      temp_feature = replace(temp_feature, order_index[1:remove_number], NA)
    }
    
    data[feature1,] = temp_feature
    
    current_na = mean(is.na(data))
    
  }
  
  return(data)
  
}

################################################################################

#Confusion matrix to evaluate results
confusionMatrix = function(ref, data){
  
  #Check row names are the same
  
  if(identical(as.numeric(row.names(ref)), 
               as.numeric(row.names(data))) == F){
    stop("Row names not the same for evaluation!")
  }
  
  group_a = grep("WT", ref[1,])
  group_b = grep("KO", ref[1,])
  
  row_index = as.numeric(row.names(ref))[-1]
  total_rows = nrow(ref) - 1
  
  #For each the reference (i.e., ground truth) and normalized feature table, 
  #extract significant features
  
  ref = ref[-1,]
  ref_sig = c()
  
  ref_numeric = sapply(ref, as.numeric)
  ref_sig = apply(ref_numeric, 1, function(w) wilcox.test(w[group_a], 
                                                          w[group_b], exact = F)$p.value)
  
  ref_sig = p.adjust(ref_sig, method = "fdr")
  ref$fdr = ref_sig
  
  ref = subset(ref, fdr < 0.05)
  ref_index = as.numeric(row.names(ref))
  
  data = data[-1,]
  data_sig = c()
  
  data_numeric = sapply(data, as.numeric)
  data_sig = apply(data_numeric, 1, function(w) wilcox.test(w[group_a], 
                                                            w[group_b], exact = F)$p.value)
  
  data_sig = p.adjust(data_sig, method = "fdr")
  data$fdr = data_sig
  
  data = subset(data, fdr < 0.05)
  data_index = as.numeric(row.names(data))
  
  #Calculate true positives, true negatives, false positives, and 
  #false negatives
  
  TP = length(intersect(ref_index, data_index))
  
  not_sig = setdiff(row_index, ref_index)
  not_sig_input = setdiff(row_index, data_index)
  
  TN = length(intersect(not_sig, not_sig_input))
  
  FP = length(setdiff(data_index, ref_index))
  
  FN = length(setdiff(ref_index, data_index))
  
  #Check rows match
  
  if(sum(c(TP, TN, FP, FN)) != total_rows){
    stop("TP, TN, FP, and FN not equal to total rows!")
  }
  
  #Calculate true positive rate, false positive rate, true negative rate,
  #and false negative rate
  
  TPR = TP/(TP + FN)
  FPR = FP/(FP + TN)
  TNR = TN/(TN + FP)
  FNR = FN/(TP + FN)
  
  results = data.frame(TPR, FPR, TNR, FNR)
  colnames(results) = c("TPR", "FPR", "TNR", "FNR")
  
  return(results)
}

################################################################################

#Set up reference intensities

reference = input$int
reference = reference[!is.na(reference)]
reference = reference[reference != 0]

#Set up a grid which will show the data structures that will be tested

grid = expand.grid(unchanged, dys, up_regulated, total_missing)
colnames(grid) = c("unchanged", "dys_regulated", "prop_up", "total_missing")
grid = subset(grid, unchanged + dys_regulated == 5000)
grid = subset(grid, prop_up + (1 - prop_up) == 1)

################################################################################

#Begin simulation 

begin_time = Sys.time()
message("Starting simulations...")

set.seed(1234)

out = foreach(i = 1:nrow(grid), .combine = rbind) %dorng% {
  
  #Create empty vectors to store results
  
  pre_TPR = c()
  pre_FPR = c()
  
  sum_TPR = c()
  sum_FPR = c()
  
  median_TPR = c()
  median_FPR = c()
  
  PQN_TPR = c()
  PQN_FPR = c()
  
  MDFC_TPR = c()
  MDFC_FPR = c()
  
  quantile_TPR = c()
  quantile_FPR = c()
  
  cquantile_TPR = c()
  cquantile_FPR = c()
  
  for(j in 1:iterations){
    
    feature_table = as.data.frame(matrix(0, nrow = grid[i,1] + grid[i,2], 
                                         ncol = replicates*2)) 
    
    #Initiate by sampling intensities from
    
    sample1 = sample(reference, size = grid[i,1] + grid[i,2], replace = T)
    feature_table[, 1] = sample1
    
    #Create additional samples using Gaussian noise
    
    for(p in 2:ncol(feature_table)){
      
      noise = rnorm(grid[i,1] + grid[i,2], 1, 0.2)
      
      while(min(noise) < 0){
        noise = rnorm(grid[i,1] + grid[i,2], 1, 0.2)
      }
      
      feature_table[,p] = feature_table[,1] * noise
      
    }
    
    #Add more significant sample-to-sample variation by dysregulating a 
    #percentage of features within each sample
    
    if(sample_imbalance){
      
      for(p in 1:ncol(feature_table)){
        
        #Randomly determine how balanced the dysregulation will be
        
        fluctuate_prob = runif(1,0.25,0.75)
        fluctuate_dir = sample(c(0,1), size = round(0.25*(grid[i,1] + grid[i,2])), 
                               prob = c(fluctuate_prob, 1 - fluctuate_prob), replace = T)
        #Select which features to dysregulate
        sample_features = sample(1:nrow(feature_table), size = length(fluctuate_dir))
        
        #Dysregulate by scaling by a factor
        bio_var = rep(0, length(sample_features))
        bio_var[fluctuate_dir == 1] = runif(sum(fluctuate_dir == 1), 1.5, 10)
        bio_var[fluctuate_dir == 0] = runif(sum(fluctuate_dir == 0), 0.1, 0.67)
        
        feature_table[sample_features,p] = feature_table[sample_features,p] *
          bio_var
        
      }
      
    }
    
    #Scale each sample by a dilution factor and save the factor
    
    conc = sample(c(0,1), size = (replicates*2) - 1, replace = T)
    
    actual_k = c()
    actual_k[1] = 1
    
    for(p in 2:ncol(feature_table)){
      if(conc[p-1] == 0){
        k = runif(1, 1, 4.5)
        feature_table[,p] = feature_table[,p] * k
      } else {
        k = runif(1, 1/4.5, 1)
        feature_table[,p] = feature_table[,p] * k
      }
      actual_k[p] = k
    }
    
    #Dysregulate a set of features in group 2; this will represent
    #features perturbed from experimental treatment
    
    dir = sample(c(0,1), size = grid[i,2], replace = T, prob = c(grid[i,3], (1 - grid[i,3])))
    fc_scale = rep(0, length(dir))
    
    if(length(dir) != 0){
      fc_scale[dir == 0] = runif(sum(dir == 0), 1.5, 10)
      fc_scale[dir == 1] = runif(sum(dir == 1), 0.1, 0.67)
    }
    
    if(length(dir) != 0){
      
      for(p in (replicates + 1):ncol(feature_table)){
        if(grid[i,1] == 0){
          unchanged_metabolites = NULL
        } else {
          unchanged_metabolites = (feature_table[,p])[1:grid[i,1]]
        }
        
        dys_metabolites = (feature_table[,p])[(grid[i,1]+1):nrow(feature_table)]
        dys_metabolites = dys_metabolites * fc_scale
        
        feature_table[,p] = c(unchanged_metabolites, dys_metabolites)
        
      }
      
    }
    
    #Inject missing values
    
    feature_table = simulate_missing(feature_table, mcar = grid[i,4]*0.1, 
                                     mar = grid[i,4]*0.1, 
                                     mnar = grid[i,4]*0.8, sim_na = grid[i,4])
    
    #Check if total missingness in feature table makes sense
    
    if(round(mean(is.na(feature_table)), 1) != grid[i,4]){
      stop("Set and simualted percent missingness do not match!")
    }
    
    #Remove features with too many missing values
    
    count_na = rowSums(is.na(feature_table))
    omit_index = count_na > replicates
    feature_table = feature_table[!omit_index,]
    
    #Find the reference sample: the sample with the least number of missing values
    #Then rescale true dilution factors using the reference
    
    normalize_by = which.min(colSums(is.na(feature_table)))
    rescaled_k = actual_k/actual_k[normalize_by]
    
    #Replace missing values with 0s
    
    feature_table[is.na(feature_table)] = 0
    
    #Reference table (i.e., ground truth)
    
    reference_table = feature_table
    
    #Normalize reference feature table with the true and rescaled dilution factors
    
    reference_table = sweep(reference_table, 2, rescaled_k, "/")
    
    #Missing value imputation
    
    for(p in 1:nrow(reference_table)){
      temp = as.numeric(reference_table[p,])
      temp[temp == 0] = min(temp[temp > 0])/5
      reference_table[p,] = temp
    }
    
    #Reformat reference feature table for benchmarking
    
    reference_table = MAFFIN_compatible(reference_table)
    
    #Evaluate! For each method, normalize and calculate TPR and FPR
    #Pre-norm
    
    raw_df = feature_table
    
    for(p in 1:nrow(raw_df)){
      temp = as.numeric(raw_df[p,])
      temp[temp == 0] = min(temp[temp > 0])/5
      raw_df[p,] = temp
    }
    
    raw_df = MAFFIN_compatible(raw_df)
    
    pre_results = confusionMatrix(ref = reference_table, data = raw_df)
    
    pre_TPR[j] = pre_results[[1]]
    pre_FPR[j] = pre_results[[2]]
    
    #Sum normalization
    
    sum_norm = sumNorm_df(feature_table)
    sum_results = confusionMatrix(ref = reference_table, data = 
                                    sum_norm[[1]])
    
    sum_TPR[j] = sum_results[[1]]
    sum_FPR[j] = sum_results[[2]]
    
    #Median normalization
    
    median_norm = medianNorm_df(feature_table) 
    median_results = confusionMatrix(ref = reference_table, data = 
                                       median_norm[[1]])
    
    median_TPR[j] = median_results[[1]]
    median_FPR[j] = median_results[[2]]
    
    #PQN
    
    PQN_norm = PQNNorm_df(feature_table)
    PQN_results = confusionMatrix(ref = reference_table, data = 
                                    PQN_norm[[1]])
    
    PQN_TPR[j] = PQN_results[[1]]
    PQN_FPR[j] = PQN_results[[2]]
    
    #MDFC normalization
    
    MDFC_norm = MDFCNorm_df(feature_table)
    MDFC_results = confusionMatrix(ref = reference_table, data = 
                                     MDFC_norm[[1]])
    
    MDFC_TPR[j] = MDFC_results[[1]]
    MDFC_FPR[j] = MDFC_results[[2]]
    
    #Quantile normalization
    
    quantile_results = confusionMatrix(ref = reference_table, 
                                       data = quantile_df(feature_table))
    
    quantile_TPR[j] = quantile_results[[1]]
    quantile_FPR[j] = quantile_results[[2]]
    
    #Class-specifc quantile normalization
    
    cquantile_results = confusionMatrix(ref = reference_table, data = 
                                          cquantile_df(feature_table))
    
    cquantile_TPR[j] = cquantile_results[[1]]
    cquantile_FPR[j] = cquantile_results[[2]]
    
  }
  
  #At the end, average TPR and FPR from all iterations
  
  results = c(mean(pre_TPR), mean(pre_FPR),
              mean(sum_TPR), mean(sum_FPR),
              mean(median_TPR), mean(median_FPR),
              mean(PQN_TPR), mean(PQN_FPR),
              mean(MDFC_TPR), mean(MDFC_FPR),
              mean(quantile_TPR), mean(quantile_FPR),
              mean(cquantile_TPR), mean(cquantile_FPR))
  
  return(results)
  
}

#Clean up

finish_time = Sys.time()
message(paste("Done. Simulations took", 
              round(as.numeric(difftime(finish_time, begin_time, units = "min")), 2), "min to finish."))
stopCluster(cl)

################################################################################

#Format results table and save

out = as.data.frame(out)
out = cbind(grid,out)
colnames(out)[5:ncol(out)] = c("pre_TPR", "pre_FPR",
                               "sum_TPR", "sum_FPR", 
                               "median_TPR", "median_FPR",
                               "PQN_TPR", "PQN_FPR", 
                               "MDFC_TPR", "MDFC_FPR",
                               "quantile_TPR", "quantile_FPR", 
                               "cquantile_TPR", "cquantile_FPR")

if(output){
  write.csv(out, "norm_benchmarking_results.csv", row.names = F)
}
