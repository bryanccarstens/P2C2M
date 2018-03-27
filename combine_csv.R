##### Combines all output CSVs from P2C2M #####

# Get csv names for each model
csv_list <- list.files(pattern = "\\.csv")

model1 <- csv_list[grep("model1+", csv_list)]
model2 <- csv_list[grep("model2+", csv_list)]
model3 <- csv_list[grep("model3+", csv_list)]
model4 <- csv_list[grep("model4+", csv_list)]

# Model1
model1_reps <- lapply(model1, function(x) strsplit(x, split = "_")[[1]][2]) # get rep numbers
model1_s <- lapply(model1, function(x) strsplit(x, split = "_")[[1]][3]) # get seg sites value
model1_rand <- lapply(model1, function(x) strsplit(strsplit(x, split = "_")[[1]][5], split = "\\.")[[1]][1]) # get randomization value

model1_df <- lapply(model1, function(x) read.csv(x, row.names=1, stringsAsFactors = FALSE)) # read in all csv files for model1

for (df in 1:length(model1_df)){ 
  new_cols <- as.data.frame(cbind(rep("model1", nrow(model1_df[[df]])), rep(model1_rand[[df]][1], nrow(model1_df[[df]])), rep(model1_s[[df]][1], nrow(model1_df[[df]])), rep(model1_reps[[df]][1], nrow(model1_df[[df]]))), stringsAsFactors = FALSE) # make columns for model, randomization, rep, and segsites values
  colnames(new_cols) <- c("Model", "Rand", "S", "Rep")
  model1_df[[df]] <- cbind(new_cols, model1_df[[df]]) # add new columns to existing dataframe
}


# Model2
model2_reps <- lapply(model2, function(x) strsplit(x, split = "_")[[1]][2])
model2_s <- lapply(model2, function(x) strsplit(x, split = "_")[[1]][3])
model2_rand <- lapply(model2, function(x) strsplit(strsplit(x, split = "_")[[1]][5], split = "\\.")[[1]][1])

model2_df <- lapply(model2, function(x) read.csv(x, row.names=1, stringsAsFactors = FALSE))

for (df in 1:length(model2_df)){
  new_cols <- as.data.frame(cbind(rep("model2", nrow(model2_df[[df]])), rep(model2_rand[[df]][1], nrow(model2_df[[df]])), rep(model2_s[[df]][1], nrow(model2_df[[df]])), rep(model2_reps[[df]][1], nrow(model2_df[[df]]))), stringsAsFactors = FALSE)
  colnames(new_cols) <- c("Model", "Rand", "S", "Rep")
  model2_df[[df]] <- cbind(new_cols, model2_df[[df]])
}


# Model3
model3_reps <- lapply(model3, function(x) strsplit(x, split = "_")[[1]][2])
model3_s <- lapply(model3, function(x) strsplit(x, split = "_")[[1]][3])
model3_rand <- lapply(model3, function(x) strsplit(strsplit(x, split = "_")[[1]][5], split = "\\.")[[1]][1])

model3_df <- lapply(model3, function(x) read.csv(x, row.names=1, stringsAsFactors = FALSE))

for (df in 1:length(model3_df)){
  new_cols <- as.data.frame(cbind(rep("model3", nrow(model3_df[[df]])), rep(model3_rand[[df]][1], nrow(model3_df[[df]])), rep(model3_s[[df]][1], nrow(model3_df[[df]])), rep(model3_reps[[df]][1], nrow(model3_df[[df]]))), stringsAsFactors = FALSE)
  colnames(new_cols) <- c("Model", "Rand", "S", "Rep")
  model3_df[[df]] <- cbind(new_cols, model3_df[[df]])
}



# Model4
model4_reps <- lapply(model4, function(x) strsplit(x, split = "_")[[1]][2])
model4_s <- lapply(model4, function(x) strsplit(x, split = "_")[[1]][3])
model4_rand <- lapply(model4, function(x) strsplit(strsplit(x, split = "_")[[1]][5], split = "\\.")[[1]][1])

model4_df <- lapply(model4, function(x) read.csv(x, row.names=1, stringsAsFactors = FALSE))

for (df in 1:length(model4_df)){
  new_cols <- as.data.frame(cbind(rep("model4", nrow(model4_df[[df]])), rep(model4_rand[[df]][1], nrow(model4_df[[df]])), rep(model4_s[[df]][1], nrow(model4_df[[df]])), rep(model4_reps[[df]][1], nrow(model4_df[[df]]))), stringsAsFactors = FALSE)
  colnames(new_cols) <- c("Model", "Rand", "S", "Rep")
  model4_df[[df]] <- cbind(new_cols, model4_df[[df]])
}

# Combine all reps for each model 
out_model1 <- as.data.frame(do.call(rbind, model1_df))
out_model2 <- as.data.frame(do.call(rbind, model2_df))
out_model3 <- as.data.frame(do.call(rbind, model3_df))
out_model4 <- as.data.frame(do.call(rbind, model4_df))

# combine all models and output new csv
out_df <- as.data.frame(rbind(out_model1, out_model2, out_model3, out_model4))
write.csv(out_df, file = "out.csv", row.names = FALSE)





