#produce cosine similarity matrix
#Figure S1

library(reshape2)
library(lsa)
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

combined_data <- read.csv("./data/long_combined_data.csv", sep=";", header = T)

################################
#DATA PREPARATION
################################
#Modify from long to wide format
combined_data$newid2 <- paste(combined_data$person_id, combined_data$month)
CLINICAL <- combined_data[combined_data$type=="clinical",]
METABOLOMICS <- combined_data[combined_data$type=="metabolomics",]

CLINICAL_wide <- dcast(CLINICAL, newid2 ~ label, value.var = "value", fun.aggregate = mean)
CLINICAL_wide$newid2 <- paste(CLINICAL_wide$newid2, "clin", sep = "_")
CLINICAL_wide <- CLINICAL_wide[, c(1,2,4:12,3)]
colnames(CLINICAL_wide) <- c("newid2", "Albumin", "ApoA1", "ApoB", "Total_C", "Creatinine", "Glucose", "HDL_C", "Clinical_LDL_C", 
                             "non_HDL_C", "Total_TG", "ApoB_by_ApoA1")

METABOLOMICS_wide <- dcast(METABOLOMICS, newid2 ~ label, value.var = "value", fun.aggregate = mean)
METABOLOMICS_wide$newid2 <- paste(METABOLOMICS_wide$newid2, "met", sep = "_")
METABOLOMICS_wide <- METABOLOMICS_wide[, c(1:4,11,7,8,9,6,10, 12,5)]

#scale to [0, 1]
# x' = (x - min(x))/(max(x) - min(x))
#combine two datasets
ALL <- rbind.data.frame(CLINICAL_wide, METABOLOMICS_wide)
ALL_scaled <- ALL
for (i in 2:length(ALL)) {
  for (j in 1:length(ALL$newid)) {
    ALL_scaled[j, i] = (ALL[j, i] - min(ALL[, i]))/(max(ALL[, i]) - min(ALL[, i]))   
  }
}

CLINICAL_wide_scaled <- ALL_scaled[1:180, ]
METABOLOMICS_wide_scaled <- ALL_scaled[181:360, ]

#transpose and calculate cosine similarity
trans_CLINICAL_wide_scaled <- as.data.frame(t(CLINICAL_wide_scaled[, -1]))
colnames(trans_CLINICAL_wide_scaled) <- CLINICAL_wide_scaled$newid

trans_METABOLOMICS_wide_scaled <- as.data.frame(t(METABOLOMICS_wide_scaled[, -1]))
colnames(trans_METABOLOMICS_wide_scaled) <- METABOLOMICS_wide_scaled$newid

trans_scaled_combined <- cbind.data.frame(trans_CLINICAL_wide_scaled, trans_METABOLOMICS_wide_scaled)
trans_scaled_combined <- trans_scaled_combined[, order(names(trans_scaled_combined))]

trans_scaled_combined2 <- as.matrix(trans_scaled_combined)
cosine_scaled_combined <- cosine(trans_scaled_combined2)
cosine_scaled_combined <- as.data.frame(cosine_scaled_combined)

#remove unnecessary columns and rows
cosine_colnames <- colnames(cosine_scaled_combined)
cosnames <- matrix(NA)
for (i in 1:length(cosine_colnames)) {
  if(grepl("met", cosine_colnames[i])){
    cosnames[i] <- i
  }
}
cosnames <- subset(cosnames, (!is.na(cosnames)))
cosine_scaled_combined <- cosine_scaled_combined[, -c(cosnames)]


cosine_scaled_combined$rownames <- rownames(cosine_scaled_combined)
for (i in 1:length(cosine_scaled_combined$`IAM01 1_clin`)) {
  if(grepl("clin", cosine_scaled_combined$rownames[i])){
    cosine_scaled_combined <- cosine_scaled_combined[-i, ]
  }
}
cosine_scaled_combined <- cosine_scaled_combined[, -181]

#Figure S1
data=as.matrix(cosine_scaled_combined)
col_fun1 = colorRamp2(c(0.4, 0.8, 1), c("#34A3DC", "white", "#F5C000"))
Heatmap(data, name = "Cosine similarity",show_row_names = FALSE, show_column_names = FALSE,
        row_order = order(rownames(data)), column_order = order(colnames(data)), col = col_fun1, 
        column_title = "Clinical", row_title = "Metabolomics") 
