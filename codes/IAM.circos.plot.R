#produce the circos plot
#Figure 3

library(circlize)
library(tidyr)
library(caret)
library(ComplexHeatmap)
library(gridBase)
# Load the dataset
data <- read.csv("./data/long_combined_data.csv", sep=";", header = T)
# Set 'label' and 'type' columns as characters to change the values. They are formatted as factors so it is not possible to change them with their defined levels.
data$label = as.character(data$label)
data$type = as.character(data$type)

# Change the names of the parameters to an appropriate format.
data = within(data, label[label == "Albumin" & type == "metabolomics"] <- as.character("Albumin:Metab"))
data = within(data, label[label == "Total_C" & type == "metabolomics"] <- as.character("Total Cholesterol:Metab"))
data = within(data, label[label == "HDL_C" & type == "metabolomics"] <- as.character("HDL Cholesterol:Metab"))
data = within(data, label[label == "Clinical_LDL_C" & type == "metabolomics"] <- as.character("LDL Cholesterol:Metab"))
data = within(data, label[label == "non_HDL_C" & type == "metabolomics"] <- as.character("Non-HDL Cholesterol:Metab"))
data = within(data, label[label == "Glucose" & type == "metabolomics"] <- as.character("Glucose:Metab"))
data = within(data, label[label == "Total_TG" & type == "metabolomics"] <- as.character("Triglyceride:Metab"))
data = within(data, label[label == "Creatinine" & type == "metabolomics"] <- as.character("Creatinine:Metab"))
data = within(data, label[label == "ApoA1" & type == "metabolomics"] <- as.character("Apolipoprotein A1:Metab"))
data = within(data, label[label == "ApoB" & type == "metabolomics"] <- as.character("Apolipoprotein B:Metab"))
data = within(data, label[label == "ApoB_by_ApoA1" & type == "metabolomics"] <- as.character("Apolipoprotein B/A1:Metab"))

data = within(data, label[label == "Albumine" & type == "clinical"] <- as.character("Albumin:Clinical"))
data = within(data, label[label == "Cholesterol" & type == "clinical"] <- as.character("Total Cholesterol:Clinical"))
data = within(data, label[label == "HDL-cholesterol" & type == "clinical"] <- as.character("HDL Cholesterol:Clinical"))
data = within(data, label[label == "LDL-cholesterol" & type == "clinical"] <- as.character("LDL Cholesterol:Clinical"))
data = within(data, label[label == "non-HDL cholesterol" & type == "clinical"] <- as.character("Non-HDL Cholesterol:Clinical"))
data = within(data, label[label == "Glucose" & type == "clinical"] <- as.character("Glucose:Clinical"))
data = within(data, label[label == "Triglyceriden" & type == "clinical"] <- as.character("Triglyceride:Clinical"))
data = within(data, label[label == "Creatinine" & type == "clinical"] <- as.character("Creatinine:Clinical"))
data = within(data, label[label == "Apolipoproteine A1" & type == "clinical"] <- as.character("Apolipoproteine A1:Clinical"))
data = within(data, label[label == "Apolipoproteine B" & type == "clinical"] <- as.character("Apolipoproteine B:Clinical"))
data = within(data, label[label == "ApoBbyA1" & type == "clinical"] <- as.character("Apolipoprotein B/A1:Clinical"))

# Set the columns again as factors.
data$label = as.factor(data$label)
data$type = as.factor(data$type)
summary(data)

# Create a subset of matrices based on the months.
months = c(1,3,5,7,9,11)
matrix = list()

# Loop over the data and select the individual ID, the parameters and their values for each month.
for (i in 1:6){
  # Select the columns from the 'molten' data.
  df = subset(data, (month == months[i]), select = c("person_id", "label", "value"))
  # Arrange the subset into a dataframe/matrix format.
  df = pivot_wider(df, names_from = label, values_from = value)
  # Remove the individual ID column.
  df$person_id = NULL
  # Normalize the parameters using min-max scaling.
  preproc = preProcess(df, method = c("range"))
  norm.df = predict(preproc, df)
  #norm.df = scale(df)
  # Transform the subset as matrix object.
  mat = as.matrix(norm.df)
  # Name the rows with indivual IDs.
  rownames(mat) = paste0("IAM", 1:30)
  # Transpose the matrix to have parameters as rows and IDs as columns.
  mat = t(mat)
  # Store the matrix into a list object.
  matrix[[i]] = mat
}

# Visualization purposes.
mat = matrix[[6]]
summary(matrix[[6]])
Heatmap(matrix[[6]])

# Define clusters to the matrix[[6]] = data from month 11 using k-means/k = 3 to set up the splits of the circos plot.
set.seed(123)
km = kmeans(matrix[[6]], centers = 3)$cluster

# Create the figure.
#svg("Circosplot.svg", width = 10, height = 10)
setEPS()
postscript("Circosplot_redim.eps", width = 10, height = 10)
plot.new()
circle_size = unit(1, "snpc")
pushViewport(viewport(x = 0, y = 1, width = circle_size, height = circle_size,
                      just = c("left", "center")))
# Define the color function for the parameters values range.
col_fun1 = colorRamp2(c(0, 0.5, 1), c("#34A3DC", "white", "#F5C000")) #Light blue, white, orange.
# Initiate the layout with the matrix from month 11. Input the clusters from k-means as splits.
circos.heatmap(matrix[[6]], col = col_fun1, rownames.side = "outside", track.height = 0.1, rownames.cex = 0.7, split = km)
circos.heatmap(matrix[[5]], col = col_fun1, track.height = 0.1) # Matrix from month 9.
circos.heatmap(matrix[[4]], col = col_fun1, track.height = 0.05) # Matrix from month 7.
circos.heatmap(matrix[[3]], col = col_fun1, track.height = 0.05) # Matrix from month 5.
circos.heatmap(matrix[[2]], col = col_fun1, track.height = 0.05) # Matrix from month 3.
circos.heatmap(matrix[[1]], col = col_fun1, track.heigh = 0.05) # Matrix from month 1.
# Set up the legend.
lgd_par = Legend(title = "Range", col_fun = col_fun1) 

# Create the adjacency matrix for the connections. This matrix is defined with from_index column (start point) and to_index column (endpoint). Values are the row indices of the parameters.
# In this case, the parameters are in sequential order.
df_link = data.frame(
  from_index = c(1,2,3,4,5,6,7,8,9,10,11),
  to_index = c(12,13,14,15,16,17,18,19,20,21,22)
)

# Connect each parameter based on its location in the plot.
for(i in seq_len(nrow(df_link))) {
  group1 = km[ df_link$from_index[i] ]
  group2 = km[ df_link$to_index[i] ]
  subset1 = get.cell.meta.data("subset", sector.index = group1)
  row_order1 = get.cell.meta.data("row_order", sector.index = group1)
  x1 = which(subset1[row_order1] == df_link$from_index[i])
  subset2 = get.cell.meta.data("subset", sector.index = group2)
  row_order2 = get.cell.meta.data("row_order", sector.index = group2)
  x2 = which(subset2[row_order2] == df_link$to_index[i])
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = rand_color(1, hue = "red"))
}
circos.clear()

upViewport()
h = dev.size()[2]
lgdlist = packLegend(lgd_par, max_height = unit(0.9*h, "inch"))
# Draw the legend.
draw(lgdlist, x = unit(240, "mm"), y = unit(120, "mm"))
dev.off()