## Plots for HEK293-AAV paper


library(dplyr)
library(readxl)
library(remotes)
library(webchem)
library(stringr)
library(readr)
library(ggplot2)

library(homals) # for gifi method
library(logisticPCA) # for logistic PCA ans logistic SVD
library(pcaMethods) # for PCA method
library(readr)
library(UpSetR)

install_version("readMzXmlData", "2.8.1")
install.packages("LipidMS")

library(LipidMS)

set.seed(42)

#############################
#### Abs. BM Composition ####
#############################

####------------------------------####
#### Convert Lipids from nm to ng ####
####------------------------------####


Results_Lipidomics_cells <- read_excel("~/data/hek/lipidome/Results_Lipidomics_cells.xlsx", sheet = "Results")
refmet <- read.csv("~/data/hek/lipidome/refmet.csv")

Results_Lipidomics_cells[8, 2] <- "Cer 34:1;O2"
Results_Lipidomics_cells[9, 2] <- "Cer 36:1;O2"
Results_Lipidomics_cells[10, 2] <- "Cer 40:1;O2"
Results_Lipidomics_cells[11, 2] <- "Cer 42:1;O2"
Results_Lipidomics_cells[12, 2] <- "Cer 42:2;O2"
Results_Lipidomics_cells[126, 2] <- "SM 32:1;O2"
Results_Lipidomics_cells[127, 2] <- "SM 32:5;O2"
Results_Lipidomics_cells[128, 2] <- "SM 34:1;O2"
Results_Lipidomics_cells[129, 2] <- "SM 34:2;O2"
Results_Lipidomics_cells[130, 2] <- "SM 36:1;O2"
Results_Lipidomics_cells[131, 2] <- "SM 36:2;O2"
Results_Lipidomics_cells[132, 2] <- "SM 36:4;O2"
Results_Lipidomics_cells[133, 2] <- "SM 38:1;O2"
Results_Lipidomics_cells[134, 2] <- "SM 38:4;O2"
Results_Lipidomics_cells[135, 2] <- "SM 40:1;O2"
Results_Lipidomics_cells[136, 2] <- "SM 40:2;O2"
Results_Lipidomics_cells[137, 2] <- "SM 42:1;O2"
Results_Lipidomics_cells[138, 2] <- "SM 42:2;O2"
Results_Lipidomics_cells[139, 2] <- "SM 42:3;O2"
Results_Lipidomics_cells[140, 2] <- "SM 42:4;O2"
Results_Lipidomics_cells[141, 2] <- "SM 44:4;O2"
Results_Lipidomics_cells[142, 2] <- "SM 44:5;O2"


colnames(refmet)[colnames(refmet) == "refmet_name"] <- "Alyte"

df_merged <- merge(Results_Lipidomics_cells, refmet, by.x = "Alyte", by.y = "Alyte", all.x = TRUE)

df_merged[29, 88] <- 1654.4 # molar mass from Pubchem for CL 86:4
df_merged[67, 88] <- 689.9 # molar mass from Pubchem for PE 32:01

df_merged_1 <- df_merged[, -c(2,3,91, 90, 89, 87, 86, 85, 84)]



rownames(df_merged_1) <- df_merged_1[, 1]
df_merged_1[, 1] <- NULL
df_merged_1[df_merged_1 == "<LOQ"] <- 0

# df_merged[8, 1] <- "Cer 34:1;O2"
df_merged_2 <- apply(df_merged_1, 2, as.numeric)



df_result <- df_merged_2[, -ncol(df_merged_2)] * df_merged_2[, ncol(df_merged_2)]

df_sum <- rbind(df_result, colSums(df_result))
rownames(df_sum)[nrow(df_sum)] <- "Column_Sum"
#rownames(df_sum) <- rownames(df_merged_1)

####----------------------------####
#### Calculate BM from raw data ####
####----------------------------####

hek_data <- read_excel("data/hek/biomass/hek_data.xlsx", sheet = "hek_data")

data <- hek_data[c(1:6), ]

data <- data[, -c(82)]

lipid_sum <- df_sum[nrow(df_sum),]

new_row <- data.frame(description = "Lipid [ng/200µl]", t(lipid_sum))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert BM to pg / cell ####

row1 <- as.numeric(data[1, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- (row1 / row2) * 10^12 

# Create a new row with the calculated values and description
new_row <- data.frame(description = "BM [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert DNA to pg / cell ####

row1 <- as.numeric(data[3, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- (row1 * 200 / 10^6) / row2 * 10^9 

# Create a new row with the calculated values and description
new_row <- data.frame(description = "DNA [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert RNA to pg / cell ####

row1 <- as.numeric(data[4, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- (row1 / 10^6 * 50 / 30 * 600 ) / row2 * 10^9 

# Create a new row with the calculated values and description
new_row <- data.frame(description = "RNA [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert Glucose to pg / cell ####

row1 <- as.numeric(data[5, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- (row1 * 250 / 30 / 1000) / row2 * 10^9 

# Create a new row with the calculated values and description
new_row <- data.frame(description = "Glucose [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert Protein to pg / cell ####

row1 <- as.numeric(data[6, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- (row1 / row2) * 10^9 

# Create a new row with the calculated values and description
new_row <- data.frame(description = "Prot [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Convert Lipids to pg / cell ####

row1 <- as.numeric(data[7, -1])
row2 <- as.numeric(data[2, -1])

# Perform the division and multiplication
new_values <- row1 / row2 / 200 * 400 * 10^3

# Create a new row with the calculated values and description
new_row <- data.frame(description = "Lipids [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

# Bind the new row to the first dataframe
data <- rbind(data, new_row)

#### Calculate remaining amount in cells in pg / cell ###

row1 <- as.numeric(data[8, -1])
row2 <- as.numeric(data[9, -1])
row3 <- as.numeric(data[10, -1])
row4 <- as.numeric(data[11, -1])
row5 <- as.numeric(data[12, -1])
row6 <- as.numeric(data[13, -1])

new_values <- row1 - row2 - row3 - row4 - row5 - row6

new_row <- data.frame(description = "Others [pg/cell]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

mean_hp <- mean(as.numeric(data[8, 2:41]))

# Calculate the mean of row 8 for columns 42 to 81
mean_lp <- mean(as.numeric(data[8, 42:81]))

# Print the results
mean_hp
mean_lp


#### Calculate BM composition in percent ####

row1 <- as.numeric(data[8, -1])
row2 <- as.numeric(data[9, -1])
row3 <- as.numeric(data[10, -1])
row4 <- as.numeric(data[11, -1])
row5 <- as.numeric(data[12, -1])
row6 <- as.numeric(data[13, -1])
row7 <- as.numeric(data[14, -1])

# DNA
new_values <- row2 / row1

new_row <- data.frame(description = "DNA [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

# RNA
new_values <- row3 / row1

new_row <- data.frame(description = "RNA [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

# Glucose
new_values <- row4 / row1

new_row <- data.frame(description = "Glucose [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

# Protein
new_values <- row5 / row1

new_row <- data.frame(description = "Protein [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

# Lipids
new_values <- row6 / row1

new_row <- data.frame(description = "Lipids [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

# Others
new_values <- row7 / row1

new_row <- data.frame(description = "Others [%]", t(new_values))
colnames(new_row) <- colnames(data)  # Ensure the column names match

data <- rbind(data, new_row)

bm_abs <- data[c(9:14),]
bm_rel <- data[c(15:20),]

bm_abs <- bm_abs[, -c(1)]
bm_rel <- bm_rel[, -c(1)]

compartment <- c("DNA", "RNA", "Glucose", "Protein", "Lipid", "Others")

bm_abs$compartment <- compartment
bm_rel$compartment <- compartment

bm_abs <- bm_abs %>%
  mutate(across(1:80, as.numeric))

bm_rel <- bm_rel %>%
  mutate(across(1:80, as.numeric))

##########################################################################
#### Calc means and std of abs and rel biomass composition + plotting ####
##########################################################################

means_list <- list()
sds_list <- list()

# Calculate means and standard deviations for each condition
for (i in seq(1, ncol(bm_abs) - 1, by = 4)) {
  means <- bm_abs %>%
    select(starts_with("Sample")) %>%
    select(all_of(i:(i + 3))) %>%
    rowwise() %>%
    summarise(mean = mean(c_across(everything())), .groups = 'drop')
  
  sds <- bm_abs %>%
    select(starts_with("Sample")) %>%
    select(all_of(i:(i + 3))) %>%
    rowwise() %>%
    summarise(sd = sd(c_across(everything())), .groups = 'drop')
  
  means$compartment <- bm_abs$compartment
  sds$compartment <- bm_abs$compartment
  
  means$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(means))
  sds$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(sds))
  
  means_list[[length(means_list) + 1]] <- means
  sds_list[[length(sds_list) + 1]] <- sds
}

# Combine means and standard deviations into dataframes
means_df <- bind_rows(means_list)
sds_df <- bind_rows(sds_list)

# Combine means and standard deviations into a single dataframe
result_df <- left_join(means_df, sds_df, by = c("compartment", "condition"), suffix = c("_mean", "_sd"))


new_conditions <- c(
  "Cond1" = "HP_TR_00",
  "Cond2" = "HP_MO_00",
  "Cond3" = "HP_TR_04",
  "Cond4" = "HP_MO_04",
  "Cond5" = "HP_TR_24",
  "Cond6" = "HP_MO_24",
  "Cond7" = "HP_TR_48",
  "Cond8" = "HP_MO_48",
  "Cond9" = "HP_TR_72",
  "Cond10" = "HP_MO_72",
  "Cond11" = "LP_TR_00",
  "Cond12" = "LP_MO_00",
  "Cond13" = "LP_TR_04",
  "Cond14" = "LP_MO_04",
  "Cond15" = "LP_TR_24",
  "Cond16" = "LP_MO_24",
  "Cond17" = "LP_TR_48",
  "Cond18" = "LP_MO_48",
  "Cond19" = "LP_TR_72",
  "Cond20" = "LP_MO_72"
)

new_class <- c(
  "DNA" = "6: DNA",
  "RNA" = "5: RNA",
  "Glucose" = "4: Glucose",
  "Protein" = "3: Protein",
  "Lipid" = "2: Lipid",
  "Others" = "1: Others"
)

# Replace condition names in the dataframe
result_df <- result_df %>%
  mutate(condition = recode(condition, !!!new_conditions))

result_df_1 <- result_df # for glucose plot later

result_df <- result_df %>%
  mutate(compartment = recode(compartment, !!!new_class))

df_combined <- result_df %>%
  group_by(condition) %>%
  dplyr::mutate(V_cum = cumsum(mean))

#desired_order <- c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
#                   "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
#                   "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", 
#                   "LP_MO_24", "LP_TR_24", "LP_MO_48", "LP_TR_48", 
#                   "LP_MO_72", "LP_TR_72")

# Modify the condition column in your dataset
#df_combined$condition <- factor(df_combined$condition, levels = desired_order)

df_combined$condition <- factor(df_combined$condition, 
                                levels = c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
                                           "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
                                           "", # Add an empty level to create space between HP and LP
                                           "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", "LP_MO_24", "LP_TR_24", 
                                           "LP_MO_48", "LP_TR_48", "LP_MO_72", "LP_TR_72"))

midpoints <- c(1.5, 3.5, 5.5, 7.5, 9.5,  # Midpoints between HP bars
               12.5, 14.5, 16.5, 18.5, 20.5)  # Midpoints between LP bars

# Custom x-axis labels for timepoints
timepoints <- c("00h", "04h", "24h", "48h", "72h",  # For HP group
                "00h", "04h", "24h", "48h", "72h")

df_combined$state <- ifelse(grepl("MO", df_combined$condition), "MO", "TR")

# Plot as before
png("abs_bm.png", width = 1200, height = 700)

ggplot(df_combined, aes(x = condition, y = mean, fill = compartment, alpha = state)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(aes(ymin = V_cum - sd, ymax = V_cum + sd), width = 0.2, position = "identity") +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Weight [pg / cell]",
       fill = "Compartment") +
  scale_x_discrete(breaks = levels(df_combined$condition), labels = rep("", length(levels(df_combined$condition))), drop = FALSE) +  # Hide default labels
  # Add custom timepoint labels at the midpoints between bars
  annotate("text", x = midpoints, y = -20, label = timepoints, size = 8, hjust = 0.5, vjust = 1) +
  scale_alpha_manual(values = c("MO" = 0.7, "TR" = 1)) +  # Set transparency for "MO" bars
  # Add group labels "HP" and "LP" above the bars
  annotate("text", x = 5.5, y = 550, label = "HP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for HP
  annotate("text", x = 16.5, y = 550, label = "LP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for LP
  theme(
    axis.text.x = element_blank(),  # Hide the default axis text to avoid overlap
    axis.ticks.x = element_blank(),  # Hide the ticks on the x-axis
    axis.text.y = element_text(size = 25),  # Increase y-axis text size
    axis.title.x = element_text(size = 25),  # Increase x-axis title size
    axis.title.y = element_text(size = 25),  # Increase y-axis title size
    legend.position = "none",
    plot.margin = unit(c(1, 1, 2, 1), "cm")  # Increase bottom margin to fit custom labels
  )

dev.off()

################
#### Rel BM ####
################

means_list <- list()
sds_list <- list()

# Calculate means and standard deviations for each condition
for (i in seq(1, ncol(bm_rel) - 1, by = 4)) {
  means <- bm_rel %>%
    select(starts_with("Sample")) %>%
    select(all_of(i:(i + 3))) %>%
    rowwise() %>%
    summarise(mean = mean(c_across(everything())), .groups = 'drop')
  
  sds <- bm_rel %>%
    select(starts_with("Sample")) %>%
    select(all_of(i:(i + 3))) %>%
    rowwise() %>%
    summarise(sd = sd(c_across(everything())), .groups = 'drop')
  
  means$compartment <- bm_rel$compartment
  sds$compartment <- bm_rel$compartment
  
  means$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(means))
  sds$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(sds))
  
  means_list[[length(means_list) + 1]] <- means
  sds_list[[length(sds_list) + 1]] <- sds
}

# Combine means and standard deviations into dataframes
means_df <- bind_rows(means_list)
sds_df <- bind_rows(sds_list)

# Combine means and standard deviations into a single dataframe
result_df <- left_join(means_df, sds_df, by = c("compartment", "condition"), suffix = c("_mean", "_sd"))


new_conditions <- c(
  "Cond1" = "HP_TR_00",
  "Cond2" = "HP_MO_00",
  "Cond3" = "HP_TR_04",
  "Cond4" = "HP_MO_04",
  "Cond5" = "HP_TR_24",
  "Cond6" = "HP_MO_24",
  "Cond7" = "HP_TR_48",
  "Cond8" = "HP_MO_48",
  "Cond9" = "HP_TR_72",
  "Cond10" = "HP_MO_72",
  "Cond11" = "LP_TR_00",
  "Cond12" = "LP_MO_00",
  "Cond13" = "LP_TR_04",
  "Cond14" = "LP_MO_04",
  "Cond15" = "LP_TR_24",
  "Cond16" = "LP_MO_24",
  "Cond17" = "LP_TR_48",
  "Cond18" = "LP_MO_48",
  "Cond19" = "LP_TR_72",
  "Cond20" = "LP_MO_72"
)

new_class <- c(
  "DNA" = "6: DNA",
  "RNA" = "5: RNA",
  "Glucose" = "4: Glucose",
  "Protein" = "3: Protein",
  "Lipid" = "2: Lipid",
  "Others" = "1: Others"
)

# Replace condition names in the dataframe
result_df <- result_df %>%
  mutate(condition = recode(condition, !!!new_conditions))

result_df <- result_df %>%
  mutate(compartment = recode(compartment, !!!new_class))

df_combined <- result_df %>%
  group_by(condition) %>%
  dplyr::mutate(V_cum = cumsum(mean))

# Modify the condition column in your dataset
#df_combined$condition <- factor(df_combined$condition, levels = desired_order)

df_combined$condition <- factor(df_combined$condition, 
                                levels = c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
                                           "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
                                           "", # Add an empty level to create space between HP and LP
                                           "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", "LP_MO_24", "LP_TR_24", 
                                           "LP_MO_48", "LP_TR_48", "LP_MO_72", "LP_TR_72"))

midpoints <- c(1.5, 3.5, 5.5, 7.5, 9.5,  # Midpoints between HP bars
               12.5, 14.5, 16.5, 18.5, 20.5)  # Midpoints between LP bars

# Custom x-axis labels for timepoints
timepoints <- c("00h", "04h", "24h", "48h", "72h",  # For HP group
                "00h", "04h", "24h", "48h", "72h")

df_combined$state <- ifelse(grepl("MO", df_combined$condition), "MO", "TR")

# Plot as before
png("rel_bm.png", width = 1200, height = 700)

ggplot(df_combined, aes(x = condition, y = mean, fill = compartment, alpha = state)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(aes(ymin = V_cum - sd, ymax = V_cum + sd), width = 0.2, position = "identity") +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Rel. composition of total biomass",
       fill = "Compartment") +
  scale_x_discrete(breaks = levels(df_combined$condition), labels = rep("", length(levels(df_combined$condition))), drop = FALSE) +  # Hide default labels
  # Add custom timepoint labels at the midpoints between bars
  annotate("text", x = midpoints, y = -0.05, label = timepoints, size = 8, hjust = 0.5, vjust = 1) +
  scale_alpha_manual(values = c("MO" = 0.7, "TR" = 1)) +  # Set transparency for "MO" bars
  # Add group labels "HP" and "LP" above the bars
  annotate("text", x = 5.5, y = 1.2, label = "HP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for HP
  annotate("text", x = 16.5, y = 1.2, label = "LP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for LP
  theme(
    axis.text.x = element_blank(),  # Hide the default axis text to avoid overlap
    axis.ticks.x = element_blank(),  # Hide the ticks on the x-axis
    axis.text.y = element_text(size = 25),  # Increase y-axis text size
    axis.title.x = element_text(size = 25),  # Increase x-axis title size
    axis.title.y = element_text(size = 25),  # Increase y-axis title size
    legend.position = "none",
    plot.margin = unit(c(1, 1, 2, 1), "cm")  # Increase bottom margin to fit custom labels
  )

dev.off()

# legend

png("leg_bm.png", width = 600, height = 200)

plot.new()

# Define the colors as seen in your image
colors <- c("#FF66CC", "#CC99FF", "#33CCFF", "#53b400", "#c49a00", "#f8766d")  # Example hex values

# Define the labels for the legend (from bottom to top)
labels <- c("DNA", "RNA", "Glucose", "Protein", "Lipid", "Others")

# Add a legend with points in a single row and multiple columns
legend("center", legend = labels, col = colors, pch = 16, cex = 1.2, bty = "n", ncol = length(labels))

dev.off()

######################
#### Glucose plot ####
######################

## continue with result_df_1

glucose_df <- result_df_1[result_df_1$compartment == "Glucose", ]

glucose_df$condition <- factor(glucose_df$condition, levels = desired_order)

glucose_df$condition <- factor(glucose_df$condition, 
                                levels = c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
                                           "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
                                           "", # Add an empty level to create space between HP and LP
                                           "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", "LP_MO_24", "LP_TR_24", 
                                           "LP_MO_48", "LP_TR_48", "LP_MO_72", "LP_TR_72"))

midpoints <- c(1.5, 3.5, 5.5, 7.5, 9.5,  # Midpoints between HP bars
               12.5, 14.5, 16.5, 18.5, 20.5)  # Midpoints between LP bars

# Custom x-axis labels for timepoints
timepoints <- c("00h", "04h", "24h", "48h", "72h",  # For HP group
                "00h", "04h", "24h", "48h", "72h")

glucose_df$state <- ifelse(grepl("MO", glucose_df$condition), "MO", "TR")

glucose_df$bar_color <- factor(c(rep(c("#3498DB", "#2ECC71"), 5),  # For HP group
                                 rep(c("#E67E22", "#D6EAF8"), 5)))  # For LP group

png("glucose.png", width = 1200, height = 700)

ggplot(glucose_df, aes(x = condition, y = mean, fill = bar_color)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  labs(x = "", y = "Glucose [pg / cell]", title = "") +
  scale_x_discrete(breaks = levels(glucose_df$condition), labels = rep("", length(levels(glucose_df$condition))), drop = FALSE) +  # Hide default labels
  annotate("text", x = midpoints, y = -0.2, label = timepoints, size = 8, hjust = 0.5, vjust = 1) +
  # Add group labels "HP" and "LP" above the bars
  annotate("text", x = 5.5, y = 8.8, label = "HP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for HP
  annotate("text", x = 16.5, y = 8.8, label = "LP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for LP
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25),  # Rotate x-axis labels and increase font size
        axis.text.y = element_text(size = 25),                         # Increase y-axis font size
        axis.title.x = element_text(size = 25),                        # Increase x-axis title font size
        axis.title.y = element_text(size = 25),                        # Increase y-axis title font size
        legend.position = "none",
        plot.title = element_text(size = 25, hjust = 0.5)) +
  scale_fill_manual(values = c("#D6EAF8", "#E67E22", "#2ECC71", "#3498DB"))


dev.off()

#### legend states

png("leg_tr.png", width = 700, height = 200)

plot.new()

# Define the colors as seen in your image
colors <- c("#E67E22", "#3498DB")  # Example hex values

# Define the labels for the legend (from bottom to top)
labels <- c("HP, TR", "LP, TR")

# Add a legend with points in a single row and multiple columns
legend("center", legend = labels, col = colors, pch = 16, cex = 1.2, bty = "n", ncol = length(labels))

dev.off()

##################################################################
#### Wilcoxon Rank test of biomass composition between states #### --> only works for glucose !
##################################################################

bm_abs_hp <- as.data.frame(bm_abs[, c(1:40)])
bm_abs_lp <- as.data.frame(bm_abs[, c(41:80)])

bm_abs_hp[] <- lapply(bm_abs_hp, function(x) as.numeric(as.character(x)))
bm_abs_lp[] <- lapply(bm_abs_lp, function(x) as.numeric(as.character(x)))

# Function to perform Wilcoxon rank-sum test for each row
perform_wilcoxon <- function(row1, row2) {
  p_values <- c()
  
  for (i in seq(1, 40, by = 4)) {
    group1 <- as.numeric(row1[i:(i+3)])
    group2 <- as.numeric(row2[i:(i+3)])
    test_result <- wilcox.test(group1, group2)
    p_values <- c(p_values, test_result$p.value)
  }
  
  return(p_values)
}

# Initialize a list to store the p-values for each row
results <- list()

# Loop through each row
for (row in 1:5) {
  row1 <- bm_abs_hp[row, ]
  row2 <- bm_abs_lp[row, ] # Assuming you have a second dataframe for comparison, replace with the actual second dataframe
  results[[row]] <- perform_wilcoxon(row1, row2)
}

# Convert the list to a dataframe for better visualization
results_df <- do.call(rbind, results)
colnames(results_df) <- paste0("Group_", 1:10)

# Apply multiple testing correction (Benjamini-Hochberg method)
adjusted_results_df <- apply(results_df, 1, function(pvals) {
  p.adjust(pvals, method = "BH")
})

# Convert adjusted results to a dataframe
adjusted_results_df <- as.data.frame(t(adjusted_results_df))
colnames(adjusted_results_df) <- paste0("Group_", 1:10)


################################
#### Abs. Lipid Composition ####
################################

lipids <- read_excel("data/hek/lipidome/Results_Lipidomics_cells.xlsx", sheet = "Lipid Classes")

lipids <- lipids[-c(16),]

lip_class <- lipids$Class

bm <- as.matrix(hek_data[1, ])
bm <- bm[, -c(82)]

# Combine dataframes ignoring column names
df1 <- rbind(as.matrix(lipids), bm)

# Convert back to dataframe (if necessary)
df1 <- as.data.frame(df1)

df1 <- df1[, -c(1)]

df1 <- data.frame(lapply(df1, as.numeric))

#df1 <- df1[, -c(1,2,3,4,5,6,7,8,41,42,43,44,45,46,47,48)]

#lipid_sum = as.list(df1[c(16),])

df1n <- df1 %>%
  mutate(across(everything(), ~ . / last(.)))

df1n <- df1n[-c(16),]

df1n <- df1n / 1000 # to make it in mmole and not nanomole

df1n$class <- lip_class

selected_classes <- c("ST", "TG", "SM", "PI", "PE", "PC")

means_list <- list()
sds_list <- list()

# Calculate means and standard deviations for each condition
for (i in seq(1, ncol(df1n) - 1, by = 4)) {
  means <- df1n %>%
    select(starts_with("Sample")) %>%
    select(i:(i + 3)) %>%
    rowwise() %>%
    summarise(mean = mean(c_across(everything())), .groups = 'drop')
  
  sds <- df1n %>%
    select(starts_with("Sample")) %>%
    select(i:(i + 3)) %>%
    rowwise() %>%
    summarise(sd = sd(c_across(everything())), .groups = 'drop')
  
  means$class <- df1n$class
  sds$class <- df1n$class
  
  means$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(means))
  sds$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(sds))
  
  means_list[[length(means_list) + 1]] <- means
  sds_list[[length(sds_list) + 1]] <- sds
}

# Combine means and standard deviations into dataframes
means_df <- bind_rows(means_list)
sds_df <- bind_rows(sds_list)

# Combine means and standard deviations into a single dataframe
result_df <- left_join(means_df, sds_df, by = c("class", "condition"), suffix = c("_mean", "_sd"))

filtered_df <- result_df %>%
  filter(class %in% selected_classes)

cumulative_df <- result_df %>%
  filter(!class %in% selected_classes) %>%
  group_by(condition) %>%
  summarise(across(starts_with("mean"), sum), across(starts_with("sd"), sum)) %>%
  mutate(class = "Other")

final_df <- bind_rows(filtered_df, cumulative_df)

new_conditions <- c(
  "Cond1" = "HP_TR_00",
  "Cond2" = "HP_MO_00",
  "Cond3" = "HP_TR_04",
  "Cond4" = "HP_MO_04",
  "Cond5" = "HP_TR_24",
  "Cond6" = "HP_MO_24",
  "Cond7" = "HP_TR_48",
  "Cond8" = "HP_MO_48",
  "Cond9" = "HP_TR_72",
  "Cond10" = "HP_MO_72",
  "Cond11" = "LP_TR_00",
  "Cond12" = "LP_MO_00",
  "Cond13" = "LP_TR_04",
  "Cond14" = "LP_MO_04",
  "Cond15" = "LP_TR_24",
  "Cond16" = "LP_MO_24",
  "Cond17" = "LP_TR_48",
  "Cond18" = "LP_MO_48",
  "Cond19" = "LP_TR_72",
  "Cond20" = "LP_MO_72"
)

new_class <- c(
  "PC" = "7: PC",
  "PE" = "6: PE",
  "PI" = "5: PI",
  "SM" = "4: SM",
  "TG" = "3: TG",
  "ST" = "2: ST",
  "Other" = "1: Other"
)

# Replace condition names in the dataframe
final_df <- final_df %>%
  mutate(condition = recode(condition, !!!new_conditions))

final_df <- final_df %>%
  mutate(class = recode(class, !!!new_class))

df_combined <- final_df %>%
  group_by(condition) %>%
  dplyr::mutate(V_cum = cumsum(mean))

#df_combined$condition <- factor(df_combined$condition, levels = desired_order)

df_combined$condition <- factor(df_combined$condition, 
                                levels = c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
                                           "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
                                           "", # Add an empty level to create space between HP and LP
                                           "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", "LP_MO_24", "LP_TR_24", 
                                           "LP_MO_48", "LP_TR_48", "LP_MO_72", "LP_TR_72"))

midpoints <- c(1.5, 3.5, 5.5, 7.5, 9.5,  # Midpoints between HP bars
               12.5, 14.5, 16.5, 18.5, 20.5)  # Midpoints between LP bars

# Custom x-axis labels for timepoints
timepoints <- c("00h", "04h", "24h", "48h", "72h",  # For HP group
                "00h", "04h", "24h", "48h", "72h")

df_combined$state <- ifelse(grepl("MO", df_combined$condition), "MO", "TR")

# Plot as before
png("abs_lip.png", width = 1200, height = 700)

ggplot(df_combined, aes(x = condition, y = mean, fill = class, alpha = state)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(aes(ymin = V_cum - sd, ymax = V_cum + sd), width = 0.2, position = "identity") +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Abs. composition of lipids [nmole]",
       fill = "class") +
  scale_x_discrete(breaks = levels(df_combined$condition), labels = rep("", length(levels(df_combined$condition))), drop = FALSE) +  # Hide default labels
  # Add custom timepoint labels at the midpoints between bars
  annotate("text", x = midpoints, y = -5, label = timepoints, size = 8, hjust = 0.5, vjust = 1) +
  scale_alpha_manual(values = c("MO" = 0.7, "TR" = 1)) +  # Set transparency for "MO" bars
  # Add group labels "HP" and "LP" above the bars
  annotate("text", x = 5.5, y = 160, label = "HP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for HP
  annotate("text", x = 16.5, y = 160, label = "LP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for LP
  theme(
    axis.text.x = element_blank(),  # Hide the default axis text to avoid overlap
    axis.ticks.x = element_blank(),  # Hide the ticks on the x-axis
    axis.text.y = element_text(size = 25),  # Increase y-axis text size
    axis.title.x = element_text(size = 25),  # Increase x-axis title size
    axis.title.y = element_text(size = 25),  # Increase y-axis title size
    legend.position = "none",
    plot.margin = unit(c(1, 1, 2, 1), "cm")  # Increase bottom margin to fit custom labels
  )

dev.off()

################################
#### Rel. Lipid Composition ####
################################

lipids <- read_excel("data/hek/lipidome/Results_Lipidomics_cells.xlsx", sheet = "Lipid Classes")

lipids <- lipids[-c(16),]

lip_class <- lipids$Class

bm <- as.matrix(hek_data[1, ])
bm <- bm[, -c(82)]

# Combine dataframes ignoring column names
df1 <- rbind(as.matrix(lipids), bm)

# Convert back to dataframe (if necessary)
df1 <- as.data.frame(df1)

df1 <- df1[, -c(1)]

df1 <- data.frame(lapply(df1, as.numeric))

#df1 <- df1[, -c(1,2,3,4,5,6,7,8,41,42,43,44,45,46,47,48)]

#lipid_sum = as.list(df1[c(16),])

df1n <- df1 %>%
  mutate(across(everything(), ~ . / last(.)))

df1n <- df1n[-c(16),]

#df1n <- data.frame(lapply(df1n, as.numeric))

df1n <- df1n %>%
  mutate(across(everything(), ~ . / sum(.)))


df1n$class <- lip_class

selected_classes <- c("ST", "TG", "SM", "PI", "PE", "PC")

means_list <- list()
sds_list <- list()

# Calculate means and standard deviations for each condition
for (i in seq(1, ncol(df1n) - 1, by = 4)) {
  means <- df1n %>%
    select(starts_with("Sample")) %>%
    select(i:(i + 3)) %>%
    rowwise() %>%
    summarise(mean = mean(c_across(everything())), .groups = 'drop')
  
  sds <- df1n %>%
    select(starts_with("Sample")) %>%
    select(i:(i + 3)) %>%
    rowwise() %>%
    summarise(sd = sd(c_across(everything())), .groups = 'drop')
  
  means$class <- df1n$class
  sds$class <- df1n$class
  
  means$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(means))
  sds$condition <- rep(paste0("Cond", (i + 3) / 4), nrow(sds))
  
  means_list[[length(means_list) + 1]] <- means
  sds_list[[length(sds_list) + 1]] <- sds
}

# Combine means and standard deviations into dataframes
means_df <- bind_rows(means_list)
sds_df <- bind_rows(sds_list)

# Combine means and standard deviations into a single dataframe
result_df <- left_join(means_df, sds_df, by = c("class", "condition"), suffix = c("_mean", "_sd"))

filtered_df <- result_df %>%
  filter(class %in% selected_classes)

cumulative_df <- result_df %>%
  filter(!class %in% selected_classes) %>%
  group_by(condition) %>%
  summarise(across(starts_with("mean"), sum), across(starts_with("sd"), sum)) %>%
  mutate(class = "Other")

final_df <- bind_rows(filtered_df, cumulative_df)

new_conditions <- c(
  "Cond1" = "HP_TR_00",
  "Cond2" = "HP_MO_00",
  "Cond3" = "HP_TR_04",
  "Cond4" = "HP_MO_04",
  "Cond5" = "HP_TR_24",
  "Cond6" = "HP_MO_24",
  "Cond7" = "HP_TR_48",
  "Cond8" = "HP_MO_48",
  "Cond9" = "HP_TR_72",
  "Cond10" = "HP_MO_72",
  "Cond11" = "LP_TR_00",
  "Cond12" = "LP_MO_00",
  "Cond13" = "LP_TR_04",
  "Cond14" = "LP_MO_04",
  "Cond15" = "LP_TR_24",
  "Cond16" = "LP_MO_24",
  "Cond17" = "LP_TR_48",
  "Cond18" = "LP_MO_48",
  "Cond19" = "LP_TR_72",
  "Cond20" = "LP_MO_72"
)

new_class <- c(
  "PC" = "7: PC",
  "PE" = "6: PE",
  "PI" = "5: PI",
  "SM" = "4: SM",
  "TG" = "3: TG",
  "ST" = "2: ST",
  "Other" = "1: Other"
)

# Replace condition names in the dataframe
final_df <- final_df %>%
  mutate(condition = recode(condition, !!!new_conditions))

final_df <- final_df %>%
  mutate(class = recode(class, !!!new_class))

df_combined <- final_df %>%
  group_by(condition) %>%
  dplyr::mutate(V_cum = cumsum(mean))

#df_combined$condition <- factor(df_combined$condition, levels = desired_order)

df_combined$condition <- factor(df_combined$condition, 
                                levels = c("HP_MO_00", "HP_TR_00", "HP_MO_04", "HP_TR_04", "HP_MO_24", "HP_TR_24", 
                                           "HP_MO_48", "HP_TR_48", "HP_MO_72", "HP_TR_72", 
                                           "", # Add an empty level to create space between HP and LP
                                           "LP_MO_00", "LP_TR_00", "LP_MO_04", "LP_TR_04", "LP_MO_24", "LP_TR_24", 
                                           "LP_MO_48", "LP_TR_48", "LP_MO_72", "LP_TR_72"))

midpoints <- c(1.5, 3.5, 5.5, 7.5, 9.5,  # Midpoints between HP bars
               12.5, 14.5, 16.5, 18.5, 20.5)  # Midpoints between LP bars

# Custom x-axis labels for timepoints
timepoints <- c("00h", "04h", "24h", "48h", "72h",  # For HP group
                "00h", "04h", "24h", "48h", "72h")

df_combined$state <- ifelse(grepl("MO", df_combined$condition), "MO", "TR")

# Plot as before
png("rel_lip.png", width = 1200, height = 700)

ggplot(df_combined, aes(x = condition, y = mean, fill = class, alpha = state)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(aes(ymin = V_cum - sd, ymax = V_cum + sd), width = 0.2, position = "identity") +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Rel. composition of lipids [nmole/nmole]",
       fill = "class") +
  scale_x_discrete(breaks = levels(df_combined$condition), labels = rep("", length(levels(df_combined$condition))), drop = FALSE) +  # Hide default labels
  # Add custom timepoint labels at the midpoints between bars
  annotate("text", x = midpoints, y = -0.05, label = timepoints, size = 8, hjust = 0.5, vjust = 1) +
  scale_alpha_manual(values = c("MO" = 0.7, "TR" = 1)) +  # Set transparency for "MO" bars
  # Add group labels "HP" and "LP" above the bars
  annotate("text", x = 5.5, y = 1.1, label = "HP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for HP
  annotate("text", x = 16.5, y = 1.1, label = "LP", size = 8, hjust = 0.5, vjust = 0) +  # Adjust x position and size for LP
  theme(
    axis.text.x = element_blank(),  # Hide the default axis text to avoid overlap
    axis.ticks.x = element_blank(),  # Hide the ticks on the x-axis
    axis.text.y = element_text(size = 25),  # Increase y-axis text size
    axis.title.x = element_text(size = 25),  # Increase x-axis title size
    axis.title.y = element_text(size = 25),  # Increase y-axis title size
    legend.position = "none",
    plot.margin = unit(c(1, 1, 2, 1), "cm")  # Increase bottom margin to fit custom labels
  )

dev.off()


png("leg_lip.png", width = 600, height = 110)

plot.new()

# Define the colors as seen in your image
colors <- c("#FF66CC", "#CC99FF", "#33CCFF", "#00c094" ,"#53b400", "#c49a00", "#f8766d")  # Example hex values

# Define the labels for the legend (from bottom to top)
labels <- c("PC", "PE", "PI", "SM", "TG", "ST", "Others")

# Add a legend with points in a single row and multiple columns
legend("center", legend = labels, col = colors, pch = 16, cex = 1.2, bty = "n", ncol = length(labels))

dev.off()

########################
#### BM vs VOL plot ####
########################

df <- as.data.frame(t(bm_vs_vol))

df$Group <- cut(as.numeric(gsub("Sample_", "", rownames(df))), 
                breaks=c(0, 20, 40, 60, 80, 83), 
                labels=c("HP_MO", "HP_TR", "LP_MO", "LP_TR", "CHO"))

df$GroupType <- ifelse(grepl("HP", df$Group), "HP",
                       ifelse(grepl("LP", df$Group), "LP", NA))

df$GroupType <- as.factor(df$GroupType)

biomass_test <- t.test(bm ~ GroupType, data = df)
print("Biomass T-Test Results:")
print(biomass_test)

volume_test <- t.test(vol ~ GroupType, data = df)
print("Volume T-Test Results:")
print(volume_test)

df$density <- df$bm / df$vol

density_test <- t.test(density ~ GroupType, data = df)
print("Volume T-Test Results:")
print(volume_test)

raw_p_values <- c(
  biomass_test$p.value,
  volume_test$p.value,
  density_test$p.value
)

# Apply multiple testing correction
# Bonferroni correction
adjusted_p_values_bonferroni <- p.adjust(raw_p_values, method = "bonferroni")

# False Discovery Rate correction
adjusted_p_values_fdr <- p.adjust(raw_p_values, method = "fdr")

# Create a summary table of the results
results <- data.frame(
  Test = c("Biomass", "Volume", "Density"),
  Raw_P_Value = raw_p_values,
  Bonferroni_Adjusted_P = adjusted_p_values_bonferroni
#  FDR_Adjusted_P = adjusted_p_values_fdr
)

# Print the summary table
print("Multiple Testing Correction Results:")
print(results)

box_colors <- c("#D6EAF8", "#E67E22", "#2ECC71","#3498DB", "red")

lighter_colors <- adjustcolor(box_colors, alpha.f = 0.5)

png("bm_vs_vol.png", width = 900, height = 700)

par(mar = c(5, 6, 4, 2) + 0.1)

boxplot(density ~ Group, data = df, main = "", 
        xlab = "Condition", ylab = "Density [gDW / µm^3]", 
        col = lighter_colors,
        cex.lab = 2,
        cex.axis = 2,
        cex.names = 2,
        cex.main = 2,
        lwd = 3)

stripchart(density ~ Group, data = df, 
           vertical = TRUE, method = "jitter", 
           pch = 1, col = "black", cex = 2, add = TRUE)

dev.off()

#########################
#### Lipid molecules ####
#########################

lipids <- read_excel("data/hek/lipidome/Results_Lipidomics_cells.xlsx", sheet = "Results")

lip_class <- lipids$Class

df <- lipids[, -c(1:3)]

df1 <- df %>%
  mutate(across(everything(), ~ ifelse(. == "<LOQ", 0, .)))

bm <- as.matrix(hek_data[1, ])
bm <- bm[, -c(1,82)]

# Combine dataframes ignoring column names
df1 <- rbind(as.matrix(df1), bm)

# Convert back to dataframe (if necessary)
df1 <- as.data.frame(df1)

df1 <- data.frame(lapply(df1, as.numeric))

df1 <- df1[, -c(1,2,3,4,5,6,7,8,41,42,43,44,45,46,47,48)]

#lipid_sum = as.list(df1[c(16),])

df1n <- df1 %>%
  mutate(across(everything(), ~ . / last(.)))

df1n <- df1n[-c(168),]


df1_t <- t(df1n)

data.pca <- prcomp(df1_t, center = TRUE, scale = TRUE)

summary(data.pca)

loadings <- as.data.frame(data.pca$rotation[, c(1,2)])
scores <- as.data.frame(data.pca$x[, c(1,2)])

sample_ids = c("HP_TR", "HP_TR", "HP_TR", "HP_TR", 
               "HP_MO", "HP_MO", "HP_MO", "HP_MO", 
               "HP_TR", "HP_TR", "HP_TR", "HP_TR", 
               "HP_MO", "HP_MO", "HP_MO", "HP_MO",
               "HP_TR", "HP_TR", "HP_TR", "HP_TR", 
               "HP_MO", "HP_MO", "HP_MO", "HP_MO",
               "HP_TR", "HP_TR", "HP_TR", "HP_TR", 
               "HP_MO", "HP_MO", "HP_MO", "HP_MO",
               "LP_TR", "LP_TR", "LP_TR", "LP_TR", 
               "LP_MO", "LP_MO", "LP_MO", "LP_MO", 
               "LP_TR", "LP_TR", "LP_TR", "LP_TR", 
               "LP_MO", "LP_MO", "LP_MO", "LP_MO",
               "LP_TR", "LP_TR", "LP_TR", "LP_TR", 
               "LP_MO", "LP_MO", "LP_MO", "LP_MO",
               "LP_TR", "LP_TR", "LP_TR", "LP_TR", 
               "LP_MO", "LP_MO", "LP_MO", "LP_MO"
)

df <- cbind(scores, sample_ids)

loadings$class <- lip_class

loadings$length <- sqrt(loadings$PC1^2 + loadings$PC2^2)

loading_results <- loadings %>%
  group_by(class) %>%
  summarize(
    loading_1 = mean(PC1),
    loading_2 = mean(PC2)
  )

loading_results$length <- sqrt(loading_results$loading_1^2 + loading_results$loading_2^2)

# Sort the dataframe by 'length' in descending order
sorted_result <- loading_results %>% arrange(desc(length))
# pathway_importance_lpca <- sorted_result[, c(1,4)]

# Take the top 10 longest vectors
top_longest_vectors <- head(sorted_result, 5)

class.labels <- c()

class.labels[c(1,2,3,4,9,10,11,12,17,18,19,20,25,26,27,28)] = "#E67E22" # 
class.labels[c(5,6,7,8,13,14,15,16,21,22,23,24,29,30,31,32)] = "#D6EAF8" # 
class.labels[c(33,34,35,36,41,42,43,44,49,50,51,52,57,58,59,60)] = "#3498DB" #  
class.labels[c(37,38,39,40,45,46,47,48,53,54,55,56,61,62,63,64)] = "#2ECC71" # 

class.shapes <- c()

class.shapes[c(1,2,3,4,5,6,7,8,33,34,35,36,37,38,39,40)] = 21 # 
class.shapes[c(9,10,11,12,13,14,15,16,41,42,43,44,45,46,47,48)] = 22 # 
class.shapes[c(17,18,19,20,21,22,23,24,49,50,51,52,53,54,55,56)] = 24 #  
class.shapes[c(25,26,27,28,29,30,31,32,57,58,59,60,61,62,63,64)] = 4 #

# Reset margins to default (optional)
par(mar = c(5, 4, 4, 2) + 0.1)

png("pca_lipidome.png", width = 1200, height = 700)

par(mar = c(8, 9, 5, 2) + 0.1)

# Adjust the mgp parameter to move the axis labels and axis titles
par(mgp = c(6, 2.5, 0))

plot(scores,  # x and y data
     #     pch=21,           # point shape
     pch=class.shapes,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=4,          # point size
     main="",     # title of plot
     xlab = "PC1 (27%)", # 
     ylab = "PC2 (21%)", # 
     cex.axis = 3,
     cex.lab = 3
)

arrows(x0 = rep(0, nrow(top_longest_vectors)), 
       y0 = rep(0, nrow(top_longest_vectors)), 
       x1 = top_longest_vectors$loading_1*100, 
       y1 = top_longest_vectors$loading_2*100, 
       col = "black", # You can change the color
       angle = 25, # Angle of the arrow head
       length = 0.3,
       lwd = 3) # Length of the arrow head

for(i in 1:nrow(top_longest_vectors)) {
  if (top_longest_vectors$loading_1[i] < 0) {
    text(x = top_longest_vectors$loading_1[i]*100, 
         y = top_longest_vectors$loading_2[i]*100, 
         labels = top_longest_vectors$class[i],
         pos = 2,  # Places the text above the point
         col = "black",  # Text color
         cex = 3  # Text size
    )
  }
  if (top_longest_vectors$loading_1[i] > 0) {
    text(x = top_longest_vectors$loading_1[i]*100, 
         y = top_longest_vectors$loading_2[i]*100, 
         labels = top_longest_vectors$class[i],
         pos = 4,  # Places the text above the point
         col = "black",  # Text color
         cex = 3  # Text size
    )
  }
}

#colors <- unique(class.labels)
#legend("bottomright", legend=c("HP, TR", "HP, MO", "LP, TR", "LP, MO"), fill=colors, cex=3)

colors <- c("#E67E22", "#D6EAF8", "#3498DB", "#2ECC71")  # Corresponding colors
shapes <- c(21, 22, 24, 4)  # Corresponding shapes

# Add combined color and shape legend
legend("topright",
       legend = c("HP, TR", "HP, MO", "LP, TR", "LP, MO"),  # Labels for colors
       fill = colors,  # Color legend
       title = "Conditions",  # Title for color legend
       cex = 2
)

# Add separate shape legend
legend("bottomright", 
       legend = c("4 HPT", "24 HPT", "48 HPT", "72 HPT"),  # Labels for shapes
       pch = shapes,  # Shape legend
       pt.bg = "black",  # Point fill color (can be adjusted if necessary)
       cex = 2,
       title = "Timepoints"  # Title for shape legend
)

dev.off()


###########################
#### Transcriptome PCA ####
###########################

df_1 <- read_csv("data/hek/transcriptome/tpm_P.csv")
df_2 <- read_csv("data/hek/transcriptome/tpm_I.csv")

df_1 <- df_1[,-1]
df_2 <- df_2[,-1]

df_3 <- df_1 %>% left_join(df_2, by = "gene")
df_3 <- df_3[, -1]
df_3[is.na(df_3)] <- 0

if(any(is.na(df_3))){
  print("Dataframe contains NA values.")
} else {
  print("Dataframe does not contain NA values.")
}

less_than_ <- df_3 < 0.5
numerical_data <- df_3[!(rowSums(less_than_) > ncol(df_3) / 2),]

df_t <- t(numerical_data)

data.pca <- prcomp(df_t, center = TRUE, scale = TRUE)

summary(data.pca)

loadings <- as.data.frame(data.pca$rotation[, c(1,2)])
scores <- as.data.frame(data.pca$x[, c(1,2)])

scores <- scores * -1 

col_names <- colnames(numerical_data)
class.labels = vector()

for (name in col_names) {
  
  # Group P, state T
  if (substr(name, 1, 1) == 'P' & substr(name, 4, 4) == 'T') {
    class.labels <- c(class.labels, "#E67E22") # orange
  }
  
  # Group P, state M
  else if (substr(name, 1, 1) == 'P' & substr(name, 4, 4) == 'M') {
    class.labels <- c(class.labels, "#D6EAF8") # turquois
  }
  
  # Group I, state T
  else if (substr(name, 1, 1) == 'I' & substr(name, 4, 4) == 'T') {
    class.labels <- c(class.labels, "#3498DB") # 
  }
  
  # Group I, state M
  else if (substr(name, 1, 1) == 'I' & substr(name, 4, 4) == 'M') {
    class.labels <- c(class.labels, "#2ECC71")
  }
}

class.shapes = c()

for (name in col_names) {
  
  # Group I, state M (handling 2-digit substrings, i.e., '13')
  if (substr(name, 6, 7) == '13') {
    class.shapes <- c(class.shapes, 24) 
  }
  
  # Another condition for 2-digit substrings (i.e., '17')
  else if (substr(name, 6, 7) == '17') {
    class.shapes <- c(class.shapes, 4) 
  }
  
  # Group P, state T
  else if (substr(name, 6, 6) == '1') {
    class.shapes <- c(class.shapes, 9) # orange
  }
  
  # Group P, state M
  else if (substr(name, 6, 6) == '5') {
    class.shapes <- c(class.shapes, 21) # turquoise
  }
  
  # Group I, state T
  else if (substr(name, 6, 6) == '9' ) {
    class.shapes <- c(class.shapes, 22) # 
  }
}

#class.shapes[c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,75)] = 9 # 
#class.shapes[c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,76)] = 21 # 
#class.shapes[c(3,8,13,18,23,28,33,38,43,48,53,58,63,68,73,77)] = 22 #  
#class.shapes[c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,78)] = 24 #
#class.shapes[c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,79)] = 4 #

png("pca_hp_lp.png", width = 1200, height = 700)

par(mar = c(8, 9, 5, 2) + 0.1)

# Adjust the mgp parameter to move the axis labels and axis titles
par(mgp = c(6, 2.5, 0))

plot(scores,  # x and y data
     pch=class.shapes,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=4,          # point size
     #     main="PCA of HEK293 strains during rAAV production",     # title of plot
     xlab = "PC1 (27%)",
     ylab = "PC2 (11%)",
     cex.lab = 3,
     cex.axis = 3
)

colors <- c("#E67E22", "#D6EAF8", "#3498DB", "#2ECC71")  # Corresponding colors
shapes <- c(9, 21, 22, 24, 4)  # Corresponding shapes

# Add combined color and shape legend
legend("top",
       legend = c("HP, TR", "HP, MO", "LP, TR", "LP, MO"),  # Labels for colors
       fill = colors,  # Color legend
       title = "Conditions",  # Title for color legend
       cex = 2
)

# Add separate shape legend
legend("bottom", 
       legend = c("0 HPT", "4 HPT", "24 HPT", "48 HPT", "72 HPT"),  # Labels for shapes
       pch = shapes,  # Shape legend
       pt.bg = "black",  # Point fill color (can be adjusted if necessary)
       cex = 2,
       title = "Timepoints"  # Title for shape legend
)

dev.off()

####################################
#### statistical test for genes ####
####################################

df_1 <- read_csv("data/hek/transcriptome/tpm_P.csv")
df_2 <- read_csv("data/hek/transcriptome/tpm_I.csv")

df_1 <- df_1[,-1]
df_2 <- df_2[,-1]

df_3 <- df_1 %>% left_join(df_2, by = "gene")

genes2test <- c("ENSG00000100644", "ENSG00000103257", "ENSG00000079215", "ENSG00000184432",
                "ENSG00000007350", "ENSG00000177156", "ENSG00000137204", "ENSG00000159423", 
                "ENSG00000138413", "ENSG00000160211")

# HIF1alpha: ENSG00000100644 
# SLC7A5: ENSG00000103257 
# SLC1A3: ENSG00000079215 
# COPB2: ENSG00000184432 
# TKTL1: ENSG00000007350 
# TALDO1: ENSG00000177156 
#-------------------------
# SLC22A7: ENSG00000137204 
# ALDH4A1: ENSG00000159423 
# IDH1: ENSG00000138413
# G6PD: ENSG00000160211 

# Initialize a data frame to store the results
results_df <- data.frame(
  gene = character(),
  t_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each gene to perform t-tests
for (gene_id in genes2test) {
  # Filter the data for the current gene
  gene_data <- df_3 %>%
    filter(gene == gene_id)
  
  if (nrow(gene_data) == 0) {
    # If the gene is not found, print a message and continue to the next gene
    message(paste("Gene", gene_id, "not found in the data."))
    next
  }
  
  # Extract expression values, excluding the 'gene' column
  expression_values <- as.numeric(gene_data[, -1])
  
  # Define indices for group 1 and group 2
  group1_indices <- 1:40
  group2_indices <- 41:length(expression_values)
  
  # Check if the expression values have enough samples
  if (length(expression_values) < max(group2_indices)) {
    message(paste("Gene", gene_id, "does not have enough samples for both groups."))
    next
  }
  
  # Split the expression values into two groups
  group1 <- expression_values[group1_indices]
  group2 <- expression_values[group2_indices]
  
  # Perform the t-test between the two groups
  t_test_result <- t.test(group1, group2)
  
  # Store the results in the results data frame
  results_df <- rbind(results_df, data.frame(
    gene = gene_id,
    t_statistic = t_test_result$statistic,
    p_value = t_test_result$p.value,
    stringsAsFactors = FALSE
  ))
}

# Apply Bonferroni correction to adjust p-values
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "bonferroni")

# Print the results
print(results_df)


##############################################
#### Analysis after addition of inhibitor ####
##############################################

data <- read_excel("data/paper4_aav/HEK_AAV/aav_inh_results_final.xlsx")

df <- data[,-c(1)]

timepoints <- as.numeric(df[1,])
values <- as.numeric(df[4,])

original_group <- rep(1:6, times = 4)
combined_group <- original_group
combined_group[original_group == 2] <- 1
combined_group[original_group == 4] <- 3
combined_group[original_group == 6] <- 5
#combined_group[original_group == 5] <- 3
#combined_group[original_group == 6] <- 3

plot_data <- data.frame(timepoints, values, group = factor(combined_group))

# Rename the levels of the group variable
levels(plot_data$group) <- c("Supplementation", "Inhibitor", "Reference")

custom_colors <- c("Supplementation" = "midnightblue", "Inhibitor" = "royalblue", "Reference" = "darkred")

mean_data <- plot_data %>%
  group_by(timepoints, group) %>%
  summarize(mean_values = mean(values))

# Save the plot as a PNG file
png("growth_validation_experiment.png", width = 1200, height = 700)

# Create the dotplot
ggplot(plot_data, aes(x = timepoints, y = values, color = group)) +
  geom_point(size = 15) +  # Plot the replicates without transparency
  geom_line(data = mean_data, aes(x = timepoints, y = mean_values, group = group), size = 2) +  # Add lines for mean values
  labs(x = "Time post transfection [h]", y = "Cell count per ml", title = "") +
  scale_color_manual(values = custom_colors, name = NULL) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.text = element_text(size = 25)   # Customize legend text size
#    legend.position.inside = c(1, 1),               # Position the legend inside the plot at top-left
#    legend.justification = c(1, 1)
)

dev.off()

####################################################
#### AAV production after addition of inhibitor ####
####################################################

df <- data[,-c(1,2,3,4,5,6,7)]

df[] <- lapply(df, function(x) as.numeric(as.character(x)))

timepoints <- as.numeric(df[1,])
values <- as.numeric(df[3,])  # Use the new row for values

# Create a group variable and combine groups 1 & 2, and 7 & 8
original_group <- rep(1:6, times = 3)
combined_group <- original_group
combined_group[original_group == 2] <- 1
combined_group[original_group == 4] <- 3
combined_group[original_group == 6] <- 5
#combined_group[original_group == 5] <- 3
#combined_group[original_group == 6] <- 3

# Combine into a new dataframe
plot_data <- data.frame(timepoints, values, group = factor(combined_group))

# Rename the levels of the group variable to specific names
levels(plot_data$group) <- c("Supplementation", "Inhibitor", "Reference")

custom_colors <- c("Supplementation" = "midnightblue", "Inhibitor" = "royalblue", "Reference" = "darkred")

# Save the plot as a PNG file
png("production_validation_experiment.png", width = 1200, height = 700)

# Create the dotplot
ggplot(plot_data, aes(x = timepoints, y = values, color = group)) +
  geom_point(size = 15) +  # Increased dot size
  labs(x = "Time post transfection [h]", y = "Viral Capsids per ml", title = "") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
#    plot.title = element_text(size = 20),
    legend.position = "none")

dev.off()

############################################
#### Specific AAV production validation ####
############################################

df <- data[,-c(1,2,3,4,5,6,7)]

df[] <- lapply(df, function(x) as.numeric(as.character(x)))

# df[3] for VG and df[2] for VC
df[5,] <- df[2,] / df[4,] / 10^6

# Reshape the data into a long format
timepoints <- as.numeric(df[1,])
values <- as.numeric(df[5,])  # Use the new row for values

# Create a group variable and combine groups 1 & 2, and 7 & 8
original_group <- rep(1:6, times = 3)
combined_group <- original_group
combined_group[original_group == 2] <- 1
combined_group[original_group == 4] <- 3
combined_group[original_group == 6] <- 5
#combined_group[original_group == 5] <- 3
#combined_group[original_group == 6] <- 3

# Combine into a new dataframe
plot_data <- data.frame(timepoints, values, group = factor(combined_group))

# Rename the levels of the group variable to specific names
levels(plot_data$group) <- c("Supplementation", "Inhibitor", "Reference")

custom_colors <- c("Supplementation" = "midnightblue", "Inhibitor" = "royalblue", "Reference" = "darkred")

# Save the plot as a PNG file
png("spec_production_VC_validation_experiment.png", width = 1200, height = 700)

# Create the dotplot
ggplot(plot_data, aes(x = timepoints, y = values, color = group)) +
  geom_point(size = 15) +  # Increased dot size
  labs(x = "Time post transfection [h]", y = "Viral capsids per cell", title = "") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    legend.position = "none")

dev.off()
#### Statistical test ####

values <- as.numeric(df[4, 13:18])

# Define the groups
group1 <- values[c(5,6)]  # Columns 17, 18, 23, 24
group2 <- values[c(3,4)]  # Columns 19, 20, 21, 22

# Perform the Wilcoxon rank-sum test
#wilcox_test <- wilcox.test(group1, group2)

t_test <- t.test(group1, group2)

# Print the results
#print(wilcox_test)
print(t_test)

#t = 15.67, df = 1.7631, p-value = 0.006727

####################################
#### Sim. vs. Exp. growth rates ####
####################################

mus <- read_csv("/home/users/lzehetner/data/hek/growth_results_w_pred.csv")

new_row <- data.frame(
  state = "HEK293F",     # First column (sample names)
  mu = 0.0281,            # Experimental growth rate
  mu_predicted = 0.0301,  # Simulated growth rate
  mu_error = NA,
  timepoint = NA,
  bm = NA,
  err = NA,
  err_perc = NA,
  mu_err_perc = NA
)

# Add the new row to the dataframe 'mus'
mus <- rbind(mus, new_row)

range_all <- range(c(mus$mu, mus$mu_predicted), na.rm = TRUE)

class.labels = c()

class.labels[c(1,2,3,4)] = "#E67E22"
class.labels[c(5,6,7,8)] = "#D6EAF8"
class.labels[c(9,10,11,12)] = "#3498DB"
class.labels[c(13,14,15,16)] = "#2ECC71"
class.labels[c(17)] = "black"

class.shapes = c()

class.shapes[c(1,5,9,13)] = 21
class.shapes[c(2,6,10,14)] = 22
class.shapes[c(3,7,11,15)] = 24
class.shapes[c(4,8,12,16)] = 4
class.shapes[c(17)] = 25


png("sim_vs_exp_mu.png", width = 1200, height = 1200)

par(mar = c(8, 9, 5, 2) + 0.1)

# Adjust the mgp parameter to move the axis labels and axis titles
par(mgp = c(6, 2.5, 0))

plot(mus$mu_predicted,
     mus$mu, # x and y data
     pch=class.shapes,           # point shape
     col=class.labels, # 
     bg=class.labels,  #
     cex=4,          # point size
     #     main="PCA of HEK293 strains during rAAV production",     # title of plot
     xlab = "Sim. Growth Rates [1/h]",
     ylab = "Exp. Growth Rates [1/h]",
     cex.lab = 3,
     cex.axis = 3,
     xlim = range_all,           # Ensure equal scaling for x-axis
     ylim = range_all 
)

arrows(mus$mu_predicted, mus$mu - mus$mu_error, 
       mus$mu_predicted, mus$mu + mus$mu_error, 
       angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)

abline(a = 0, b = 1, lty = 1, col = "black", lwd = 2)  # Solid black line for y = x (perfect agreement)

# Add dashed lines for the error bounds (+10%, -10%, +25%, -25%)
abline(a = 0, b = 1.10, lty = 2, col = "black", lwd = 2)  # +10% line
abline(a = 0, b = 0.90, lty = 2, col = "black", lwd = 2)  # -10% line
abline(a = 0, b = 1.25, lty = 2, col = "black", lwd = 2)  # +25% line
abline(a = 0, b = 0.75, lty = 2, col = "black", lwd = 2)  # -25% line

# Extend the yellow shaded region all the way to the edge of the plot
polygon(c(par("usr")[1], par("usr")[2], par("usr")[2], par("usr")[1]),   # x coordinates (from plot limits)
        c(par("usr")[3] * 0.90, par("usr")[4] * 0.90, par("usr")[4] * 1.10, par("usr")[3] * 1.10),  # y coordinates (from plot limits)
        col = rgb(1, 1, 0, 0.2), border = NA)  # Light yellow shaded region

# Add labels to the right side and top for the error bounds
text(x = 0.035, y = 0.04, labels = "-10%", pos = 4, cex = 2)  # +10% 
text(x = 0.03, y = 0.04, labels = "-25%", pos = 4, cex = 2)  # +25% 
text(x = 0.041, y = 0.035, labels = "+10%", pos = 2, cex = 2)  # -10% 
text(x = 0.041, y = 0.03, labels = "+25%", pos = 2, cex = 2)  # -25% 

colors <- c("#E67E22", "#D6EAF8", "#3498DB", "#2ECC71", "black")  # Corresponding colors
shapes <- c(21, 22, 24, 4,25)  # Corresponding shapes

# Add combined color and shape legend
legend("topleft",
       legend = c("HP, TR", "HP, MO", "LP, TR", "LP, MO", "HEK293F"),  # Labels for colors
       fill = colors,  # Color legend
       title = "Conditions",  # Title for color legend
       cex = 2
)

# Add separate shape legend
legend("bottomright", 
       legend = c("4 HPT", "24 HPT", "48 HPT", "72 HPT", "Exp.Phase (HEK293F)"),  # Labels for shapes
       pch = shapes,  # Shape legend
       pt.bg = "black",  # Point fill color (can be adjusted if necessary)
       cex = 2,
       title = "Timepoints"  # Title for shape legend
)

dev.off()

####################
#### Upset plot ####
####################

upset(gpr_human1, nsets = 78, point.size = 1.5, line.size = 1, mb.ratio = c(0.15, 0.85), order.by = "freq", 
      mainbar.y.label = "Rxn Intersection", sets.x.label = "Reactions Per Model")