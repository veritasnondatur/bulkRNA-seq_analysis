############## COMPARATIVE ANALYSIS OF EMBRYONIC HINDLIMB DATA #################
# Input:.xslx file with Pbx1/2 and Hand2 mutant and ctrl E10.5 hindlimb
# Output: .xslx file with Up/Downregulated genes (with DESeq2)
# Last update of code: 05/20/2025
# AUTHOR: Vera Laub

# Load necessary libraries
library(readxl)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

# Step 1: Load bulk RNA-seq data of E10.5 hindlimb (Pbx1/2 and Hand2 mutants vs ctrl)
raw_data <- read_excel("~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/RNA-seq_Pbx1+Hand2mutant_ctrl.xlsx", 
                       sheet = "pbx1+hand2_mutants_all")

# Step 2: Prepare the data separately for all conditions
gene_data_hand2 <- raw_data %>%
  select(gene_name, starts_with("hand2"))
gene_data_pbx <- raw_data %>%
  select(gene_name, starts_with("pbx"))

# Step 3: Create a DESeq2-compatible format
# Extract count data for all samples
count_data_hand2 <- gene_data_hand2 %>%
  select(starts_with("hand2")) %>%
  as.matrix()
count_data_pbx <- gene_data_pbx %>%
  select(starts_with("pbx")) %>%
  as.matrix()

# Sample information
col_data_hand2 <- data.frame(
  condition = c("control", "control", "control", "control", 
                "Hand2_mutant", "Hand2_mutant", "Hand2_mutant"),
  row.names = colnames(count_data_hand2)
)
col_data_pbx <- data.frame(
  condition = c("control", "control", "control", 
                "Pbx1/2_mutant", "Pbx1/2_mutant", "Pbx1/2_mutant"),
  row.names = colnames(count_data_pbx)
)

# Step 4: DESeq2 analysis
dds_hand2 <- DESeqDataSetFromMatrix(countData = count_data_hand2,
                                    colData = col_data_hand2, 
                                    design = ~ condition)
dds_pbx <- DESeqDataSetFromMatrix(countData = count_data_pbx,
                                  colData = col_data_pbx, 
                                  design = ~ condition)

# Run DESeq2
dds_hand2 <- DESeq(dds_hand2)
dds_pbx <- DESeq(dds_pbx)

# Step 5: Extract results from DESeq2
res_hand2 <- results(dds_hand2, contrast = c("condition", "Hand2_mutant", "control"))
res_pbx <- results(dds_pbx, contrast = c("condition", "Pbx1/2_mutant", "control"))

# Convert DESeq2 results to a data frame
res_df_hand2 <- as.data.frame(res_hand2)
gene_info <- gene_data_hand2 %>% select(gene_name)
res_df_hand2 <- cbind(gene_info, res_df_hand2)

res_df_pbx <- as.data.frame(res_pbx)
gene_info <- gene_data_pbx %>% select(gene_name)
res_df_pbx <- cbind(gene_info, res_df_pbx)

# Adjust p-values for multiple comparisons using FDR (False Discovery Rate)
res_df_hand2$padj <- p.adjust(res_df_hand2$pvalue, method = "BH")
res_df_pbx$padj <- p.adjust(res_df_pbx$pvalue, method = "BH")

# Step 6: Combine raw counts with DESeq2 results
# Extract raw counts for E14 samples (will be used for comparison)
raw_counts_hand2 <- gene_data_hand2 %>% select(starts_with("hand2"))
raw_counts_pbx <- gene_data_pbx %>% select(starts_with("pbx"))

# Combine DESeq2 results first, then raw counts
final_results_hand2 <- cbind(res_df_hand2, raw_counts_hand2)
final_results_pbx <- cbind(res_df_pbx, raw_counts_pbx)

# Step 7: Filter significantly dysregulated genes (optional)
significant_genes_hand2 <- final_results_hand2 %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)
significant_genes_pbx <- final_results_pbx %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

# Step 8: Output results
# Save the combined table (DESeq2 results first, then raw counts)
write_xlsx(final_results_hand2, 
          path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Hand2mutant_analysis.xlsx")
write_xlsx(final_results_pbx, 
           path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Pbx1+2mutant_analysis.xlsx")

# Save the filtered significant genes (optional)
write_xlsx(significant_genes_hand2, 
           path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Hand2mutant_analysis_significant.xlsx")
write_xlsx(significant_genes_pbx, 
           path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Pbx1+2mutant_analysis_significant.xlsx")


# Step 9: Visualize results
# Plot a volcano plot to visualize upregulated vs downregulated genes
# Extract top 10 upregulated and downregulated genes based on log2FoldChange

### Plot for Hand2 mutant
top_upregulated_hand2 <- res_df_hand2 %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_downregulated_hand2 <- res_df_hand2 %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combine top upregulated and downregulated genes for labeling
top_genes <- bind_rows(
  mutate(top_upregulated_hand2, direction = "Upregulated"),
  mutate(top_downregulated_hand2, direction = "Downregulated")
)

# Create the volcano plot with gene names displayed for the top 10 upregulated and downregulated genes
volcano_plot_hand2 <- ggplot(res_df_hand2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(
    values = c("black", "red"),
    labels = c("FALSE", "TRUE")  # Custom legend labels
  ) +
  labs(
    title = "Volcano Plot of Differential Expression (Hand2 mutant, RNA-seq E10.5)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  geom_text(
    data = top_genes, 
    aes(x = log2FoldChange, y = -log10(padj), label = gene_name),
    vjust = 1, hjust = 1, size = 3, color = "blue", fontface = "italic"
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),  # White background for the plot area
    plot.background = element_rect(fill = "white"),   # White background for the whole plot
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),  # Optional: light grid lines
    panel.grid.minor = element_blank(),  # Optional: remove minor grid lines
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position.inside = c(0.8, 0.8),  # Move legend inside the plot at a specific location (top-right)
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: add a border around the legend
    legend.title = element_text(face = "bold", size = 10)  # Title for the legend, optional styling
  ) +
  guides(
    color = guide_legend(title = "Differential expression")  # Add the legend title
  )

# Save the plot to the specified folder
# Define the file path
output_file_path <- "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Hand2mutant_volcano_plot.png"

# Save the plot as a PNG file
ggsave(output_file_path, plot = volcano_plot_hand2, width = 17, height = 10, dpi = 300)

# Print the output file path
cat("Volcano plot saved to:", output_file_path, "\n")


### Plot for Pbx1/2 mutant
top_upregulated_pbx <- res_df_pbx %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top_downregulated_pbx <- res_df_pbx %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combine top upregulated and downregulated genes for labeling
top_genes <- bind_rows(
  mutate(top_upregulated_pbx, direction = "Upregulated"),
  mutate(top_downregulated_pbx, direction = "Downregulated")
)

# Create the volcano plot with gene names displayed for the top 10 upregulated and downregulated genes
volcano_plot_pbx <- ggplot(res_df_pbx, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(
    values = c("black", "red"),
    labels = c("FALSE", "TRUE")  # Custom legend labels
  ) +
  labs(
    title = "Volcano Plot of Differential Expression (Pbx1/2 mutant, RNA-seq E10.5)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  geom_text(
    data = top_genes, 
    aes(x = log2FoldChange, y = -log10(padj), label = gene_name),
    vjust = 1, hjust = 1, size = 3, color = "blue", fontface = "italic"
  ) +
  theme_minimal() +  # Use a minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white"),  # White background for the plot area
    plot.background = element_rect(fill = "white"),   # White background for the whole plot
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),  # Optional: light grid lines
    panel.grid.minor = element_blank(),  # Optional: remove minor grid lines
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position.inside = c(0.8, 0.8),  # Move legend inside the plot at a specific location (top-right)
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: add a border around the legend
    legend.title = element_text(face = "bold", size = 10)  # Title for the legend, optional styling
  ) +
  guides(
    color = guide_legend(title = "Differential expression")  # Add the legend title
  )

# Save the plot to the specified folder
# Define the file path
output_file_path <- "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Pbx1+2mutant_volcano_plot.png"

# Save the plot as a PNG file
ggsave(output_file_path, plot = volcano_plot_pbx, width = 17, height = 10, dpi = 300)

# Print the output file path
cat("Volcano plot saved to:", output_file_path, "\n")


################################################################################
#### OLD CODE WITHOUT DESEQ2

# Load required libraries
library(readxl)   # For reading Excel files
library(dplyr)    # For data manipulation
library(tidyr)    # For reshaping data
library(ggplot2)  # For plotting
library(writexl)  # To write output to excel

# === 1. Load the data ===
data <- read_excel("~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/RNA-seq_Pbx1+Hand2mutant_ctrl.xlsx", 
                   sheet = "pbx1+hand2_mutants_all")

# Ensure expression columns are numeric
data <- data %>%
  mutate(across(-1, as.numeric))

# === 2. Reshape data to long format ===
long_data <- data %>%
  pivot_longer(cols = -1,
               names_to = c("Condition", "Replicate"),
               names_sep = "_",
               values_to = "Expression") %>%
  rename(Gene = ...1)

# === 3. Summary stats ===
summary_data <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(
    mean_expression = mean(Expression, na.rm = TRUE),
    se = sd(Expression, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# === 4. Function to compute t-test p-values ===
perform_t_test <- function(data) {
  ctrl_data <- filter(data, Replicate == "ctrl")
  mutant_data <- filter(data, Replicate == "mutant")
  
  if (nrow(ctrl_data) > 1 && nrow(mutant_data) > 1) {
    if (sd(ctrl_data$Expression) > 0 && sd(mutant_data$Expression) > 0) {
      ttest_result <- t.test(ctrl_data$Expression, mutant_data$Expression)
      return(ttest_result$p.value)
    }
  }
  return(NA_real_)
}

# === 5. Compute p-values for each gene and condition ===
p_values <- long_data %>%
  group_by(Gene, Condition) %>%
  summarise(p_value = perform_t_test(cur_data()), .groups = 'drop')

# === 6. Compute log2 fold change per gene per condition ===
logfc_data <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(mean_expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Replicate, values_from = mean_expression) %>%
  mutate(log2FC = log2(mutant / ctrl)) %>%
  select(Gene, Condition, log2FC)

# === 7. Join p-values and logFC to summary_data ===
summary_data <- summary_data %>%
  left_join(p_values, by = c("Gene", "Condition")) %>%
  left_join(logfc_data, by = c("Gene", "Condition"))

# === 8. Prepare full expression + log2FC + p-value tables ===

# Get average expression values for ctrl and mutant
expression_summary <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Replicate, values_from = mean_expr) %>%
  rename(ctrl = ctrl, mutant = mutant)

# Combine with log2FC and p-values
full_results <- expression_summary %>%
  left_join(logfc_data, by = c("Gene", "Condition")) %>%
  left_join(p_values, by = c("Gene", "Condition")) %>%
  select(Gene, Condition, log2FC, p_value, ctrl, mutant)

# Subset into up/downregulated gene tables with p-value < 0.05
pbx_up <- full_results %>%
  filter(Condition == "pbx", log2FC > 0, p_value < 0.05) %>%
  arrange(desc(log2FC))

pbx_down <- full_results %>%
  filter(Condition == "pbx", log2FC < 0, p_value < 0.05) %>%
  arrange(log2FC)

hand2_up <- full_results %>%
  filter(Condition == "hand2", log2FC > 0, p_value < 0.05) %>%
  arrange(desc(log2FC))

hand2_down <- full_results %>%
  filter(Condition == "hand2", log2FC < 0, p_value < 0.05) %>%
  arrange(log2FC)


# === 9. Write everything to an Excel workbook ===
write_xlsx(
  list(
    Summary = summary_data,
    pbx_Upregulated = pbx_up,
    pbx_Downregulated = pbx_down,
    hand2_Upregulated = hand2_up,
    hand2_Downregulated = hand2_down
  ),
  path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Pbx1+Hand2_mutVSctrl_analysis.xlsx"
)



# Filter for the gene of interest
gene_of_interest <- "Hand1"  # Replace with your actual gene name
gene_data <- summary_data %>%
  filter(Gene == gene_of_interest)

# Check if gene_data is empty
if (nrow(gene_data) == 0) {
  stop("No data available for the specified gene of interest.")
}

# Create a new variable to combine Replicate and Condition
gene_data <- gene_data %>%
  mutate(Group = interaction(Replicate, Condition))

# Create the bar plot for the specific gene
ggplot(gene_data, aes(x = Group, y = mean_expression, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +  # Black border around bars
  geom_errorbar(aes(ymin = mean_expression - se, ymax = mean_expression + se), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(title = paste("Expression of", gene_of_interest, "in E10.5 hindlimb"),
       x = "Biological context",
       y = "Bulk RNA-seq mean expression [AU]") +
  theme_minimal() +
  scale_fill_manual(values = c("hand2" = "lightblue", "pbx" = "purple"), 
                    labels = c("hand2" = "Hand2", "pbx" = "Pbx")) +  # Custom colors and labels
  geom_text(data = gene_data %>% filter(Replicate == "mutant"),
            aes(label = ifelse(!is.na(p_value), paste("p =", round(p_value, 5)), "")), 
            position = position_dodge(0.5), vjust = -2, size = 3) +  # Display p-values only for mutant condition
  theme(legend.title = element_text(size = 10),  # Customize legend title size
        legend.text = element_text(size = 10),   # Customize legend text size
        legend.position = "right") +              # Position the legend on the right
  guides(fill = guide_legend(title = "Condition"))  # Add legend title

