############################################################
##              reading the data set                     ##
############################################################
rm(list = ls(all.names = TRUE))
gc()

#install.packages("tidyverse")
library(tidyverse)

base_dir <- "/Users/dabkek/Library/CloudStorage/Box-Box/1_Primary_Recurrent_Project/DATA/Akoya_data_images_91123/"

df <- read.table(paste0(base_dir, "AKOYA_Data_ROI_specific_92223/summarized_results/merged_cell_summary_akoya.txt"), sep = "\t", header = F) 
colnames(df) <- c("cell_type", "raw_count", "norm_count", "accession_number")

metadata <- read_csv(paste0(base_dir, "Akoya_sample_info_25slides_072023.csv"))
metadata <- as.data.frame(metadata)

setdiff(metadata$accession_number, df$accession_number)
setdiff(df$accession_number, metadata$accession_number)

# S13-32894-G1 is labelled as S13-32894_G1 in the merged file
# changing this specific accession number in the metadata file

metadata[grep("S13-32894-G1", metadata$accession_number),1] <- "S13-32894_G1"
setdiff(metadata$accession_number, df$accession_number)
setdiff(df$accession_number, metadata$accession_number)

df <- left_join(df, metadata)
sum(is.na(df))

### converting normalized counts (cell count/ total number of cells) to percentages
df$norm_count <- df$norm_count*100

##########################################################################
### Samples did not have any tumor specific ROI: 05S-14310_D; S09-38573_C1
##########################################################################

##########################################################################
### basic cell type summary
##########################################################################
which(table(df$cell_type) == 23)
# all samples have the following cell types: Cytotoxic T cells, Endothelial cells, Helper T, M1, M2, other, smooth muscle, Tregs, Tumor
which(table(df$cell_type) == 1)
# cell type only identified in 1 sample: T cells double negative (24487_1923_Lt_ovary)

#########################################################
### bar plots of cell type percentages
#########################################################

cell_type_col <-  c("B cells" ="#0000FF", "Cytotoxic T cells" ="#FF0000", "Endothelial cells" ="#00FF00",
                    "Fibroblasts" ="#000033","Helper T cells"= "#FF00B6", "M1 Macrophages"= "#005300",
                    "M2 Macrophages" ="#FFD300","Necrotic Tumor cells"= "#009FFF", "other"="#9A4D42",
                    "Smooth muscle cells"="#00FFBE", "Tregs"="#783FC1","Tumor cells"= "#1F9698",
                    "Dendritic cells"="#FFACFD", "NK cells"="#B1CC71","Monocytes"= "#F1085C",
                    "T cells double negative"="#FE8F42", "Epithelial cells"="#DD00FF")

png("./Figures/quant23_bar_plots.png", units = "in", height = 8, width = 12, res = 600)
ggplot(df, aes(x = POCROC_labels, y = norm_count, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


df$prim_BRCA <- paste(df$prim_recur, df$germline_BRCA, sep = "-")

df$HRD_cat <- "BRCAmut"
df$HRD_cat <- ifelse(df$germline_BRCA == "BRCAwt", "BRCAwt",df$HRD_cat)

df$HRD_cat <- ifelse(df$HRD_status == "HRD" & df$germline_BRCA == "BRCAwt",
                                 "other_HRD", df$HRD_cat)
df$prim_HRD <- paste(df$prim_recur, df$HRD_status, sep = "-")

df <- df[order(df$prim_HRD),]
##########################################################################
### cell type ratios calculation
##########################################################################

sample <- unique(df$accession_number)

temp <- df

for (i in 1:23) {
  temp[grep(sample[i], temp$accession_number),3] <- 
    (filter(df, accession_number == sample[i]))$norm_count/ (filter(df, cell_type == "Tumor cells") %>% filter(., accession_number == sample[i]))$norm_count
}

df_tumor_ratio <- temp

sum(is.na(df_tumor_ratio))

##########################################################################
### total percentage cell type comparisons (bassem recommended not to do this)
##########################################################################
group4_comparison <- list( c("Primary-HRD", "Primary-HRP"),
                           c("Recurrent-HRD", "Recurrent-HRP"))

group4_comparison_germline <- list( c("Primary-BRCAmut", "Primary-BRCAwt"),
                                    c("Recurrent-BRCAmut", "Recurrent-BRCAwt"))

group2_HRD <- list( c("HRD", "HRP"))
group2_prim_recur <- list( c("Primary", "Recurrent"))


group4_color <- c("#AA323F",'#286D97', "#E99DA5",'#96C2DD')
HRD_color <- c("#AA323F",'#286D97')
recur_color <- c("#f97c03",'#5e4fa2')

library(ggpubr)

plot_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("Cytotoxic T cells", "Endothelial cells", "Helper T cells",
                                         "M1 Macrophages", "M2 Macrophages", "Smooth muscle cells",
                                         "Tregs", "Tumor cells")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 100) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,110))
  
}

png("./Figures/prim_HRD_quant23_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_function(df, "prim_HRD", group4_color, group4_comparison)
dev.off()


png("./Figures/prim_HRD_quant23_recuronly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_function(df, "prim_recur", recur_color, group2_prim_recur)
dev.off()

png("./Figures/prim_HRD_quant23_HRDonly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_function(df, "HRD_status", HRD_color, group2_HRD)
dev.off()


plot_ratio_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("Cytotoxic T cells", "Endothelial cells", "Helper T cells",
                                           "M1 Macrophages", "M2 Macrophages", "Smooth muscle cells",
                                           "Tregs")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 1) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("cell type Ratio") +
    labs(fill = "")+ ylim(c(0.0000,2))
  
}

png("./Figures/prim_HRD_quant23_cell_ratio_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_ratio_function(df_tumor_ratio, "prim_HRD",group4_color, group4_comparison)
dev.off()


png("./Figures/prim_HRD_quant23_cell_ratio_Recuronly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_ratio_function(df_tumor_ratio, "prim_recur",recur_color, group2_prim_recur)
dev.off()

png("./Figures/prim_HRD_quant23_cell_ratio_HRDonly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_ratio_function(df_tumor_ratio, "HRD_status",HRD_color, group2_HRD)
dev.off()



##########################################################################
### germline comparisons
##########################################################################


plot_germline_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("Cytotoxic T cells", "Endothelial cells", "Helper T cells",
                                           "M1 Macrophages", "M2 Macrophages", "Smooth muscle cells",
                                           "Tregs", "Tumor cells")) %>% filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 100) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,110))
  
}


png("./Figures/prim_HRD_quant23_germline_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_germline_function(df, "prim_BRCA", group4_color, group4_comparison_germline)
dev.off()


plot_germline_ratio_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("Cytotoxic T cells", "Endothelial cells", "Helper T cells",
                                           "M1 Macrophages", "M2 Macrophages", "Smooth muscle cells",
                                           "Tregs")) %>% filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 1) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("cell type Ratio") +
    labs(fill = "")+ ylim(c(0.0000,2))
  
}

png("./Figures/prim_HRD_quant23_cell_ratio_germline_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_germline_ratio_function(df_tumor_ratio, "prim_BRCA",group4_color, group4_comparison_germline)
dev.off()


##########################################################################
### low quant cell types
##########################################################################
plot_lowquant_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("B cells", "Dendritic cells", "Fibroblasts",
                                           "Monocytes", "NK cells")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 50) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,60))
  
}

png("./Figures/prim_HRD_quant23_lowQuant_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_function(df, "prim_HRD", group4_color, group4_comparison)
dev.off()


png("./Figures/prim_HRD_quant23_lowQuant_recuronly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_function(df, "prim_recur", recur_color, group2_prim_recur)
dev.off()

png("./Figures/prim_HRD_quant23_lowquant_HRDonly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_function(df, "HRD_status", HRD_color, group2_HRD)
dev.off()


plot_lowquant_ratio_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("B cells", "Dendritic cells", "Fibroblasts",
                                           "Monocytes", "NK cells")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 9) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Cell type ratio") +
    labs(fill = "")+ ylim(c(0.0000,10))
  
}

png("./Figures/prim_HRD_quant23_lowquant_cell_ratio_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_ratio_function(df_tumor_ratio, "prim_HRD",group4_color, group4_comparison)
dev.off()


png("./Figures/prim_HRD_quant23_lowquant_cell_ratio_Recuronly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_ratio_function(df_tumor_ratio, "prim_recur",recur_color, group2_prim_recur)
dev.off()

png("./Figures/prim_HRD_quant23_lowquant_cell_ratio_HRDonly_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_ratio_function(df_tumor_ratio, "HRD_status",HRD_color, group2_HRD)
dev.off()

##########################################################################
### germline comparisons; low quant cell types
##########################################################################


plot_lowquant_germline_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("B cells", "Dendritic cells", "Fibroblasts",
                                           "Monocytes", "NK cells")) %>% 
    filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 50) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,60))
  
}


png("./Figures/prim_HRD_quant23_lowquant_germline_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_germline_function(df, "prim_BRCA", group4_color, group4_comparison_germline)
dev.off()


plot_lowquant_germline_ratio_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("B cells", "Dendritic cells", "Fibroblasts",
                                           "Monocytes", "NK cells")) %>% 
    filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 8) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("cell type Ratio") +
    labs(fill = "")+ ylim(c(0.0000,10))
  
}

png("./Figures/prim_HRD_quant23_lowquant_cell_ratio_germline_ttest.png", units = "in", height = 11, width = 11, res = 600)
plot_lowquant_germline_ratio_function(df_tumor_ratio, "prim_BRCA",group4_color, group4_comparison_germline)
dev.off()


#### have not plotted the following 4 cells types: Epithelial cells, Necrotic Tumor cells, T cell double negative, and other

##########################################################################
### boxplot of broad categories; total perce
##########################################################################

df$broad_cells <- "Immune_cells"
df$broad_cells <- ifelse(df$cell_type == "Endothelial cells", "Stroma", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "Fibroblasts", "Stroma", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "Necrotic Tumor cells", "Necrotic Tumor cells", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "other", "other", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "Smooth muscle cells", "Stroma", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "Tumor cells", "Tumor", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "T cells double negative", "T cells double negative", df$broad_cells)
df$broad_cells <- ifelse(df$cell_type == "Epithelial cells", "Epithelial cells", df$broad_cells)

unique(df$broad_cells)


plot_broad_cat <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., broad_cells %in% c("Immune_cells", "Stroma", "Tumor",
                                           "other")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "broad_cells")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 100) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,110))
  
}

png("./Figures/quant23_broadcells_ttest.png", units = "in", height = 8, width = 8, res = 600)
plot_broad_cat(df, "prim_HRD",group4_color, group4_comparison)
dev.off()

##########################################################################
### boxplot of broad categories; cell ratios
##########################################################################
plot_cellratio_broad_cat <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., broad_cells %in% c("Immune_cells", "Stroma",
                                             "other")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "broad_cells")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Cell type ratio")
    #labs(fill = "")+ ylim(c(0.0000,2))
  
}

png("./Figures/quant23_broadcells_cellRatio_HRDonly_ttest.png", units = "in", height = 8, width = 8, res = 600)
plot_cellratio_broad_cat(df_tumor_ratio, "HRD_status",HRD_color, group2_HRD)
dev.off()

##########################################################################
### boxplot of broad categories; total percentage: germline
##########################################################################

plot_broad_germline_cat <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., broad_cells %in% c("Immune_cells", "Stroma", "Tumor",
                                             "other")) %>% filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "broad_cells")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 100) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0.0000,110))
  
}

png("./Figures/quant23_broadcells_germline_ttest.png", units = "in", height = 8, width = 8, res = 600)
plot_broad_germline_cat(df, "prim_BRCA",group4_color, group4_comparison_germline)
dev.off()

plot_cellratio_broad_germline_cat <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., broad_cells %in% c("Immune_cells", "Stroma",
                                             "other")) %>% filter(., !HRD_cat == "other_HRD") %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "broad_cells")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Cell type ratio")
    #labs(fill = "")+ ylim(c(0.0000,110))
  
}

png("./Figures/quant23_cellratio_broadcells_germline_ttest.png", units = "in", height = 8, width = 8, res = 600)
plot_cellratio_broad_germline_cat(df_tumor_ratio, "prim_BRCA",group4_colsor, group4_comparison_germline)
dev.off()

#########################################################
### ERAP1 and sting intensity box plots
#########################################################
# subsetted out the tumor cells and their erap1 and sting expression
# calculated the average of erap1 and sting intensity of all tumor cells in each sample
# so each point is a erap1 or sting average intensity of all cells in one sample

merged_data <- read.csv("merged_ERAP1_STING_avg_intensity.csv", header = F)
colnames(merged_data) <- c("ERAP1", "STING", "Sample")
merged_data <- separate(merged_data, Sample, into = c("one", "two", "three"),
                        sep = "_")

merged_data$accession_number <- paste(merged_data$one, merged_data$two, sep = "-")
merged_data$accession_number <- paste(merged_data$accession_number, merged_data$three, sep = "_")

setdiff(metadata$accession_number, merged_data$accession_number)
setdiff(merged_data$accession_number, metadata$accession_number)


##########################################################################
### removing the spleen sample: 05S-14310_D; S09-38573_C1
##########################################################################
merged_data <- na.omit(merged_data)


merged_data <- left_join(merged_data, metadata)

merged_data$prim_HRD <- paste(merged_data$prim_recur, merged_data$HRD_status, sep = "-")


merged_data <- merged_data[order(merged_data$prim_HRD),]

p <- ggboxplot(merged_data, x = "prim_HRD", y = "ERAP1",
               fill = "prim_HRD", 
               palette = group4_color,
               add = "jitter")

png("./Figures/prim_HRD_quant23_ERAP1_ttest.png", units = "in", height = 6, width = 8, res = 600)
p + stat_compare_means(comparisons = group4_comparison, 
                       method = "t.test",  size = 8, ylab = 1) +
  theme(text = element_text(size = 15),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  #panel.border = element_rect(colour = "black", fill=NA, size=0.5) +
  ylab("ERAP1 avg intensity") +
  labs(fill = "")+ ylim(c(-1,2))
dev.off()


##########################################################################
### re-plot for Michelle's presentation; M1 and helper T cells
##########################################################################


plot_function <- function(data, x_axis, color, comparisons){
  p <- data %>% filter(., cell_type %in% c("M1 Macrophages")) %>% 
    ggboxplot(., x = paste0(x_axis), y = "norm_count",
              fill = paste0(x_axis), 
              palette = color,
              add = "jitter",
              facet.by = "cell_type")
  # Use only p.format as label. Remove method name.
  p + stat_compare_means(comparisons = comparisons, 
                         method = "t.test",  size = 8, ylab = 11) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Total percentage cell type") +
    labs(fill = "")+ ylim(c(0, 12))
  
}

png("./Figures/prim_HRD_quant23_M1_ttest.png", units = "in", height = 6, width = 8, res = 600)
plot_function(df, "prim_HRD", group4_color, group4_comparison)
dev.off()


##########################################################################
### intensity plot- average cell type
##########################################################################

avg_celltype_intensity <- read.csv("merged_celltype_intensity_summary_Tumoronly.txt", header = F)

colnames(avg_celltype_intensity) <- c("cell_type", "avg_norm_intensity", "Sample")
avg_celltype_intensity <- separate(avg_celltype_intensity, Sample, into = c("one", "two", "three"),
                        sep = "_")

avg_celltype_intensity$accession_number <- paste(avg_celltype_intensity$one, avg_celltype_intensity$two, sep = "-")
avg_celltype_intensity$accession_number <- paste(avg_celltype_intensity$accession_number, avg_celltype_intensity$three, sep = "_")

setdiff(metadata$accession_number, merged_data$accession_number)
setdiff(merged_data$accession_number, metadata$accession_number)

avg_celltype_intensity <- left_join(avg_celltype_intensity, metadata)

distinct_colors <-  c("#0000FF", "#FF0000", "#00FF00","#000033", "#FF00B6",  "#005300",
                    "#FFD300","#009FFF", "#9A4D42","#00FFBE", "#783FC1", "#1F9698",
                    "#FFACFD", "#B1CC71","#F1085C", "#FE8F42", "#DD00FF","#204254",
                    "#720055", "#766C95", "#02AD24","#C8FF00", "#886C00","#FFB79F",
                    "#858567", "#A10300", "#a76b56")
names(distinct_colors) <- unique(avg_celltype_intensity$cell_type)


avg_celltype_intensity$prim_HRD <- paste(avg_celltype_intensity$prim_recur, avg_celltype_intensity$HRD_status, sep = "-")
avg_celltype_intensity <- avg_celltype_intensity[order(avg_celltype_intensity$prim_HRD),]

plot_list <- list()
for (i in 1:27) {
  p <- avg_celltype_intensity %>% filter(., cell_type == names(distinct_colors)[i]) %>% 
    ggboxplot(., x = "prim_HRD", y = "avg_norm_intensity",
              fill = "prim_HRD", 
              palette = group4_color,
              add = "jitter") + 
    stat_compare_means(comparisons = group4_comparison, 
                         method = "t.test", size = 8) +
    theme(text = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "top",
          plot.margin = margin(1, 1, 1, 1, "cm")) +
    #panel.border = linewidth(colour = "black", fill=NA, size=0.5) +
    ylab(paste0(names(distinct_colors)[i], " Avg intensity(Tumor)"))+
  labs(fill = "")
  plot_list[[i]] <- p
}

pdf("plots.pdf")
for (i in 1:27) {
  print(plot_list[[i]])
}
dev.off()



p <- avg_celltype_intensity %>% filter(., cell_type == names(distinct_colors)[3]) %>% 
  ggboxplot(., x = "prim_recur", y = "avg_norm_intensity",
            fill = "prim_recur", 
            palette = recur_color,
            add = "jitter") + 
  stat_compare_means(comparisons = group2_prim_recur, 
                     method = "t.test", size = 8) +
  theme(text = element_text(size = 15),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  #panel.border = linewidth(colour = "black", fill=NA, size=0.5) +
  ylab(paste0(names(distinct_colors)[3], " Avg intensity(Tumor)"))+
  labs(fill = "")

















