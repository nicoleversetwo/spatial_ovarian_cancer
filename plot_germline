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
