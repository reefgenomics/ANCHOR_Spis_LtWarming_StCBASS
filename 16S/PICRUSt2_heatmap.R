
setwd("your working directory")
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)
library(ggdendro)



metadata <- read_delim("meta.txt", delim = "\t", col_names = TRUE, trim_ws = TRUE)
df <- read.table(gzfile("pred_metagenome_unstrat.tsv.gz"))
colnames(df) <- df[1,1:203]
df <- df[-1,]

df_function <- read.table(gzfile("path_abun_unstrat.tsv.gz"))
colnames(df_function) <- df_function[1,1:203]
df_function <- df_function[-1,]


names_columns <- colnames(df_function)[2:203]
metadata$sample_name <- names_columns
metadata <- metadata[,c(length(metadata),1:length(metadata)-1)]
meta_RNA <- filter(metadata, type == "RNA")
meta_DNA <- filter(metadata, type == "DNA")
meta_DNA_LT25 <- filter(meta_DNA, LT_temp == 25)
meta_RNA_T1T2 <- filter(meta_RNA,Timepoint == "T1" | Timepoint == "T2")
meta_RNA_3036 <- filter(meta_RNA, CBASS_temp == "30" |CBASS_temp == "36")
meta_RNA_30 <- filter(meta_RNA, CBASS_temp == "30")
meta_RNA_LT25 <- filter(meta_RNA_3036,LT_temp == "25")
meta_RNA_LT27 <- filter(meta_RNA_3036,LT_temp == "27")
meta_RNA_LT30 <- filter(meta_RNA_3036,LT_temp == "30")

meta_RNA_T1 <- filter(meta_RNA_3036,Timepoint == "T1")
meta_RNA_T2 <- filter(meta_RNA_3036,Timepoint == "T2")
meta_RNA_T3 <- filter(meta_RNA_3036,Timepoint == "T3")
meta_RNA_T3_30 <- filter(meta_RNA_30,Timepoint == "T3")
meta_RNA_T2_30 <- filter(meta_RNA_30,Timepoint == "T2")
meta_RNA_T1_30 <- filter(meta_RNA_30,Timepoint == "T1")




abundances_DNA_function <- df_function[,grep("DNA|pathway", x=names(df_function))]
abundances_RNA_function <- df_function[,grep("RNA|pathway", x=names(df_function))]
abundances_DNA_LT25_function <- abundances_DNA_function[,grep("LT25|pathway", x=names(abundances_DNA_function))]
abundances_RNA_function_3036 <- abundances_RNA_function[,grep("_30_|_36_|pathway", x=names(abundances_RNA_function))]
abundances_RNA_function_30 <- abundances_RNA_function[,grep("_30_|pathway", x=names(abundances_RNA_function))]
abundances_RNA_LT25_function <- abundances_RNA_function_3036[,grep("LT25|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_LT27_function <- abundances_RNA_function_3036[,grep("LT27|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_LT30_function <- abundances_RNA_function_3036[,grep("LT30|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_T1T2 <- abundances_RNA[,grep("T1S|T2S|function", x=names(abundances_RNA))]

abundances_RNA_T1_function <- abundances_RNA_function_3036[,grep("T1|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_T2_function <- abundances_RNA_function_3036[,grep("T2S|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_T3_function <- abundances_RNA_function_3036[,grep("T3S|pathway", x=names(abundances_RNA_function_3036))]
abundances_RNA_T3_30_function <- abundances_RNA_function_30[,grep("T3S|pathway", x=names(abundances_RNA_function_30))]
abundances_RNA_T2_30_function <- abundances_RNA_function_30[,grep("T2S|pathway", x=names(abundances_RNA_function_30))]
abundances_RNA_T1_30_function <- abundances_RNA_function_30[,grep("T1S|pathway", x=names(abundances_RNA_function_30))]

i <- c(2:length(abundances_RNA))
abundances_RNA[ , i] <- apply(abundances_RNA[ , i], 2,            # Specify own function within apply
                              function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_DNA_function))
abundances_DNA_function[ , i] <- apply(abundances_DNA_function[ , i], 2,            # Specify own function within apply
                                       function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_DNA_LT25_function))
abundances_DNA_LT25_function[ , i] <- apply(abundances_DNA_LT25_function[ , i], 2,            # Specify own function within apply
                                            function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_LT25_function))
abundances_RNA_LT25_function[ , i] <- apply(abundances_RNA_LT25_function[ , i], 2,            # Specify own function within apply
                                            function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_LT27_function))
abundances_RNA_LT27_function[ , i] <- apply(abundances_RNA_LT27_function[ , i], 2,            # Specify own function within apply
                                            function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_LT30_function))
abundances_RNA_LT30_function[ , i] <- apply(abundances_RNA_LT30_function[ , i], 2,            # Specify own function within apply
                                            function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T1_function))
abundances_RNA_T1_function[ , i] <- apply(abundances_RNA_T1_function[ , i], 2,            # Specify own function within apply
                                          function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T2_function))
abundances_RNA_T2_function[ , i] <- apply(abundances_RNA_T2_function[ , i], 2,            # Specify own function within apply
                                          function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T3_function))
abundances_RNA_T3_function[ , i] <- apply(abundances_RNA_T3_function[ , i], 2,            # Specify own function within apply
                                          function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T3_30_function))
abundances_RNA_T3_30_function[ , i] <- apply(abundances_RNA_T3_30_function[ , i], 2,            # Specify own function within apply
                                             function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T2_30_function))
abundances_RNA_T2_30_function[ , i] <- apply(abundances_RNA_T2_30_function[ , i], 2,            # Specify own function within apply
                                             function(x) as.numeric(as.character(x)))
i <- c(2:length(abundances_RNA_T1_30_function))
abundances_RNA_T1_30_function[ , i] <- apply(abundances_RNA_T1_30_function[ , i], 2,            # Specify own function within apply
                                             function(x) as.numeric(as.character(x)))

#write_delim(abundances_DNA, "output directory", delim = "\t")
#write_delim(abundances_RNA, "output directory", delim = "\t")
#write_delim(abundances_RNA_T1T2, "output directory", delim = "\t")
#write_delim(abundances_RNA_LT25_function, "output directory", delim = "\t")

#write_delim(meta_RNA, "output directory", delim = "\t")
#write_delim(meta_DNA, "output directory", delim = "\t")
#write_delim(meta_RNA_T1T2, "output directory", delim = "\t")







#### Heatmap
class(abundances_RNA_T3_30_function)
row.names(abundances_RNA_T3_30_function) <- NULL
daa_results_df_RNA_T3_30_fdr <- pathway_daa(abundance = abundances_RNA_T3_30_function %>% column_to_rownames("pathway"), metadata = meta_RNA_T3_30, group = "LT_temp", daa_method = "DESeq2", p.adjust = "fdr")
daa_annotated_results_df_RNA_T3_30_fdr <- pathway_annotation(pathway = "MetaCyc", daa_results_df = daa_results_df_RNA_T3_30_fdr, ko_to_kegg = FALSE)

feature_with_p_0.01_RNA_T3_30_fdr <- daa_annotated_results_df_RNA_T3_30_fdr %>% filter(p_adjust < 0.01)


test_df_3 <- feature_with_p_0.01_RNA_T3_30_fdr |> 
  mutate(Combinations = ifelse(group1 == 25 & group2 ==27, "Combination_1", ifelse(group1 == 25 & group2 ==30, "Combination_2", "Combination_3" )))



Df1_T3 <- abundances_RNA_T3_30_function
Df2_T3 <- test_df_T3


for (i in 1:nrow(Df2_T3)) {
  feature_val <- Df2_T3$feature[i]
  combination_val <- Df2_T3$Combinations[i]
  other_rows <- Df2_T3[Df2_T3$feature == feature_val, ]
  unique_combinations <- unique(other_rows$Combinations)
  
  if (length(unique_combinations) == 1) {
    if (combination_val == "Combination_1") {
      pattern <- "LT30_"
    } else if (combination_val == "Combination_2") {
      pattern <- "LT27_"
    } else if (combination_val == "Combination_3") {
      pattern <- "LT25_"
    }
    Df1_T3[Df1_T3$pathway == feature_val, grep(pattern, names(Df1_T3))] <- NA
  }
}


pathway_heatmap_2 <- function(abundance, metadata, group, description_df) {
  # Heatmaps use color changes to visualize changes in values. However, if the
  # data for plotting the heat map are too different, for example, if the heat
  # map is plotted using gene expression data, gene1 is expressed above 1000 in
  # all samples and gene2 is expressed between 1-10 in all samples, it is
  # difficult to plot the heat map with small changes in the expression of two
  # genes in different samples by the colors to reflect. Therefore, when
  # plotting a heat map, we usually normalize the gene expression data, that
  # is, we subtract the mean value of each gene expression from the expression
  # of this gene in all samples and divide it by its standard deviation, and
  # this normalization is called standard normalization or Z-score processing.
  # The processed values are reduced equally, and the expression of each gene
  # in all samples becomes a set of values with a mean of 0 and a standard
  # deviation of 1. At this point, the plotted heat map gives a good indication
  # of the variation in expression of all genes across samples.
  
  # Check that 'group' is a column in 'metadata'
  if (!group %in% colnames(metadata)) {
    stop(paste("group:", group, "must be a column in metadata"))
  }
  
  # Find the column in metadata that matches the column names of abundance
  sample_name_col <- colnames(metadata)[sapply(colnames(metadata), function(x) all(colnames(abundance) %in% metadata[[x]]))]
  metadata$sample_name <- metadata %>% select(all_of(c(sample_name_col))) %>% pull()
  
  if (!all(colnames(abundance) %in% metadata$sample_name)) {
    stop("Samples in abundance and metadata must match")
  }
  
  # Now sample_name_col contains the column name in metadata that stores the sample names
  
  z_abundance <- t(apply(abundance, 1, scale))
  colnames(z_abundance) <- colnames(abundance)
  
  # Convert the abundance matrix to a data frame
  z_df <- as.data.frame(z_abundance)
  metadata <- metadata %>% as.data.frame()
  
  # Order the samples based on the environment information
  ordered_metadata <- metadata[order(metadata[, group]),]
  ordered_sample_names <- ordered_metadata$sample_name
  order <- ordered_metadata$sample_name
  ordered_group_levels <- ordered_metadata %>% select(all_of(c(group))) %>% pull()
  
  
  # Convert the abundance data frame to a long format
  long_df <- z_df %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(cols = -rowname,
                        names_to = "Sample",
                        values_to = "Value") %>% 
    left_join(metadata %>% select(all_of(c("sample_name",group))), by = c("Sample" = "sample_name")) %>%
    mutate(CBASS_temp = ifelse(grepl("_30_", Sample), "30", NA),
           # Timepoint = ifelse(grepl("T1S", Sample), "T1", ifelse(grepl("T2S", Sample), "T2", "T3")),
           LT_temp = ifelse(grepl("LT25", Sample), "25", ifelse(grepl("LT27", Sample), "27", "30"))) |> 
    group_by(LT_temp, rowname) |> 
    mutate(Mean_Value = mean(Value),  
           Value = NULL) |> 
    ungroup() 
  
  # Set the order of the samples in the heatmap
  long_df$Sample <- factor(long_df$Sample, levels = order)
  long_df <- long_df %>%
    left_join(select(description_df,feature,  description), by = c("rowname" = "feature"))
  long_df <- long_df %>%
    mutate(Sample = NULL) |> 
    distinct()
  long_df <- arrange(long_df, rowname) 
  long_df
  
  
  
  # Dendogram
  
  long_df_scaled <- long_df
  long_df_scaled[, 4] <- scale(long_df[, 4])
  long_df_matrix <- as.matrix(dist(long_df_scaled[, 4], method = "euclidean"))
  rownames(long_df_matrix) <- long_df_scaled$description
  long_df_matrix[is.na(long_df_matrix)] <- 0
  long_df_dendro <- as.dendrogram(hclust(d = dist(x = long_df_matrix)))
  long_df_order <- order.dendrogram(long_df_dendro)
  long_df$description <- factor(
    x = long_df$description,
    levels = unique(long_df_scaled$description[long_df_order]),
    ordered = TRUE
  )
  

  
  # Compute breaks from the data
  breaks <- range(long_df$Value, na.rm = TRUE)
  
  colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")
  
  # Create the heatmap using ggplot
  p <-
    ggplot2::ggplot(data = long_df,
                    mapping = ggplot2::aes(x = LT_temp, y = description, fill = Mean_Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "darkslateblue", mid = "floralwhite", high = "darkred", midpoint = 0, limits = c(-1.5, 1.5)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_y_discrete(expand = c(0, 0), position = "left") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    # Customize the appearance of the heatmap
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      panel.spacing = unit(0, "lines"),
      legend.title = ggplot2::element_text(size = 12, color = "black",face = "bold"),
      legend.text = ggplot2::element_text(size = 12, color = "black",face = "bold"),
      panel.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(l = 0, unit = "cm"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    # Add a color bar to the heatmap
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        direction = "vertical",
        reverse = F,
        barwidth = unit(0.6, "cm"),
        barheight = unit(10, "cm"),
        title = "Z Score",
        title.position = "top",
        title.hjust = -1,
        ticks = TRUE,
        label = TRUE
      )
    ) + ggh4x::facet_nested(cols = vars(!!sym(group)), space = "free", scale = "free", switch = NULL, strip =strip_nested(background_x = elem_list_rect(color = "black")))
  # Print the ordered sample names and group levels (if needed)
  cat("The Sample Names in order from left to right are:\n")
  cat(ordered_sample_names, sep = ", ")
  cat("\n")
  
  cat("The Group Levels in order from left to right are:\n")
  cat(ordered_group_levels, sep = ", ")
  cat("\n")
  
  return(p)
}





T3 <- pathway_heatmap_2(abundance = Df1_T3 %>% filter(pathway %in% feature_with_p_0.01_RNA_T3_30_fdr$feature) %>% 
                          column_to_rownames("pathway"), 
                        metadata = meta_RNA_T3_30, 
                        group = "CBASS_temp", 
                        description_df = daa_annotated_results_df_RNA_T3_30_fdr)


T1 + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title = element_text(size = 20),
        panel.spacing = unit(1, "lines"),
        strip.background =element_rect(fill="white")) 


