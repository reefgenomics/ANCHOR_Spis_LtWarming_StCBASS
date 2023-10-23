##############################
### Transforming table phyloseq object
##############################
setwd("your working directory")
library(ggthemes)
library(microViz)
library(cowplot)
library(patchwork)
library(devtools)
library(stringr)
library(pheatmap)
library(openxlsx)
library(microViz)
library(DESeq2)
library(microbiome)
library(grafify)
library(phylosmith)
library(devtools)
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(pairwiseAdonis)
library(rstatix)
library(ggVennDiagram)
library(biomformat)
##############################
### load file
##############################
ASV_500 <- read.table("ASV_500.txt", header = TRUE, sep = "\t") # read in the data
colnames(ASV_500) <- gsub('_4wkLT25_', '_._', colnames(ASV_500), fixed=TRUE)
colnames(ASV_500) <- gsub('_._', '_', colnames(ASV_500), fixed=TRUE)

##############################
### seperate column names
##############################
col_asv_names <- colnames(ASV_500[2:204]) # get the labels of the samples
sep_names <- strsplit(col_asv_names, split = "_") # split the labels based on "_"

##############################
### turn list of separated names into data frame
##############################
df_sep_names <- as.data.frame(do.call(rbind, sep_names)) # turns the list into a data frame

##############################
### use regex to transform df
##############################
Timepoint <- str_extract(df_sep_names$V1, "T\\d") # get the timep oint based on the pattern T followed by any digit
df_sep_names$Timepoint <- Timepoint #create a column called timepoint in the data frame
df_sep_names$V1 <- str_remove(df_sep_names$V1, "T\\d") # remove the time point variable based on the same pattern as before 
df_sep_names$LT_temp <- str_extract(df_sep_names$V2, "\\d\\d\\b") 
# create a ne column called LT_temp with the values extracted from column V2 based on the pattern 
#"any digit (\\d)","any digit (\\d)", "found at the end of the word (\\b)"
df_sep_names$V2 <- str_remove(df_sep_names$V2, "\\d\\d\\b") # remove the LT temperature from column V2 based on the pattern before
df_sep_names$V2 <- str_remove(df_sep_names$V2, "wkLT") # remove wkLT from every row in column V2
df_sep_names$V4 <- str_remove(df_sep_names$V4, "16S") # remove 16S from every row in column V4
sorted_names <- c("Genotype", "Weeks", "CBASS_temp", "type", "ID", "Timepoint", "LT_temp") # create new column names
colnames(df_sep_names) <- sorted_names # exchange column names
df_sep_names[141:203,2] <- str_replace(df_sep_names[141:203,2], "24", "28") 
# change the value 24 to 28 in column weeks for timepoint T§

##############################
### change column order + add original labels as saftey check
##############################
col_order <- c("Timepoint", "Genotype", "Weeks", "LT_temp", "CBASS_temp", "type", "ID") 
# realize that the order is not as the order in the original sample and the maintain the order switch th column names
df_sep_names <- df_sep_names[, col_order]
rownames(df_sep_names) <- col_asv_names



##############################
### remove everything except df_sep_names and ASV_500
##############################
rm(list=setdiff(ls(), c("df_sep_names", "ASV_500", "col_asv_names")))



########################################################################
########################################################################
########################################################################
##### ASV TBALE AS PHYLOSEQ
########################################################################
########################################################################
########################################################################

########################################################################
### create table for ASVs only
########################################################################

df <- data.frame(ASV_500[,1:204]) #only abundances of ASVs
df <- t(df) #transpose (switch rows and columns)
df <- data.frame(df) # as data frame
colnames(df) <- df[1, ] # first row as column names
df <- df[-1,] # delete first row
df <- as.data.frame(lapply(df,as.numeric)) # change everything in ASV abundance table into numeric
rownames(df) <- col_asv_names

########################################################################
### create taxonomy table
########################################################################
taxon <- ASV_500[,c("Row.names", "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species")] # pick columns with taxonomy
taxon2 <- taxon[,-1] #create taxon2 without the first column
rownames(taxon2) <- taxon[,1] # add first column of taxon as row names of taxon2
taxon <- as.matrix(taxon2) # taxon2 becomes taxon in form of a matrix

ASV <- otu_table(df, taxa_are_rows = FALSE) # create otu_table (here ASV) for phyloseq with taxa as columns
TAX <- tax_table(taxon) # create taxonomy table for phyloseq
samples <- sample_data(df_sep_names) # create sample table for phyloseq
phy <- phyloseq(ASV, TAX, samples) # put all together to create a phyloseq object
phy <- prune_samples(sample_sums(phy) >= 5000, phy) # remove samples with less than 5000 reads
rm(list=setdiff(ls(), c("phy")))
View(otu_table(phy))
########################################################################
####### Rarefy ######
########################################################################

rare_phy <- rarefy_even_depth(phy, sample.size = min(sample_sums(phy)), rngseed=123)
rm(list=setdiff(ls(), c("phy", "rare_phy")))



########################################################################
### be happy that it worked and try stuff
########################################################################


########################################################################
### create a color palette for the 20 most abundant families per type and time point
########################################################################


### Get a vector with the 20 most relative abundant families for 16SrDNA
phy.type <- subset_samples(rare_phy, type == "DNA")
phy.type.relative <- transform_sample_counts(phy.type, function(x) x / sum(x))
phy.family <- tax_glom(phy.type.relative, taxrank = "Family")
N <- 20
top <- names(sort(taxa_sums(phy.family), decreasing = TRUE))[1:N]
phy.family.relative.top <- prune_taxa(top, phy.family)
top20.tax <- as.data.frame(phy.family.relative.top@tax_table)
top20.families.DNA <- top20.tax$Family


### Get a vector for each time point (T1, T2, T3) of the 20 most relative abundant families for 16srRNA
phy.type <- subset_samples(rare_phy, type == "RNA")
phy.type.time <- subset_samples(phy.type, Timepoint == "T1") 
# change the time point T1, T2 or T3
phy.type.relative <- transform_sample_counts(phy.type.time, function(x) x / sum(x) )
phy.family <- tax_glom(phy.type.relative, taxrank = "Family")
top <- names(sort(taxa_sums(phy.family), decreasing = TRUE))[1:N]
test.top.RNA <- prune_taxa(top, phy.family)
T1.top.20 <- as.data.frame(test.top.RNA@tax_table) 
# change the time point T1, T2 or T3
T1_family <- T1.top.20$Family
T2_family <- T2.top.20$Family
T3_family <- T3.top.20$Family

overall.top20.families <- unique(c(T1_family, T2_family, T3_family, top20.families.DNA))
# create a vector with all unique families from DNA and all time points
col_vector_families <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                         "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                         "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                         "#E5C494", "#B3B3CC", "#8DD3C7", "#BEBADA", "#FDCDAC",
                         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#AEBBD2",
                         "#FFBE7D", "#C6DBEF", "#FDB462", "#FCCDE5")

### create a palette with the chosen colors ###
palette_family <- col_vector_families
names(palette_family) <- overall.top20.families
palette_family


########################################################################
### create a color palette for the 20 most abundant ASVs per type and time point
########################################################################

### Get a vector with the 20 most relative abundant ASVs for 16S rDNA
phy.type <- subset_samples(rare_phy, type == "DNA")
phy.type.relative <- transform_sample_counts(phy.type, function(x) x / sum(x))
N <- 20
top <- names(sort(taxa_sums(phy.type.relative), decreasing = TRUE))[1:N]
phy.family.relative.top <- prune_taxa(top, phy.type.relative)
top20.tax <- as.data.frame(phy.family.relative.top@tax_table)
top20.ASVs.DNA <- row.names(top20.tax)

### Get a vector for each time point (T1, T2, T3) of the 20 most relative abundant ASVS for 16S rRNA
phy.type <- subset_samples(rare_phy, type == "RNA")
phy.type.time <- subset_samples(phy.type, Timepoint == "T3") 
# change the time point T1, T2 or T3
phy.type.relative <- transform_sample_counts(phy.type.time, function(x) x / sum(x) )
top <- names(sort(taxa_sums(phy.type.relative), decreasing = TRUE))[1:N]
test.top.RNA <- prune_taxa(top, phy.type.relative)
T3.top.20 <- as.data.frame(test.top.RNA@tax_table) 
# change the time point T1, T2 or T3
T1_ASV <- row.names(T1.top.20)
T2_ASV <- row.names(T2.top.20)
T3_ASV <- row.names(T3.top.20)


overall.top20.ASVs <- unique(c(T1_ASV, T2_ASV, T3_ASV, top20.ASVs.DNA))
length(overall.top20.ASVs)
col_vector_ASVs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                     "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                     "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                     "#E5C494", "#B3B3CC", "#8DD3C7", "#BEBADA", "#FDCDAC",
                     "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#AEBBD2",
                     "#FFBE7D", "#C6DBEF", "#FDB462", "#FCCDE5", "#FF0000", 
                     "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF",
                     "#800000", "#008000", "#000080", "#808000", "#800080", 
                     "#008080", "#FFA500", "#FFC0CB", "#800080", "#008080", 
                     "#FFD700", "#001080", "#FF4500", "#FF1493", "#483D8B", 
                     "#2E8B57", "#00CED1", "#00FA9A", "grey" )
palette_ASV <- col_vector_ASVs
names(palette_ASV) <- overall.top20.ASVs
palette_ASV
rm(list=setdiff(ls(), c("phy", "palette_family", "palette_ASV", "rare_phy")))


########################################################################
### function to easy filter and then plot the 20 most abundant families
########################################################################

find.top.family.CBASS <- function(data, filter1, filter2){ 
  # data = original phyloseq object, type = RNA or DNA, temperature = LT_temp. But they can be anything you want to filter for
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1))) 
  # I don't know why it didn't work with just subset_samples(data, parameter) but this workaround works
  # do.call constructs and executes a function call from a name or a function and a list of arguments to be passed to it.
  # the list here is, the phyloseq object followed by the filter condition
  # Subset all samples of your wished type from the original phyloseq object. Here DNA or RNA
  phy.type.time <- do.call("subset_samples", list(quote(phy.type), substitute(filter2))) 
  # Subset all samples from the already filtered phyloseq object. Here Timepoint = timepoint
  phy.type.relative <- transform_sample_counts(phy.type.time, function(x) x / sum(x) )
  # transform_sample_counts takes as arguments a phyloseq object and an R function and 
  # returns a phyloseq-object in which the abundance values have been transformed
  # Here we transform the abundance values into relative abundance
  # If you want it in percent function(x) 100*x / sum(x)
  #phy.type.time <- do.call("subset_samples", list(quote(phy.type.relative), substitute(timepoint)))
  # Subset all samples from the already filtered phyloseq object. Here Timepoint = timepoint
  phy.family <- tax_glom(phy.type.relative, taxrank = "Family")
  # This method merges species that have the same taxonomy at a certain taxonomic rank. Here we use the family level
  # Genus and Species become NA
  N <- 20
  # We only want the 20 most abundant families in the RNA phyloseq object
  top <- names(sort(taxa_sums(phy.family), decreasing = TRUE))[1:N]
  # Creates a vector with the 20 most abundant families in decreasing order from the RNA phyloseq object!!!
  phy.family.relative.top <- prune_taxa(top, phy.family)
  # subsets the pyhloseq object with relative abundances by the 20 most abundant families
  plot_bar(phy.family.relative.top, x = "ID", fill = "Family") + 
    facet_grid(LT_temp~ CBASS_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ylab("Relative abundance") +
    scale_fill_manual(values = palette_family) +
    ggtitle(substitute(filter2))
  
  # Creates a barplot of the 20 most abundant families with relative abundance
  # Some of the families are multiple times in the top 20 in the plot they are put together. 
  # For that reason only 9 families are shown.
}

find.top.family.CBASS(rare_phy, type == "RNA", Timepoint == "T1")




find.top.taxa.CBASS_30 <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.type <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  phy.type.relative <- transform_sample_counts(phy.type, function(x) x / sum(x) )
  phy.family <- tax_glom(phy.type.relative, taxrank = "Family")
  N <- 20
  top <- names(sort(taxa_sums(phy.family), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy.family)
  plot_bar(phy.family.relative.top,x = "ID", fill = "Family") + 
    facet_grid(LT_temp ~Timepoint, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ylab("Relative abundance") +
    scale_fill_manual(values = palette_family) +
    ggtitle(substitute(filter1),substitute(filter2))
}
find.top.taxa.CBASS_30(rare_phy, type == "DNA")


########################################################################
### General top 20 families over genotypes
########################################################################




# Concatenate unique values in a vector
concat_unique <- function(vec){
  uniq <- unique(as.character(vec))
  return(paste(uniq, collapse = "/"))
}

# Like psmelt, but only uses the otu_table and sample_data
ps_semi_melt <- function(ps){
  otu_table(ps) %>%
    data.frame(taxid = row.names(.)) %>%
    rename_with(function(x){gsub("X", "", x)}) %>%
    pivot_longer(!taxid, names_to = "sample_id", values_to = "abundance") %>%
    left_join(sample_data(ps) %>%
                data.frame(sample_id = row.names(.)),
              by = "sample_id")
}

# Function that summarizes a vector based on its class
summarise_vec <- function(vec){
  if(class(vec) %in% c("numeric", "integer", "logical")){
    return(mean(vec, na.rm = T))
  } else if (class(vec) %in% c("factor", "character")){
    return(concat_unique(vec))
  } else {
    stop("Error: unknown column type")
  }
}

# Converts a summary df to an otu_table
summ_to_otu_tbl <- function(summ){
  summ %>% 
    select(taxid, sample_id, abundance) %>% 
    pivot_wider(names_from = "sample_id", values_from = "abundance") %>%
    column_to_rownames('taxid') %>%
    as.matrix() %>%
    otu_table(, taxa_are_rows = TRUE)
}

# Converts a summary df to sample_data
summ_to_sample_dat <- function(summ){
  summ %>% 
    select(!c(taxid, abundance)) %>% 
    unique() %>%
    column_to_rownames('sample_id') %>%
    sample_data()
}

# Function that merges phyloseq samples based on the names of one or more grouping factors
# present in sample_data(ps)
merge_ps_samples <- function(ps, grouping){
  
  # Make sure taxa are rows
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_table(ps) <- phyloseq::otu_table(t(otu_table(ps)), taxa_are_rows = T)
  }
  
  # Convert to long format
  ps_long <- ps_semi_melt(ps)
  
  # Summarise all columns
  summ <- ps_long %>%
    group_by(across(all_of(!!grouping))) %>%
    group_by(taxid, .add = T) %>%
    summarise(across(everything(), summarise_vec)) %>%
    ungroup()
  
  # Convert to otu_table and sample_data
  otu_tbl <- summ_to_otu_tbl(summ)
  sample_dat <- summ_to_sample_dat(summ)
  
  # Create new physeq object
  new_ps <- phyloseq(otu_tbl, sample_dat, tax_table(ps))
  return(new_ps)
}


find.top.family_overall <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type_T <- do.call("subset_samples", list(quote(phy_type), substitute(filter2)))
  phy_type_T_r <- transform_sample_counts(phy_type_T, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  phy_family <- tax_glom(general_phy, taxrank = "Family")
  N <- 20
  top <- names(sort(taxa_sums(phy_family), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy_family)
  plot_bar(phy.family.relative.top, fill = "Family") +
    facet_grid(~ LT_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_family) +
    scale_x_discrete(labels = phy.family.relative.top@sam_data$CBASS_temp) +
  theme(plot.title = element_text(size = 20, face="bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15))+
    ggtitle(substitute(filter2)) 
}
  
find.top.family_overall(rare_phy, type == "RNA", Timepoint == "T3", group = c("LT_temp", "CBASS_temp"))


find.top.family_overall_2 <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type <- do.call("subset_samples", list(quote(phy_type), substitute(filter1)))
  phy_type_T_r <- transform_sample_counts(phy_type, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  phy_family <- tax_glom(general_phy, taxrank = "Family")
  N <- 20
  top <- names(sort(taxa_sums(phy_family), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy_family)
  plot_bar(phy.family.relative.top, x = "Timepoint", fill = "Family") +
    facet_grid(LT_temp ~ CBASS_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_family) +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          legend.margin = margin(t = -10))+
    ggtitle(substitute(filter1))
  
}

find.top.family_overall_2(rare_phy, type == "RNA", group = c("Timepoint","LT_temp", "CBASS_temp"))


find.top.family_overall_CBASS30 <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type <- do.call("subset_samples", list(quote(phy_type), substitute(filter1)))
  phy_type_T_r <- transform_sample_counts(phy_type, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  phy_family <- tax_glom(general_phy, taxrank = "Family")
  N <- 20
  top <- names(sort(taxa_sums(phy_family), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy_family)
  plot_bar(phy.family.relative.top, x = "Timepoint", fill = "Family") +
    facet_grid(~ LT_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_family) +
    theme(plot.title = element_text(size = 20, face="bold"),
          plot.subtitle = element_text(size = 15, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          plot.background = element_rect(colour = "black"),
          legend.margin = margin(t = -10))+
    ggtitle(substitute(filter1),substitute(filter2))
  
}
top20_DNA_CBASS30 <- find.top.family_overall_CBASS30(rare_phy, type == "DNA", CBASS_temp == 30, group = c("Timepoint", "LT_temp"))
top20_RNA_CBASS30 <- find.top.family_overall_CBASS30(rare_phy, type == "RNA", CBASS_temp == 30, group = c("Timepoint", "LT_temp"))



nested_top20 <- (top20_DNA_CBASS30/top20_RNA_CBASS30)+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested_top20

########################################################################
### function to easy filter and then plot the 20 most abundant ASVs
########################################################################


find.top.ASV <- function(data, filter1, filter2){ 
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.type.time <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  phy.type.relative <- transform_sample_counts(phy.type.time, function(x) x / sum(x) )
  N <- 20
  top <- names(sort(taxa_sums(phy.type.relative), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy.type.relative)
  phy.family.relative.top <- tax_names2rank(phy.family.relative.top, colname = "ASV")
  plot_bar(phy.family.relative.top,x = "ID", fill = "ASV") + 
    facet_grid(LT_temp~ CBASS_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ylab("Relative abundance") +
    scale_fill_manual(values = palette_ASV) +
    ggtitle(substitute(filter2))
}

find.top.ASV(rare_phy, type == "RNA", Timepoint == "T3")


find.top.ASV.CBASS_30 <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.type <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  phy.type.relative <- transform_sample_counts(phy.type, function(x) x / sum(x) )
  N <- 20
  top <- names(sort(taxa_sums(phy.type.relative), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, phy.type.relative)
  phy.family.relative.top <- tax_names2rank(phy.family.relative.top, colname = "ASV")
  plot_bar(phy.family.relative.top,x = "ID", fill = "ASV") + 
    facet_grid(LT_temp ~Timepoint, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ylab("Relative abundance") +
    scale_fill_manual(values = palette_ASV) +
    ggtitle(substitute(filter1),substitute(filter2))
}
find.top.ASV.CBASS_30(rare_phy, type == "RNA", CBASS_temp == 30)


########################################################################
### General top 20 ASVs over genotypes
########################################################################

find.top.ASV_overall <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type_T <- do.call("subset_samples", list(quote(phy_type), substitute(filter2)))
  phy_type_T_r <- transform_sample_counts(phy_type_T, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  N <- 20
  top <- names(sort(taxa_sums(general_phy), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, general_phy)
  phy.family.relative.top <- tax_names2rank(phy.family.relative.top, colname = "ASV")
  plot_bar(phy.family.relative.top, fill = "ASV") +
    facet_grid(~ LT_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_ASV) +
    scale_x_discrete(labels = phy.family.relative.top@sam_data$CBASS_temp) +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ggtitle(substitute(filter2))
}

find.top.ASV_overall(rare_phy, type == "RNA", CBASS_temp == "30", group = c("LT_temp", "CBASS_temp"))


find.top.ASV_overall_2 <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type <- do.call("subset_samples", list(quote(phy_type), substitute(filter1)))
  phy_type_T_r <- transform_sample_counts(phy_type, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  N <- 20
  top <- names(sort(taxa_sums(general_phy), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, general_phy)
  phy.family.relative.top <- tax_names2rank(phy.family.relative.top, colname = "ASV")
  plot_bar(phy.family.relative.top, x = "Timepoint", fill = "ASV") +
    facet_grid(LT_temp ~ CBASS_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_ASV) +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15)) +
    ggtitle(substitute(filter1))
}

find.top.ASV_overall_2(rare_phy, type == "RNA", group = c("Timepoint", "LT_temp", "CBASS_temp"))



find.top.ASV_overall_CBASS30 <- function(data, filter1, filter2, group){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_type <- do.call("subset_samples", list(quote(phy_type), substitute(filter1)))
  phy_type_T_r <- transform_sample_counts(phy_type, function(x) x / sum(x) )
  general_phy <- merge_ps_samples(phy_type_T_r, group)
  otu_table(general_phy) <- otu_table(t(otu_table(general_phy)), taxa_are_rows = F)
  N <- 20
  top <- names(sort(taxa_sums(general_phy), decreasing = TRUE))[1:N]
  phy.family.relative.top <- prune_taxa(top, general_phy)
  phy.family.relative.top <- tax_names2rank(phy.family.relative.top, colname = "ASV")
  plot_bar(phy.family.relative.top, x = "Timepoint", fill = "ASV") +
    facet_grid(~ LT_temp, scales = "free") +
    theme_bw() +
    ylim(0, 1) +
    scale_fill_manual(values = palette_ASV) +
    theme(plot.title = element_text(size = 20, face="bold"),
          plot.subtitle = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          plot.background = element_rect(colour = "black"),
          legend.margin = margin(t = -10)) +
    ggtitle(substitute(filter1),substitute(filter2))
}



top20_DNA_CBASS30_ASV <- find.top.ASV_overall_CBASS30(rare_phy, type == "DNA", CBASS_temp == 30, group = c("Timepoint", "LT_temp"))
top20_RNA_CBASS30_ASV <- find.top.ASV_overall_CBASS30(rare_phy, type == "RNA", CBASS_temp == 30, group = c("Timepoint", "LT_temp"))



nested_top20_ASV <- (top20_DNA_CBASS30_ASV/top20_RNA_CBASS30_ASV)+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested_top20_ASV
########################################################################
### building a function to easy filter and then plot the alpha diversity
########################################################################


######### Shannon (Diversity Index)

find.alpha.diversity.CBASS <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for LT_temp
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  anno_df <- compare_means(Shannon ~ CBASS_temp, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(CBASS_temp = group1,
           ypos = rep(c(5.5, 5.3, 5.1, 4.9, 4.7, 4.5), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.05, "*", NA)
  anno_df_filter <- filter(anno_df, p < 0.05)
  plot <- plot_richness(phy.temp, x = "CBASS_temp", measures = c("Shannon"), color = "CBASS_temp") +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, stackratio = 1, binwidth = 0.1, 
                 aes(fill = CBASS_temp),
                 position=position_dodge(0.75)) +
    facet_wrap(~ LT_temp) +
    ylim(2.5, 5.5) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = group1, 
                    xmax = group2, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.background = element_rect(colour = "black")) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(substitute(filter2))
  plot$layers <- plot$layers[-1]
  # plot the alpha diversity based on Shannon diversity
  out <- list(plot,anno_df)
  return(out)
}


find.alpha.diversity.CBASS(rare_phy, type == "RNA", Timepoint == "T3")



find.alpha.diversity.noCBASS <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for RNA or DNA
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  anno_df <- compare_means(Shannon ~ Timepoint, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(Timepoint = group1,
           ypos = rep(c(5.5, 5.2, 5.0), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.05, "*", NA)
  anno_df_filter <- filter(anno_df, p < 0.05)
  plot <- plot_richness(phy.temp, x = "Timepoint", measures = c("Shannon"), color = "Timepoint") +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, stackratio = 1, binwidth = 0.1,
                 aes(fill = Timepoint),
                 position=position_dodge(0.75)) +
    facet_wrap(~ LT_temp) +
    ylim(2.5, 5.5) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = group1, 
                    xmax = group2, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          plot.subtitle = element_text(size = 15, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.background = element_rect(colour = "black")) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(substitute(filter1),substitute(filter2))
  plot$layers <- plot$layers[-1]
  # plot the alpha diversity based on Shannon diversity
  out <- list(plot,anno_df)
  return(out)
  
}


RNA_alpha_CBASS30 <- find.alpha.diversity.noCBASS(rare_phy, type == "RNA", CBASS_temp == 30)
DNA_alpha_CBASS30 <- find.alpha.diversity.noCBASS(rare_phy, type == "DNA", CBASS_temp == 30)
nested_alpha <- (DNA_alpha_CBASS30[[1]]/RNA_alpha_CBASS30[[1]])+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested_alpha



find.alpha.diversity.noCBASS_Timepoint <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for RNA or DNA
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  anno_df <- compare_means(Shannon ~ LT_temp, group.by = "Timepoint", data = together_df, p.adjust.method = "bonf") |> 
    mutate(LT_temp = group1,
           ypos = rep(c(5.5, 5.2, 5.0), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.05, "*", NA)
  anno_df_filter <- filter(anno_df, p < 0.05)
  plot <- plot_richness(phy.temp, x = "LT_temp", measures = c("Shannon"), color = "LT_temp") +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.5, stackratio = 1, binwidth = 0.1,
                 aes(fill = LT_temp),
                 position=position_dodge(0.75)) +
    facet_wrap(~ Timepoint) +
    ylim(2.5, 5.5) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = group1, 
                    xmax = group2, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(substitute(filter1),substitute(filter2))
  plot$layers <- plot$layers[-1]
  # plot the alpha diversity based on Shannon diversity
  out <- list(plot,anno_df)
  return(out)
  
}





######### Chao1 (Richness)

find.chao <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for LT_temp
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  anno_df <- compare_means(Chao1 ~ CBASS_temp, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(CBASS_temp = group1,
           ypos = rep(c(420, 405, 390, 375, 360, 345), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.5, "*", NA)
  anno_df_filter <- filter(anno_df, p < 0.05)
  plot <- plot_richness(phy.temp, x = "CBASS_temp", measures = c("Chao1"), color = "CBASS_temp") +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 7, stackratio = 1, binwidth = 1, 
                 aes(fill = CBASS_temp),
                 position=position_dodge(0.75)) +
    facet_wrap(~ LT_temp) +
    ylim(75, 420) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = group1, 
                    xmax = group2, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(substitute(filter2))
  plot$layers <- plot$layers[-c(1:2)]
  # plot the alpha diversity based on Shannon diversity
  out <- list(plot,anno_df)
  return(out)
  
}

find.chao(rare_phy, type == "RNA", Timepoint == "T1")

find.chao2 <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for LT_temp
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  anno_df <- compare_means(Chao1 ~ Timepoint, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(Timepoint = group1,
           ypos = rep(c(420, 390, 360), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.5, "*", NA)
  anno_df_filter <- filter(anno_df, p < 0.05)
  plot <- plot_richness(phy.temp, x = "Timepoint", measures = c("Chao1"), color = "Timepoint") +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 7, stackratio = 1, binwidth = 1, 
                 aes(fill = Timepoint),
                 position=position_dodge(0.75)) +
    facet_wrap(~ LT_temp) +
    ylim(75, 420) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = group1, 
                    xmax = group2, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ggtitle(substitute(filter1),substitute(filter2))
  plot$layers <- plot$layers[-c(1:2)]
  # plot the alpha diversity based on Shannon diversity
  out <- list(plot,anno_df)
  return(out)
  
}

find.chao2(rare_phy, type == "DNA")
######### Observed (Richness)

find.observed <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  test_1 <- together_df %>%
    group_by(LT_temp,CBASS_temp) %>%
    arrange(LT_temp,CBASS_temp) %>%
    filter(row_number()==1) |> 
    mutate(group1 = CBASS_temp,
           ID_start = ID) |> 
    ungroup()
  test_2 <- together_df %>%
    group_by(LT_temp,CBASS_temp) %>%
    arrange(LT_temp,CBASS_temp) %>%
    filter(row_number()==n())|> 
    mutate(group2 = CBASS_temp,
           ID_end = ID) |> 
    ungroup()
  anno_df <- compare_means(Observed ~ CBASS_temp, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(CBASS_temp = group1,
           ypos = rep(c(420, 400, 380, 360, 340, 320), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.5, "*", NA)
  zusammen <- left_join(anno_df, test_1[, c("LT_temp", "group1", "ID_start")], by = c("LT_temp", "group1"))
  zusammen <- left_join(zusammen, test_2[, c("LT_temp", "group2", "ID_end")], by = c("LT_temp", "group2"))
  anno_df_filter <- filter(zusammen, p < 0.05)
  # subset for LT_temp
  plot <- plot_richness(phy.temp, x = "ID", measures = c("Observed"), color = "CBASS_temp") +
    geom_bar(aes(fill = CBASS_temp), color = "black",
             stat='identity') +
    facet_grid(cols = vars(LT_temp), scale = "free") +
    ylim(0, 420) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = ID_start, 
                    xmax = ID_end, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    ggtitle(substitute(filter2))
  # plot the alpha diversity based on Shannon diversity
  plot$layers <- plot$layers[-1]
  out <- list(plot, anno_df, test_1, zusammen)
  return(out)
  
}
find.observed(rare_phy, type == "RNA", Timepoint == "T3")


find.observed2 <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for LT_temp
  es_t1 <- estimate_richness(phy.temp)
  test_df <- as.data.frame(sample_data(phy.temp))
  together_df <- cbind(es_t1, test_df)
  test_1 <- together_df %>%
    group_by(LT_temp,Timepoint) %>%
    arrange(LT_temp,Timepoint) %>%
    filter(row_number()==1) |> 
    mutate(group1 = Timepoint,
           ID_start = ID) |> 
    ungroup()
  test_2 <- together_df %>%
    group_by(LT_temp,Timepoint) %>%
    arrange(LT_temp,Timepoint) %>%
    filter(row_number()==n())|> 
    mutate(group2 = Timepoint,
           ID_end = ID) |> 
    ungroup()
  anno_df <- compare_means(Observed ~ Timepoint, group.by = "LT_temp", data = together_df, p.adjust.method = "bonf") |> 
    mutate(Timepoint = group1,
           ypos = rep(c(400, 380, 360), times = 3))
  anno_df$p.adjust.signif <- ifelse(anno_df$p.adj <= 0.5, "*", NA)
  zusammen <- left_join(anno_df, test_1[, c("LT_temp", "group1", "ID_start")], by = c("LT_temp", "group1"))
  zusammen <- left_join(zusammen, test_2[, c("LT_temp", "group2", "ID_end")], by = c("LT_temp", "group2"))
  anno_df_filter <- filter(zusammen, p < 0.05)
  plot <- plot_richness(phy.temp, x = "ID", measures = c("Observed"), color = "Timepoint") +
    geom_bar(aes(fill = Timepoint), color = "black",
             stat='identity') +
    facet_grid(cols = vars(LT_temp), scale = "free") +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 420) +
    geom_signif(data=anno_df_filter, 
                aes(xmin = ID_start, 
                    xmax = ID_end, 
                    annotations = p.adjust.signif,
                    y_position = ypos), 
                manual= TRUE,
                textsize = 7,
                show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    ggtitle(substitute(filter1),substitute(filter2))
  # plot the alpha diversity based on Shannon diversity
  plot$layers <- plot$layers[-1]
  out <- list(plot, anno_df, zusammen)
  return(out)
  
}

find.observed2(rare_phy, type == "RNA", CBASS_temp == 30)






########################################################################
### building a function to easy filter and then do a PCoA plot
########################################################################

pcoa.plot.CBASS <- function(data, filter1, filter2){
  pcoa.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA
  pcoa.temp <- do.call("subset_samples", list(quote(pcoa.type), substitute(filter2)))
  # subset for LT_temp
  pcoa.ord <- ordinate(pcoa.temp, method = "PCoA", distance = "bray")
  # calculate the distances, based on bray-curtis
  
  plot <- plot_ordination(pcoa.temp,pcoa.ord, color = "CBASS_temp", shape = "LT_temp")  + 
    geom_point(size = 5, alpha = 1) +
    ggtitle(substitute(filter2)) +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")# +
   # stat_ellipse(aes(color = CBASS_temp, group = CBASS_temp),type = "t")
  #PCoA plot 
  plot
  
}

pcoa.plot.CBASS(rare_phy, type == "RNA", Timepoint == "T2" )


pcoa.plot.2 <- function(data, filter1, filter2){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  # subset for RNA or DNA (or what ever you like)
  phy.temp <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  # subset for RNA or DNA
  pcoa.ord <- ordinate(phy.temp, method = "PCoA", distance = "bray")
  # calculate the distances, based on bray-curtis
  
  plot <- plot_ordination(phy.temp,pcoa.ord, shape = "LT_temp", color = "Timepoint")  + 
    geom_point(size = 8, alpha = 1) +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 20, face="bold"),
          plot.subtitle = element_text(size = 15, face="bold"),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.background = element_rect(colour = "black")) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")+
    ggtitle(substitute(filter1),substitute(filter2))
  #PCoA plot 
  plot
  
}

RNA_CBASS30_pcoa <- pcoa.plot.2(rare_phy, type == "RNA", CBASS_temp == 30)
DNA_CBASS30_pcoa <- pcoa.plot.2(rare_phy, type == "DNA", CBASS_temp == 30)

nested_pcoa <- (DNA_CBASS30_pcoa|RNA_CBASS30_pcoa)+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested_pcoa









########################################################################
### DESeq2 and Differentially abundant ASVs (DAAs)
########################################################################



DAA_CBASS_temp <- function(data, filter1, filter2, filter3, filter4, filter5, filter6, filter7){
  phy.type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy.lt <- do.call("subset_samples", list(quote(phy.type), substitute(filter2)))
  phy.time <- do.call("subset_samples", list(quote(phy.lt), substitute(filter3)))
  phy.temp_30 <- do.call("subset_samples", list(quote(phy.time), substitute(filter4)))
  phy.temp_33 <- do.call("subset_samples", list(quote(phy.time), substitute(filter5)))
  phy.temp_33 <- merge_phyloseq(phy.temp_30, phy.temp_33)
  des_33 <- phyloseq_to_deseq2(phy.temp_33, ~ CBASS_temp)
  des_33 <-  DESeq(des_33, test="Wald", fitType="parametric")
  res_33 <-  results(des_33)
  res_33 <-  res_33[order(res_33$padj, na.last=NA), ]
  alpha <-  0.05
  df_33 <-  res_33[(res_33$padj < alpha), ]
  df_33 <-  cbind(as(df_33, "data.frame"), as(tax_table(phy.temp_33)[rownames(df_33), ], "matrix"))
  x_33 <- sum(df_33$log2FoldChange > 0, na.rm=TRUE)
  y_33 <- sum(df_33$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_33 <- "33 vs 30"
  phy.temp_36 <- do.call("subset_samples", list(quote(phy.time), substitute(filter6)))
  phy.temp_36 <- merge_phyloseq(phy.temp_30, phy.temp_36)
  des_36 <- phyloseq_to_deseq2(phy.temp_36, ~ CBASS_temp)
  des_36 <-  DESeq(des_36, test="Wald", fitType="parametric")
  res_36 <-  results(des_36)
  res_36 <-  res_36[order(res_36$padj, na.last=NA), ]
  alpha <-  0.05
  df_36 <-  res_36[(res_36$padj < alpha), ]
  df_36 <- cbind(as(df_36, "data.frame"), as(tax_table(phy.temp_36)[rownames(df_36), ], "matrix"))
  x_36 <- sum(df_36$log2FoldChange > 0, na.rm=TRUE)
  y_36 <- sum(df_36$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_36 <- "36 vs 30"
  phy.temp_39 <- do.call("subset_samples", list(quote(phy.time), substitute(filter7)))
  phy.temp_39 <- merge_phyloseq(phy.temp_30, phy.temp_39)
  des_39 <- phyloseq_to_deseq2(phy.temp_39, ~ CBASS_temp)
  des_39 <-  DESeq(des_39, test="Wald", fitType="parametric")
  res_39 <-  results(des_39)
  res_39 <-  res_39[order(res_39$padj, na.last=NA), ]
  alpha <-  0.05
  df_39 <-  res_39[(res_39$padj < alpha), ]
  df_39 <-  cbind(as(df_39, "data.frame"), as(tax_table(phy.temp_39)[rownames(df_39), ], "matrix"))
  x_39 <- sum(df_39$log2FoldChange > 0, na.rm=TRUE)
  y_39 <- sum(df_39$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_39 <- "39 vs 30"
  ext_t <- as.data.frame(sample_data(phy.temp_30))
  ext_ti <- unique(ext_t$Timepoint)
  ext_LT <- unique(ext_t$LT_temp)
  mal_schauen <- data.frame(DAAs = c(x_33, y_33, x_36, y_36, x_39, y_39),
                            Comparison = rep(c(comp_33, comp_36, comp_39), each = 2),
                            Timepoint = ext_ti,
                            LT_temp = ext_LT)
  in_33 <- filter(df_33, log2FoldChange > 0)
  de_33 <- filter(df_33, log2FoldChange < 0)
  in_36 <- filter(df_36, log2FoldChange > 0)
  de_36 <- filter(df_36, log2FoldChange < 0)
  in_39 <- filter(df_39, log2FoldChange > 0)
  de_39 <- filter(df_39, log2FoldChange < 0)
  list_data <- list(mal_schauen, 
                    df_33, in_33, de_33, des_33,
                    df_36, in_36, de_36, des_36,
                    df_39, in_39, de_39, des_39)
  names(list_data) <- c("Log2FC", 
                        "CBASS_33_vs_30", "increase_CBASS_33_vs_30", "decrease_CBASS_33_vs_30", "dds_33",
                        "CBASS_36_vs_30", "increase_CBASS_36_vs_30", "decrease_CBASS_36_vs_30", "dds_36",
                        "CBASS_39_vs_30", "increase_CBASS_39_vs_30", "decrease_CBASS_39_vs_30", "dds_39")
  list_data
}



DAA_timepoint <- function(data, filter1, filter2, filter3, filter4, filter5, filter6){
  phy_type <- do.call("subset_samples", list(quote(data), substitute(filter1)))
  phy_lt <- do.call("subset_samples", list(quote(phy_type), substitute(filter2)))
  phy_CBASS <- do.call("subset_samples", list(quote(phy_lt), substitute(filter3)))
  phy_T1 <- do.call("subset_samples", list(quote(phy_CBASS), substitute(filter4)))
  phy_T2 <- do.call("subset_samples", list(quote(phy_CBASS), substitute(filter5)))
  phy_T2 <- merge_phyloseq(phy_T1, phy_T2)
  des_T2 <- phyloseq_to_deseq2(phy_T2, ~ Timepoint)
  des_T2 <-  DESeq(des_T2, test="Wald", fitType="parametric")
  res_T2 <-  results(des_T2)
  res_T2 <-  res_T2[order(res_T2$padj, na.last=NA), ]
  alpha <-  0.05
  df_T2 <-  res_T2[(res_T2$padj < alpha), ]
  df_T2 <-  cbind(as(df_T2, "data.frame"), as(tax_table(phy_T2)[rownames(df_T2), ], "matrix"))
  x_T2 <- sum(df_T2$log2FoldChange > 0, na.rm=TRUE)
  y_T2 <- sum(df_T2$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_T2 <- "T2 vs T1"
  phy_T3 <- do.call("subset_samples", list(quote(phy_CBASS), substitute(filter6)))
  phy_T3 <- merge_phyloseq(phy_T1, phy_T3)
  des_T3 <- phyloseq_to_deseq2(phy_T3, ~ Timepoint)
  des_T3 <-  DESeq(des_T3, test="Wald", fitType="parametric")
  res_T3 <-  results(des_T3)
  res_T3 <-  res_T3[order(res_T3$padj, na.last=NA), ]
  alpha <-  0.05
  df_T3 <-  res_T3[(res_T3$padj < alpha), ]
  df_T3 <- cbind(as(df_T3, "data.frame"), as(tax_table(phy_T3)[rownames(df_T3), ], "matrix"))
  x_T3 <- sum(df_T3$log2FoldChange > 0, na.rm=TRUE)
  y_T3 <- sum(df_T3$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_T3 <- "T3 vs T1"
  phy_T3_T2 <- merge_phyloseq(phy_T2, phy_T3)
  des_T3_T2 <- phyloseq_to_deseq2(phy_T3_T2, ~ Timepoint)
  des_T3_T2 <-  DESeq(des_T3_T2, test="Wald", fitType="parametric")
  res_T3_T2 <-  results(des_T3_T2)
  res_T3_T2 <-  res_T3_T2[order(res_T3_T2$padj, na.last=NA), ]
  alpha <-  0.05
  df_T3_T2 <-  res_T3_T2[(res_T3_T2$padj < alpha), ]
  df_T3_T2 <-  cbind(as(df_T3_T2, "data.frame"), as(tax_table(phy_T3_T2)[rownames(df_T3_T2), ], "matrix"))
  x_T3_T2 <- sum(df_T3_T2$log2FoldChange > 0, na.rm=TRUE)
  y_T3_T2 <- sum(df_T3_T2$log2FoldChange < 0, na.rm=TRUE) * -1
  comp_T3_T2 <- "T3 vs T2"
  ext_t <- as.data.frame(sample_data(phy_T1))
  ext_CBASS <- unique(ext_t$CBASS_temp)
  ext_LT <- unique(ext_t$LT_temp)
  mal_schauen <- data.frame(DAAs = c(x_T2, y_T2, x_T3, y_T3, x_T3_T2, y_T3_T2),
                            Comparison = rep(c(comp_T2, comp_T3, comp_T3_T2), each = 2),
                            CBASS_temp = ext_CBASS,
                            LT_temp = ext_LT)
  in_T2 <- filter(df_T2, log2FoldChange > 0)
  de_T2 <- filter(df_T2, log2FoldChange < 0)
  in_T3 <- filter(df_T3, log2FoldChange > 0)
  de_T3 <- filter(df_T3, log2FoldChange < 0)
  in_T3_T2 <- filter(df_T3_T2, log2FoldChange > 0)
  de_T3_T2 <- filter(df_T3_T2, log2FoldChange < 0)
  list_data <- list(mal_schauen, 
                    df_T2, in_T2, de_T2, des_T2,
                    df_T3, in_T3, de_T3, des_T3,
                    df_T3_T2, in_T3_T2, de_T3_T2, des_T3_T2)
  names(list_data) <- c("Log2FC", 
                        "Timepoint_T2_vs_T1", "increase_Timepoint_T2_vs_T1", "decrease_Timepoint_T2_vs_T1", "dds_T2",
                        "Timepoint_T3_vs_T1", "increase_Timepoint_T3_vs_T1", "decrease_Timepoint_T3_vs_T1", "dds_T3",
                        "Timepoint_T3_vs_T2", "increase_Timepoint_T3_vs_T2", "decrease_Timepoint_T3_vs_T2", "dds_T3_T2")
  list_data
}

DNA_LT_25_CBASS_30 <- DAA_timepoint(phy, 
             type == "DNA", 
             LT_temp == 25,
             CBASS_temp == 30, 
             Timepoint == "T1", 
             Timepoint == "T2", 
             Timepoint == "T3")
DNA_LT_27_CBASS_30 <- DAA_timepoint(phy, 
                                type == "DNA", 
                                LT_temp == 27,
                                CBASS_temp == 30, 
                                Timepoint == "T1", 
                                Timepoint == "T2", 
                                Timepoint == "T3")
DNA_LT_30_CBASS_30 <- DAA_timepoint(phy, 
                                type == "DNA", 
                                LT_temp == 30,
                                CBASS_temp == 30, 
                                Timepoint == "T1", 
                                Timepoint == "T2", 
                                Timepoint == "T3")


DNA_together_CBASS30 <- bind_rows(DNA_LT_25_CBASS_30$Log2FC, 
                              DNA_LT_27_CBASS_30$Log2FC, 
                              DNA_LT_30_CBASS_30$Log2FC)

DNA_neg_CBASS30_LT25 <- filter(DNA_LT_25_CBASS_30$Log2FC, DAAs < 0)
DNA_pos_CBASS30_LT25 <- filter(DNA_LT_25_CBASS_30$Log2FC, DAAs >= 0)


DNA_CBASS30_LT25_T1vs <- ggplot(DNA_LT_25_CBASS_30$Log2FC, aes(x = Comparison, y = DAAs, fill = Comparison))+
  geom_bar(stat = "identity")+
  geom_text(DNA_neg_CBASS30_LT25, mapping = aes(label = DAAs), nudge_y = -20, size =7)+
  geom_text(DNA_pos_CBASS30_LT25, mapping = aes(label = DAAs), nudge_y = 20, size =7)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-175, 175), breaks=c(-150, -100, -50, 0, 50, 100, 150)) +
  theme_clean() +
  theme(plot.title = element_text(size = 20, face="bold"),
        plot.subtitle = element_text(size = 15, face="bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.background = element_rect(colour = "black")) +
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("DNA", "CBASS 30, LT = 25")



t1_25 <- DAA_CBASS_temp(phy, 
    type == "RNA", 
    LT_temp == 25,
    Timepoint == "T1",
    CBASS_temp == 30, 
    CBASS_temp == 33, 
    CBASS_temp == 36, 
    CBASS_temp == 39)
t2_25 <- DAA_CBASS_temp(phy, 
          type == "RNA", 
          LT_temp == 25,
          Timepoint == "T2",
          CBASS_temp == 30, 
          CBASS_temp == 33, 
          CBASS_temp == 36, 
          CBASS_temp == 39)
t3_25 <- DAA_CBASS_temp(phy, 
          type == "RNA", 
          LT_temp == 25,
          Timepoint == "T3",
          CBASS_temp == 30, 
          CBASS_temp == 33, 
          CBASS_temp == 36, 
          CBASS_temp == 39)
t1_27 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 27,
             Timepoint == "T1",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)
t2_27 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 27,
             Timepoint == "T2",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)
t3_27 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 27,
             Timepoint == "T3",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)
t1_30 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 30,
             Timepoint == "T1",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)
t2_30 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 30,
             Timepoint == "T2",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)
t3_30 <- DAA_CBASS_temp(phy, 
             type == "RNA", 
             LT_temp == 30,
             Timepoint == "T3",
             CBASS_temp == 30, 
             CBASS_temp == 33, 
             CBASS_temp == 36, 
             CBASS_temp == 39)


allt <- bind_rows(t1_25$Log2FC, t1_27$Log2FC ,t1_30$Log2FC,
          t2_25$Log2FC ,t2_27$Log2FC ,t2_30$Log2FC,
          t3_25$Log2FC ,t3_27$Log2FC ,t3_30$Log2FC)

neg_allt <- filter(allt, DAAs < 0)
pos_allt <- filter(allt, DAAs >= 0)

ggplot(allt, aes(x = Comparison, y = DAAs, fill = Comparison))+
  geom_bar(stat = "identity")+
  geom_text(neg_allt, mapping = aes(label = DAAs), nudge_y = -20)+
  geom_text(pos_allt, mapping = aes(label = DAAs), nudge_y = 20)+
  geom_hline(yintercept = 0)+
  facet_grid(LT_temp ~ Timepoint) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face="bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))



######################### CBASS_temp 36 vs 30



allt_30vs36 <- filter(allt, Comparison == "36 vs 30")
neg_allt_30vs36 <- filter(allt_30vs36, DAAs < 0)
pos_allt_30vs36 <- filter(allt_30vs36, DAAs >= 0)

all_LT_Time <- ggplot(allt_30vs36, aes(x = LT_temp, y = DAAs, fill = Timepoint))+
  geom_bar(stat = "identity")+
  geom_text(neg_allt_30vs36, mapping = aes(label = DAAs), nudge_y = -20, size = 9)+
  geom_text(pos_allt_30vs36, mapping = aes(label = DAAs), nudge_y = 20, size = 9)+
  geom_hline(yintercept = 0)+
  facet_wrap(~ Timepoint) +
  theme_clean() +
  theme(plot.title = element_text(size = 35, face="bold"),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        strip.text = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30)) +
  scale_fill_brewer(palette = "Dark2") 
all_LT_Time


# Identify log fold change of ASVs
### 
L2FC_t1_25 <- t1_25$CBASS_36_vs_30 |> mutate(ASVs = row.names(t1_25$CBASS_36_vs_30 ),
                                             LT = 25) 
row.names(L2FC_t1_25) <- NULL
L2FC_t1_27 <- t1_27$CBASS_36_vs_30 |> mutate(ASVs = row.names(t1_27$CBASS_36_vs_30 ), 
                                             LT = 27) 
row.names(L2FC_t1_27) <- NULL
L2FC_t1_30 <- t1_30$CBASS_36_vs_30 |> mutate(ASVs = row.names(t1_30$CBASS_36_vs_30 ),
                                             LT = 30) 
row.names(L2FC_t1_30) <- NULL
log_FC_T1 <- bind_rows(L2FC_t1_25,
                       L2FC_t1_27,
                       L2FC_t1_30)
L2FC_t2_25 <- t2_25$CBASS_36_vs_30 |> mutate(ASVs = row.names(t2_25$CBASS_36_vs_30 ),
                                             LT = 25) 
row.names(L2FC_t2_25) <- NULL
L2FC_t2_27 <- t2_27$CBASS_36_vs_30 |> mutate(ASVs = row.names(t2_27$CBASS_36_vs_30 ),
                                             LT = 27) 
row.names(L2FC_t2_27) <- NULL
L2FC_t2_30 <- t2_30$CBASS_36_vs_30 |> mutate(ASVs = row.names(t2_30$CBASS_36_vs_30 ),
                                             LT = 30) 
row.names(L2FC_t2_30) <- NULL
log_FC_T2 <- bind_rows(L2FC_t2_25,
                       L2FC_t2_27,
                       L2FC_t2_30)
L2FC_t3_25 <- t3_25$CBASS_36_vs_30 |> mutate(ASVs = row.names(t3_25$CBASS_36_vs_30 ),
                                             LT = 25) 
row.names(L2FC_t3_25) <- NULL
L2FC_t3_27 <- t3_27$CBASS_36_vs_30 |> mutate(ASVs = row.names(t3_27$CBASS_36_vs_30 ),
                                             LT = 27) 
row.names(L2FC_t3_27) <- NULL
L2FC_t3_30 <- t3_30$CBASS_36_vs_30 |> mutate(ASVs = row.names(t3_30$CBASS_36_vs_30 ),
                                             LT = 30) 
row.names(L2FC_t3_30) <- NULL
log_FC_T3 <- bind_rows(L2FC_t3_25,
                       L2FC_t3_27,
                       L2FC_t3_30)
L2FC_t3_30 <- L2FC_t3_30 |> select(ASVs,LT, everything()) |> 
  arrange(ASVs)

write.xlsx(L2FC_t3_30,"output directory",colNames = TRUE)




####### Create tables to identify DAA numbers
Intersection_ASVS <- read.table("ASVs_27_and_30.txt", sep = "\t" , header = T )
tax_names <- as.data.frame(tax_table(rare_phy))
tax_names$ASVs <- row.names(tax_names)
######### T1
df_T1 <- data.frame(ASVs = Intersection_ASVS$T1.LT27_LT30[1:11])
df_T1 <- merge(df_T1, tax_names, by = "ASVs")
unique(df_T1$Family)
length(unique(df_T1$Family))
test_T1 <- sort(table(df_T1$Family, useNA = "ifany"), decreasing = T)
test_T1
sum(test_T1)

df_T1_LT27_unique <- data.frame(ASVs = Intersection_ASVS$T1.LT27_unique[1:70])
df_T1_LT27_unique <- merge(df_T1_LT27_unique, tax_names, by = "ASVs")
unique(df_T1_LT27_unique$Family)
length(unique(df_T1_LT27_unique$Family))
T1_LT27_unique <- sort(table(df_T1_LT27_unique$Family, useNA = "ifany"), decreasing = T)
T1_LT27_unique
sum(T1_LT27_unique)

df_T1_LT30_unique <- data.frame(ASVs = Intersection_ASVS$T1.LT30_unique[1:70])
df_T1_LT30_unique <- merge(df_T1_LT30_unique, tax_names, by = "ASVs")
unique(df_T1_LT30_unique$Family)
length(unique(df_T1_LT30_unique$Family))
T1_LT30_unique <- sort(table(df_T1_LT30_unique$Family, useNA = "ifany"), decreasing = T)
T1_LT30_unique
sum(T1_LT30_unique)





###### T2
df_T2 <- data.frame(ASVs = Intersection_ASVS$T2.LT27_LT30[1:43])
df_T2 <- merge(df_T2, tax_names, by = "ASVs")
unique(df_T2$Family)
length(unique(df_T2$Family))
test_T2 <- sort(table(df_T2$Family, useNA = "ifany"), decreasing = F)
test_T2
sum(test_T2)

df_T2_LT27_unique <- data.frame(ASVs = Intersection_ASVS$T2.LT27_unique[1:70])
df_T2_LT27_unique <- merge(df_T2_LT27_unique, tax_names, by = "ASVs")
unique(df_T2_LT27_unique$Family)
length(unique(df_T2_LT27_unique$Family))
T2_LT27_unique <- sort(table(df_T2_LT27_unique$Family, useNA = "ifany"), decreasing = F)
T2_LT27_unique
sum(T2_LT27_unique)

df_T2_LT30_unique <- data.frame(ASVs = Intersection_ASVS$T2.LT30_unique[1:70])
df_T2_LT30_unique <- merge(df_T2_LT30_unique, tax_names, by = "ASVs")
unique(df_T2_LT30_unique$Family)
length(unique(df_T2_LT30_unique$Family))
T2_LT30_unique <- sort(table(df_T2_LT30_unique$Family, useNA = "ifany"), decreasing = F)
T2_LT30_unique
sum(T2_LT30_unique)

########## T3
df_T3 <- data.frame(ASVs = Intersection_ASVS$T3.LT27_LT30[1:107])
df_T3 <- merge(df_T3, tax_names, by = "ASVs")
unique(df_T3$Family)
length(unique(df_T3$Family))
test_T3 <- sort(table(df_T3$Family, useNA = "ifany"), decreasing = F)
test_T3
sum(test_T3)

df_T3_LT27_unique <- data.frame(ASVs = Intersection_ASVS$T3.LT27_unique[1:129])
df_T3_LT27_unique <- merge(df_T3_LT27_unique, tax_names, by = "ASVs")
unique(df_T3_LT27_unique$Family)
length(unique(df_T3_LT27_unique$Family))
T3_LT27_unique <- sort(table(df_T3_LT27_unique$Family, useNA = "ifany"), decreasing = F)
T3_LT27_unique
sum(T3_LT27_unique)

df_T3_LT30_unique <- data.frame(ASVs = Intersection_ASVS$T3.LT30_unique[1:48])
df_T3_LT30_unique <- merge(df_T3_LT30_unique, tax_names, by = "ASVs")
unique(df_T3_LT30_unique$Family)
length(unique(df_T3_LT30_unique$Family))
T3_LT30_unique <- sort(table(df_T3_LT30_unique$Family, useNA = "ifany"), decreasing = F)
T3_LT30_unique
sum(T3_LT30_unique)

filter(metadata, Timepoint == "T3" & LT_temp == 25 & CBASS_temp == 36)


######## T2vsT1_T3vsT1 LT 25 RNA
df_T2vsT1_T3vsT1 <- data.frame(ASVs = ASVs_T2vsT1_and_T3vsT1)
df_T2vsT1_T3vsT1 <- merge(df_T2vsT1_T3vsT1, tax_names, by = "ASVs")
length(unique(df_T2vsT1_T3vsT1$Family))
count(df_T2vsT1_T3vsT1$Family)
test_vector <- sort(table(df_T2vsT1_T3vsT1$Family, useNA = "ifany"), decreasing = T)
sum(test_vector)
write.xlsx(df_T1,"/Users/thomasmaisch/Desktop/Master/Intersection/T1_27_and_30.xlsx",colNames = TRUE)

######## T2vsT1_T3vsT1 LT 25 DNA
DNA_test_ASV_function
DNA_ASVs_T2vsT1_complete <- sort(DNA_test_ASV_function$T2vsT1)
DNA_ASVs_T3vsT1_complete <- sort(DNA_test_ASV_function$T3vsT1)

DNA_ASVs_T2vsT1_and_T3vsT1 <- sort(DNA_ASVs_T2vsT1_complete[DNA_ASVs_T2vsT1_complete %in% DNA_ASVs_T3vsT1_complete])
DNA_df_T2vsT1_T3vsT1 <- data.frame(ASVs = DNA_ASVs_T2vsT1_and_T3vsT1)
DNA_df_T2vsT1_T3vsT1 <- merge(DNA_df_T2vsT1_T3vsT1, tax_names, by = "ASVs")
length(unique(DNA_df_T2vsT1_T3vsT1$Family))
test_vector <- sort(table(DNA_df_T2vsT1_T3vsT1$Family, useNA = "ifany"), decreasing = T)
sum(test_vector)




#### Identify log2FC of shared 27 and 30 response

LT27_LT30_T1_l2FC <- merge(log_FC_T1, df_T1, by = c("ASVs","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
LT27_LT30_T2_l2FC <- merge(log_FC_T2, df_T2, by = c("ASVs","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
LT27_LT30_T3_l2FC <- merge(log_FC_T3, df_T3, by = c("ASVs","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

LT25_df_T1_LT27_unique


taxonomy <- as.data.frame(tax_table(rare_phy))
taxonomy <- taxonomy |> mutate(ASV = row.names(taxonomy))
taxonomy <- taxonomy |> select(ASV, everything()) |> 
  arrange(ASV)
write.xlsx(taxonomy,"your output directory",colNames = TRUE)




######Unique T1
all_DAA <- read.table("Unique_and_everywhere.txt", sep = "\t" , header = T )
Unique <- all_DAA[,c(1,2,3,5,6,7,9,10,11)]
Everywhere <- all_DAA[,c(4,8,12)]
unique_LT25_T1 <- as.data.frame(Unique$unique_LT25_T1)
unique_LT25_T1 <- as.data.frame(Unique$unique_LT25_T1) |> 
  transmute(ASVs = unique_LT25_T1[,1])
unique_LT25_T1 <- merge(unique_LT25_T1, L2FC_t1_25, by = "ASVs")


unique_LT27_T1 <- as.data.frame(Unique$unique_LT27_T1) 
unique_LT27_T1 <- as.data.frame(Unique$unique_LT27_T1) |> 
  transmute(ASVs = unique_LT27_T1[,1])
unique_LT27_T1 <- merge(unique_LT27_T1, L2FC_t1_27, by = "ASVs")

unique_LT30_T1 <- as.data.frame(Unique$unique_LT30_T1)
unique_LT30_T1 <- as.data.frame(Unique$unique_LT30_T1) |> 
  transmute(ASVs = unique_LT30_T1[,1])
unique_LT30_T1 <- merge(unique_LT30_T1, L2FC_t1_30, by = "ASVs")



####Unique T2

unique_LT25_T2 <- as.data.frame(Unique$unique_LT25_T2)
unique_LT25_T2 <- as.data.frame(Unique$unique_LT25_T2) |> 
  transmute(ASVs = unique_LT25_T2[,1])
unique_LT25_T2 <- merge(unique_LT25_T2, L2FC_t2_25, by = "ASVs")


unique_LT27_T2 <- as.data.frame(Unique$unique_LT27_T2) 
unique_LT27_T2 <- as.data.frame(Unique$unique_LT27_T2) |> 
  transmute(ASVs = unique_LT27_T2[,1])
unique_LT27_T2 <- merge(unique_LT27_T2, L2FC_t2_27, by = "ASVs")

unique_LT30_T2 <- as.data.frame(Unique$unique_LT30_T2)
unique_LT30_T2 <- as.data.frame(Unique$unique_LT30_T2) |> 
  transmute(ASVs = unique_LT30_T2[,1])
unique_LT30_T2 <- merge(unique_LT30_T2, L2FC_t2_30, by = "ASVs")

#### Unique T3

unique_LT25_T3 <- as.data.frame(Unique$unique_LT25_T3)
unique_LT25_T3 <- as.data.frame(Unique$unique_LT25_T3) |> 
  transmute(ASVs = unique_LT25_T3[,1])
unique_LT25_T3 <- merge(unique_LT25_T3, L2FC_t3_25, by = "ASVs")


unique_LT27_T3 <- as.data.frame(Unique$unique_LT27_T3) 
unique_LT27_T3 <- as.data.frame(Unique$unique_LT27_T3) |> 
  transmute(ASVs = unique_LT27_T3[,1])
unique_LT27_T3 <- merge(unique_LT27_T3, L2FC_t3_27, by = "ASVs")

unique_LT30_T3 <- as.data.frame(Unique$unique_LT30_T3)
unique_LT30_T3 <- as.data.frame(Unique$unique_LT30_T3) |> 
  transmute(ASVs = unique_LT30_T3[,1])
unique_LT30_T3 <- merge(unique_LT30_T3, L2FC_t3_30, by = "ASVs")

Everywhere$LT25_LT27_LT30_T2

Everywhere_T1 <- as.data.frame(Everywhere$LT25_LT27_LT30_T1)
Everywhere_T1 <- Everywhere_T1 |> 
  transmute(ASVs = Everywhere_T1[,1])
Everywhere_T1 <- merge(Everywhere_T1, log_FC_T1, by = "ASVs")

Everywhere_T2 <- as.data.frame(Everywhere$LT25_LT27_LT30_T2)
Everywhere_T2 <- Everywhere_T2 |> 
  transmute(ASVs = Everywhere_T2[,1])
Everywhere_T2 <- merge(Everywhere_T2, log_FC_T2, by = "ASVs")

Everywhere_T3 <- as.data.frame(Everywhere$LT25_LT27_LT30_T3)
Everywhere_T3 <- Everywhere_T3 |> 
  transmute(ASVs = Everywhere_T3[,1])
Everywhere_T3 <- merge(Everywhere_T3, log_FC_T3, by = "ASVs")





######################### T1 und T3 CBASS_temp 36 vs 30


allt_T1vsT3 <- filter(allt, Timepoint == "T1" | Timepoint == "T3")
allt_T1vsT3_36 <- filter(allt_T1vsT3, Comparison == "36 vs 30")
neg_allt_T1vsT3_36 <- filter(allt_T1vsT3_36, DAAs < 0)
pos_allt_T1vsT3_36 <- filter(allt_T1vsT3_36, DAAs >= 0)

ggplot(allt_T1vsT3_36, aes(x = Timepoint, y = DAAs, fill = Timepoint))+
  geom_bar(stat = "identity")+
  geom_text(neg_allt_T1vsT3_36, mapping = aes(label = DAAs), nudge_y = -20, size = 7)+
  geom_text(pos_allt_T1vsT3_36, mapping = aes(label = DAAs), nudge_y = 20, size = 7)+
  geom_hline(yintercept = 0)+
  facet_wrap(~ LT_temp) +
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face="bold"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

########################################################################
######################### Venn diagram
########################################################################

theme_clean2 <- function(base_size = 12,
                        base_family = "sans") {
  (
    theme_foundation(base_size = base_size,
                     base_family = base_family) + theme(
                       axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                     axis.text = element_text(size = ceiling(base_size * 0.7), colour = "black"),
                     axis.title = element_text(size = ceiling(base_size * 0.8)),
                     panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_blank(),
                      strip.background = element_rect(linetype = 0),
                      strip.text = element_text(),
                      strip.text.x = element_text(vjust = 0.5),
                      strip.text.y = element_text(angle = -90),
                      legend.text = element_text(size = ceiling(base_size * 0.9), family = "sans"),
                      legend.title = element_text(
                         size = base_size,
                         face = "bold",
                         family = "sans"
                       ),
                       legend.position = "right",
                       legend.key = element_rect(fill = "white", colour = NA),
                       legend.background = element_rect(colour = "black"),
                       plot.background = element_rect(colour = "black"),
                       plot.title = element_text(size = ceiling(base_size * 1.1), face = "bold"),
                       plot.subtitle = element_text(size = ceiling(base_size * 1.05))
                     )
  )
}




venn_maker3000 <- function(data){
  
  
  # for T1 high = "#003300", low = "#99FFCC"
  # for T2 high = "#993300", low = "#FFCC99"
  # for T3 high = "#330066", low = "#CC99FF"
  data <- data[c("LT_27", "LT_25", "LT_30")]
  venn <- Venn(data)
  data_venn <- process_data(venn)
  ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = count), data = venn_region(data_venn), show.legend = TRUE) +
    # 2. set edge layer
    geom_sf(color = "black", data = venn_setedge(data_venn), show.legend = FALSE, size = 4) +
    # 3. set label layer
    geom_sf_text(aes(label = c(27, 25, 30)), size = 10, data = venn_setlabel(data_venn),show.legend = FALSE) +
    # 4. region label layer
    geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                  data = venn_region(data_venn),
                  size = 5.5) +
    theme_clean2()+
    scale_fill_continuous(high = "#003300", low = "#99FFCC") +
    theme(legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))
}


venn_maker3001 <- function(data){
  data <- data[c("T2vsT1", "T3vsT1")]
  venn <- Venn(data)
  data_venn <- process_data(venn)
  ggplot() +
    # 1. region count layer
    geom_sf(aes(fill = count), data = venn_region(data_venn), show.legend = TRUE) +
    # 2. set edge layer
    geom_sf(color = "black", data = venn_setedge(data_venn), show.legend = FALSE, size = 4) +
    # 3. set label layer
    geom_sf_text(aes(label = c("T2vsT1", "T3vsT1")), size = 10, data = venn_setlabel(data_venn),show.legend = FALSE) +
    # 4. region label layer
    geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                  data = venn_region(data_venn),
                  size = 7) +
                  #position = "jitter") +
    theme_void()+
    scale_fill_continuous(high = "#003366", low = "#99CCFF") +
    theme(legend.text = element_text(size = 20),
          legend.title = element_text(size = 15),
          plot.background = element_rect(colour = "black"))
          #legend.key.height= unit(1, 'cm'))
}


# function to extract DAAs from a certain LT_temp from different time point from the prior created DAA tables (t1_25, t1_27, ...) 
per_LT_temp <- function(Timepoint, LT_temp, CBASS_temp, tax_type, x) {
  result_loop <- list()
  result_list <- list()
  for (t in Timepoint) {
    if (tax_type == "ASV"){
      if (t == 1) {
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be '25', '27' or '30'.")
        }
      }
      else if (t == 2){
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be '25', '27' or '30'.")
        }
      }
      else if (t == 3){
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be '25', '27' or '30'.")
        }
      }
      else{
        stop("Invalid value of 't'. It should be '1', '2' or '3'.")
      }
    }
    else if (tax_type =="Species" |
             tax_type =="Genus" | 
             tax_type =="Family" |
             tax_type =="Order" |
             tax_type =="Class" |
             tax_type =="Phylum" |
             tax_type =="Kingdom") {
      if (t == 1) {
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_25$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_25$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_25$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_25$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_25$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_25$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_25$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_25$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_25$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_27$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_27$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_27$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_27$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_27$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_27$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_27$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_27$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_27$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_30$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_30$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_30$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_30$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_30$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_30$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t1_30$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t1_30$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1_30$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be 'LT_25', 'LT_27' or 'LT_30'.")
        }
      }
      else if (t == 2){
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_25$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_25$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_25$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_25$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_25$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_25$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_25$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_25$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_25$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_27$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_27$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_27$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_27$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_27$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_27$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_27$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_27$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_27$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_30$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_30$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_30$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_30$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_30$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_30$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t2_30$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t2_30$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t2_30$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be 'LT_25', 'LT_27' or 'LT_30'.")
        }
      }
      else if (t == 3){
        if ( LT_temp == "LT_25"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_25$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_25$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_25$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_25$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_25$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_25$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_25$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_25$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t1t3_25_30$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_27"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_27$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_27$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_27$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_27$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_27$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_27$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_27$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_27$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_27$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( LT_temp == "LT_30"){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_30$increase_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_30$increase_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_30$increase_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_30$decrease_CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_30$decrease_CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_30$decrease_CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- t3_30$CBASS_33_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- t3_30$CBASS_36_vs_30[[tax_type]]
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- t3_30$CBASS_39_vs_30[[tax_type]]
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else{
          stop("Invalid value of 'LT_temp'. It should be 'LT_25', 'LT_27' or 'LT_30'.")
        }
      }
      else{
        stop("Invalid value of 't'. It should be '1', '2' or '3'.")
      }
    }
    else{
      stop("Invalid value of 'tax_type'. 
           It should be 'ASV', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum' or 'Kingdom'.")
    }
    
    
    
    result_list[[t]] <- data_swapped
  }
  #result_loop <- append(result_loop, result_list)
  result_loop <- result_list
  if (t == 1){
    names(result_loop) <- "T1"
  }
  else if (t == 2){
    names(result_loop) <- c("T1", "T2")
  }
  else if (t == 3){
    names(result_loop) <- c("T1", "T2", "T3")
  }
  print(result_loop)
}

T1_T2_T3_36 <- per_LT_temp(Timepoint = c(1, 2, 3),
            LT_temp = "LT_25",
            CBASS_temp = 36,
            tax_type = "ASV",
            x = "all")

venn_maker3000(T1_T2_T3_36)

# function to extract DAAs from a certain time point from different long term temperatures from the prior created DAA tables (t1_25, t1_27, ...) 
per_Timepoint <- function(Timepoint, LT_temp,CBASS_temp, tax_type, x) {
  result_loop <- list()
  result_list <- list()
  
  
  for (z in LT_temp) {
    if ( tax_type == "ASV"){
      if (z == "LT_25") {
        if ( Timepoint == 1){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase', 'decrease' or 'all'.")
          }
        }
        else if ( Timepoint == 2){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else if ( Timepoint == 3){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_25$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_25$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_25$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else{
          stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
        }
      }
      else if (z == "LT_27"){
        if ( Timepoint == 1){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else if ( Timepoint == 2){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else if ( Timepoint == 3){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_27$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_27$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_27$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else{
          stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
        }
      }
      else if (z == "LT_30"){
        if ( Timepoint == 1){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t1_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t1_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t1_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else if ( Timepoint == 2){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t2_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t2_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t2_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else if ( Timepoint == 3){
          if (x == "increase") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$increase_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$increase_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$increase_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if (x == "decrease") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$decrease_CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$decrease_CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$decrease_CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          } 
          else if(x == "all") {
            if (CBASS_temp == 33) {
              data_swapped <- row.names(t3_30$CBASS_33_vs_30)
            } 
            else if (CBASS_temp == 36) {
              data_swapped <- row.names(t3_30$CBASS_36_vs_30)
            } 
            else if (CBASS_temp == 39) {
              data_swapped <- row.names(t3_30$CBASS_39_vs_30)
            } 
            else {
              stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
            }
          }
          else {
            stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
          }
        }
        else{
          stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
        }
      }
      else{
        stop("Invalid value of 'z'. It should be 'LT_25', 'LT_27' or 'LT_30'.")
      }
    }
    else if (tax_type =="Species" |
             tax_type =="Genus" | 
             tax_type =="Family" |
             tax_type =="Order" |
             tax_type =="Class" |
             tax_type =="Phylum" |
             tax_type =="Kingdom") {
    if (z == "LT_25") {
      if ( Timepoint == 1){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_25$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_25$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_25$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_25$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_25$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_25$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_25$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_25$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_25$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 2){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_25$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_25$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_25$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_25$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_25$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_25$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_25$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_25$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_25$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 3){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_25$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_25$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_25$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_25$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_25$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_25$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_25$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_25$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_25$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else{
        stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
      }
    }
    else if (z == "LT_27"){
      if ( Timepoint == 1){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_27$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_27$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_27$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_27$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_27$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_27$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_27$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_27$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_27$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 2){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_27$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_27$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_27$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_27$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_27$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_27$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_27$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_27$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_27$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 3){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_27$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_27$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_27$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_27$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_27$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_27$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_27$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_27$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_27$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else{
        stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
      }
    }
    else if (z == "LT_30"){
      if ( Timepoint == 1){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_30$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_30$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_30$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_30$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_30$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_30$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t1_30$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t1_30$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t1_30$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 2){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_30$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_30$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_30$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_30$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_30$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_30$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t2_30$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t2_30$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t2_30$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else if ( Timepoint == 3){
        if (x == "increase") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_30$increase_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_30$increase_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_30$increase_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if (x == "decrease") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_30$decrease_CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_30$decrease_CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_30$decrease_CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        } 
        else if(x == "all") {
          if (CBASS_temp == 33) {
            data_swapped <- t3_30$CBASS_33_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 36) {
            data_swapped <- t3_30$CBASS_36_vs_30[[tax_type]]
          } 
          else if (CBASS_temp == 39) {
            data_swapped <- t3_30$CBASS_39_vs_30[[tax_type]]
          } 
          else {
            stop("Invalid value of 'CBASS_temp'. It should be 33, 36, or 39.")
          }
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase' or 'decrease'.")
        }
      }
      else{
        stop("Invalid value of 'Timepoint'. It should be '1', '2' or '3'.")
      }
    }
    else{
      stop("Invalid value of 'z'. It should be 'LT_25', 'LT_27' or 'LT_30'.")
    }
    }
    else{
      stop("Invalid value of 'tax_type'. 
           It should be 'ASV', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum' or 'Kingdom'.")
    }
    
    
    
    if (z == "LT_25"){
      number <- 1
    }
    else if (z == "LT_27" & length(LT_temp) == 1){
      number <- 1
    }
    else if (z == "LT_27" & length(LT_temp) == 2){
      number <- 2
    }
    else if (z == "LT_27" & length(LT_temp) == 3){
      number <- 2
    }
    else if (z == "LT_30" & length(LT_temp) == 1){
      number <- 1
    }
    else if (z == "LT_30" & length(LT_temp) == 2){
      number <- 2
    }
    else if (z == "LT_30" & length(LT_temp) == 3){
      number <- 3
    }

    
    
    
    result_list[[number]] <- data_swapped
  }

  result_loop <- result_list
  if (z == "LT_25"){
    names(result_loop) <- "LT_25"
  }
  else if (z == "LT_27" & length(LT_temp) == 1){
    names(result_loop) <- LT_temp
  }
  else if (z == "LT_27" & length(LT_temp) == 2){
    names(result_loop) <- LT_temp
  }
  else if (z == "LT_27" & length(LT_temp) == 3){
    names(result_loop) <- LT_temp
  }
  else if (z == "LT_30" & length(LT_temp) == 1){
    names(result_loop) <- LT_temp
  }
  else if (z == "LT_30" & length(LT_temp) == 2){
    names(result_loop) <- LT_temp
  }
  else if (z == "LT_30" & length(LT_temp) == 3){
    names(result_loop) <- LT_temp
  }
  print(result_loop)
}

T1_36_ASV <- per_Timepoint(Timepoint = 1, 
                              LT_temp = c("LT_25","LT_27", "LT_30"), 
                              CBASS_temp = 36, 
                              tax_type = "ASV", 
                              x = "all")

Venn_T3 <- venn_maker3000(T3_36_ASV)
Venn_T2 <- venn_maker3000(T2_36_ASV)
Venn_T1 <- venn_maker3000(T1_36_ASV)

DNA_LT_25_CBASS_30_function <- function(Timepoint, tax_type, x) {
  result_loop <- list()
  result_list <- list()
  for (t in Timepoint) {
    if (tax_type == "ASV"){
      if (t == "T2vsT1") {
          if (x == "increase") {
              data_swapped <- row.names(DNA_LT_25_CBASS_30$increase_Timepoint_T2_vs_T1)
          } 
        else if (x == "decrease") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$decrease_Timepoint_T2_vs_T1)
        } 
        else if(x == "all") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$Timepoint_T2_vs_T1)
            }
        else {
            stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
          }
        }
      else if (t == "T3vsT1"){
        if (x == "increase") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$increase_Timepoint_T3_vs_T1)
        } 
        else if (x == "decrease") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$decrease_Timepoint_T3_vs_T1)
        } 
        else if(x == "all") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$Timepoint_T3_vs_T1)
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
        }
      }
      else if (t == "T3vsT2"){
        if (x == "increase") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$increase_Timepoint_T3_vs_T2)
        } 
        else if (x == "decrease") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$decrease_Timepoint_T3_vs_T2)
        } 
        else if(x == "all") {
          data_swapped <- row.names(DNA_LT_25_CBASS_30$Timepoint_T3_vs_T2)
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
        }
      }
      else{
        stop("Invalid value of 't'. It should be 'T2vsT1', 'T3vsT1' or 'T3vsT2'.")
      }
    }
    else if (tax_type =="Species" |
             tax_type =="Genus" | 
             tax_type =="Family" |
             tax_type =="Order" |
             tax_type =="Class" |
             tax_type =="Phylum" |
             tax_type =="Kingdom") {
      if (t == "T2vsT1") {
        if (x == "increase") {
          data_swapped <- DNA_LT_25_CBASS_30$increase_Timepoint_T2_vs_T1[[tax_type]]
        } 
        else if (x == "decrease") {
          data_swapped <- DNA_LT_25_CBASS_30$decrease_Timepoint_T2_vs_T1[[tax_type]]
        } 
        else if(x == "all") {
          data_swapped <- DNA_LT_25_CBASS_30$Timepoint_T2_vs_T1[[tax_type]]
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
        }
      }
      else if (t == "T3vsT1"){
        if (x == "increase") {
          data_swapped <- DNA_LT_25_CBASS_30$increase_Timepoint_T3_vs_T1[[tax_type]]
        } 
        else if (x == "decrease") {
          data_swapped <- DNA_LT_25_CBASS_30$decrease_Timepoint_T3_vs_T1[[tax_type]]
        } 
        else if(x == "all") {
          data_swapped <- DNA_LT_25_CBASS_30$Timepoint_T3_vs_T1[[tax_type]]
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
        }
      }
      else if (t == "T3vsT2"){
        if (x == "increase") {
          data_swapped <- DNA_LT_25_CBASS_30$increase_Timepoint_T3_vs_T2[[tax_type]]
        } 
        else if (x == "decrease") {
          data_swapped <- DNA_LT_25_CBASS_30$decrease_Timepoint_T3_vs_T2[[tax_type]]
        } 
        else if(x == "all") {
          data_swapped <- DNA_LT_25_CBASS_30$Timepoint_T3_vs_T2[[tax_type]]
        }
        else {
          stop("Invalid value of 'x'. It should be 'increase','decrease' or 'all'.")
        }
      }
      else{
        stop("Invalid value of 't'. It should be 'T2vsT1', 'T3vsT1' or 'T3vsT2'.")
      }
    }
    else{
      stop("Invalid value of 'tax_type'. 
           It should be 'ASV', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum' or 'Kingdom'.")
    }
    
    
    
  if (t == "T2vsT1"){
    number <- 1
  }
  else if (t == "T3vsT1" & length(Timepoint) == 1){
    number <- 1
  }
  else if (t == "T3vsT1" & length(Timepoint) == 2){
    number <- 2
  }
  else if (t == "T3vsT1" & length(Timepoint) == 3){
    number <- 2
  }
  else if (t == "T3vsT2" & length(Timepoint) == 1){
    number <- 1
  }
  else if (t == "T3vsT2" & length(Timepoint) == 2){
    number <- 2
  }
  else if (t == "T3vsT2" & length(Timepoint) == 3){
    number <- 3
  }
  
  
  
  
  result_list[[number]] <- data_swapped
  }
  #result_loop <- append(result_loop, result_list)
  
    result_loop <- result_list
  if (t == "T2vsT1"){
    names(result_loop) <- "T2vsT1"
  }
  else if (t == "T3vsT1" & length(Timepoint) == 1){
    names(result_loop) <- Timepoint
  }
  else if (t == "T3vsT1" & length(Timepoint) == 2){
    names(result_loop) <- Timepoint
  }
  else if (t == "T3vsT1" & length(Timepoint) == 3){
    names(result_loop) <- Timepoint
  }
  else if (t == "T3vsT2" & length(Timepoint) == 1){
    names(result_loop) <- Timepoint
  }
  else if (t == "T3vsT2" & length(Timepoint) == 2){
    names(result_loop) <- Timepoint
  }
  else if (t == "T3vsT2" & length(Timepoint) == 3){
    names(result_loop) <- Timepoint
  }
  print(result_loop)

}



DNA_test_ASV_function <- DNA_LT_25_CBASS_30_function(c("T2vsT1", "T3vsT1"), "ASV", "all")
test_DNA <- venn_maker3001(DNA_test_ASV_function)
DNA_CBASS30_LT25_T1vs
nested <- (DNA_CBASS30_LT25_T1vs/test_DNA)+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested #view multi-panel figure


ASVs_T2vsT1_complete <- sort(test_ASV_function$T2vsT1)
ASVs_T3vsT1_complete <- sort(test_ASV_function$T3vsT1)


ASVs_T2vsT1_only <- sort(ASVs_T2vsT1_complete[!(ASVs_T2vsT1_complete %in% ASVs_T3vsT1_complete)])
ASVs_T3vsT1_only <- sort(ASVs_T3vsT1_complete[!(ASVs_T3vsT1_complete %in% ASVs_T2vsT1_complete)])
ASVs_T2vsT1_and_T3vsT1 <- sort(ASVs_T2vsT1_complete[ASVs_T2vsT1_complete %in% ASVs_T3vsT1_complete])

max_ln <- max(c(length(ASVs_T2vsT1_complete),
                length(ASVs_T3vsT1_complete),
                length(ASVs_T2vsT1_only),
                length(ASVs_T3vsT1_only),
                length(ASVs_T2vsT1_and_T3vsT1)))

df_T2vsT1_and_T3vsT1 <- data.frame(unique_T2vsT1 = c(ASVs_T2vsT1_only, rep(NA, max_ln -length(ASVs_T2vsT1_only))),
                                   unique_T3vsT1 = c(ASVs_T3vsT1_only, rep(NA, max_ln -length(ASVs_T3vsT1_only))),
                                   total_T2vsT1 = c(ASVs_T2vsT1_complete, rep(NA, max_ln -length(ASVs_T2vsT1_complete))),
                                   total_T3vsT1 = c(ASVs_T3vsT1_complete, rep(NA, max_ln -length(ASVs_T3vsT1_complete))),
                                   shared = c(ASVs_T2vsT1_and_T3vsT1, rep(NA, max_ln -length(ASVs_T2vsT1_and_T3vsT1))))



LT25_T3_30vs36_complete <- sort(rownames(t3_25$CBASS_36_vs_30))
LT27_T3_30vs36_complete <- sort(rownames(t3_27$CBASS_36_vs_30))
LT30_T3_30vs36_complete <- sort(rownames(t3_30$CBASS_36_vs_30))
LT25_T3_30vs36_only <- sort(LT25_T3_30vs36_complete[!(LT25_T3_30vs36_complete %in% c(LT27_T3_30vs36_complete, LT30_T3_30vs36_complete))])
LT27_T3_30vs36_only <- sort(LT27_T3_30vs36_complete[!(LT27_T3_30vs36_complete %in% c(LT25_T3_30vs36_complete, LT30_T3_30vs36_complete))])
LT30_T3_30vs36_only <- sort(LT30_T3_30vs36_complete[!(LT30_T3_30vs36_complete %in% c(LT25_T3_30vs36_complete, LT27_T3_30vs36_complete))])


LT25andLT27_T3_30vs36 <- sort(LT25_T3_30vs36_complete[LT25_T3_30vs36_complete %in% LT27_T3_30vs36_complete & !(LT25_T3_30vs36_complete %in% LT30_T3_30vs36_complete)])
LT25andLT30_T3_30vs36 <- sort(LT25_T3_30vs36_complete[LT25_T3_30vs36_complete %in% LT30_T3_30vs36_complete & !(LT25_T3_30vs36_complete %in% LT27_T3_30vs36_complete)])
LT27andLT30_T3_30vs36 <- sort(LT27_T3_30vs36_complete[LT27_T3_30vs36_complete %in% LT30_T3_30vs36_complete & !(LT27_T3_30vs36_complete %in% LT25_T3_30vs36_complete)])

shared <- sort(LT25_T3_30vs36_complete[LT25_T3_30vs36_complete %in% LT27_T3_30vs36_complete & LT25_T3_30vs36_complete %in% LT30_T3_30vs36_complete])

max_ln <- max(c(length(LT25_T3_30vs36_complete),
                length(LT27_T3_30vs36_complete),
                length(LT30_T3_30vs36_complete),
                length(LT25_T3_30vs36_only),
                length(LT27_T3_30vs36_only),
                length(LT30_T3_30vs36_only),
                length(LT25andLT27_T3_30vs36),
                length(LT25andLT30_T3_30vs36),
                length(LT27andLT30_T3_30vs36),
                length(shared)))

df_T3 <- data.frame(unique_LT25 = c(LT25_T3_30vs36_only, rep(NA, max_ln -length(LT25_T3_30vs36_only))),
                    unique_LT27 = c(LT27_T3_30vs36_only, rep(NA, max_ln -length(LT27_T3_30vs36_only))),
                    unique_LT30 = c(LT30_T3_30vs36_only, rep(NA, max_ln -length(LT30_T3_30vs36_only))),
                    total_LT25 = c(LT25_T3_30vs36_complete, rep(NA, max_ln -length(LT25_T3_30vs36_complete))),
                    total_LT27 = c(LT27_T3_30vs36_complete, rep(NA, max_ln -length(LT27_T3_30vs36_complete))),
                    total_LT30 = c(LT30_T3_30vs36_complete, rep(NA, max_ln -length(LT30_T3_30vs36_complete))),
                    LT25_LT27 = c(LT25andLT27_T3_30vs36, rep(NA, max_ln -length(LT25andLT27_T3_30vs36))),
                    LT25_LT30 = c(LT25andLT30_T3_30vs36, rep(NA, max_ln -length(LT25andLT30_T3_30vs36))),
                    LT27_LT30 = c(LT27andLT30_T3_30vs36, rep(NA, max_ln -length(LT27andLT30_T3_30vs36))),
                    LT25_LT27_LT30 = c(shared, rep(NA, max_ln -length(shared)))
)

write.xlsx(df_T3,"output directory",colNames = TRUE)


######################### Data frame of ASVs with taxonomic annotations
LT25_T1 <- as.data.frame(t1_25[["CBASS_36_vs_30"]])






a <- as.data.frame(LT_25_CBASS_30$Timepoint_T2_vs_T1)[,7:13]
a <- tibble::rownames_to_column(a, "ASV")
b <- as.data.frame(LT_25_CBASS_30$Timepoint_T3_vs_T1)[,7:13]
b <- tibble::rownames_to_column(b, "ASV")
c <- as.data.frame(LT_25_CBASS_30$Timepoint_T3_vs_T2)[,7:13]
c <- tibble::rownames_to_column(c, "ASV")
e <- rbind(a,b,c)
e_unique <- unique(e)
e_unique <- e_unique %>% arrange(ASV)
colnames(e_unique)[1] <- "DAA"
schnittmenge_a_b <- semi_join(a, b, by = "ASV")
schnittmenge_a_b <- schnittmenge_a_b %>% arrange(ASV)
colnames(schnittmenge_a_b)[1] <- "DAA"
write.csv(schnittmenge_a_b, "output directory", row.names=T)


########################################################################
### Patchwork plot
########################################################################
all_LT_Time
Venn_T1
Venn_T2
Venn_T3
nested <- (all_LT_Time/(Venn_T1|Venn_T2|Venn_T3))+
  plot_annotation(tag_levels = 'A') & #add figure labels
  theme(plot.tag = element_text(size = 25))#change tag size
nested #view multi-panel figure

plot_grid(all_LT_Time, Venn_T1, Venn_T2,Venn_T3 )
?plot_grid

########################################################################
### Permanova
########################################################################



rare_RNA <- subset_samples(rare_phy, type == "RNA")
relative_rare <- microbiome::transform(rare_RNA, 'compositional')

ps_dist_matrix <- phyloseq::distance(relative_rare, method ="bray")
metadata <- as(sample_data(rare_RNA), "data.frame")
z <- adonis2(ps_dist_matrix ~ LT_temp*Timepoint+Genotype+CBASS_temp,
        data = metadata)
z <- cbind(" "=rownames(z), z)
write.xlsx(z,"output directory",colNames = TRUE)


rare_RNA <- subset_samples(rare_phy, type == "RNA")
rare_RNA_CBASS_30 <- subset_samples(rare_RNA, CBASS_temp == 30)
relative_rare_30 <- microbiome::transform(rare_RNA_CBASS_30, 'compositional')

ps_dist_matrix_30 <- phyloseq::distance(relative_rare_30, method ="bray")
metadata_30 <- as(sample_data(rare_RNA_CBASS_30), "data.frame")
z_30 <- adonis2(ps_dist_matrix_30 ~ LT_temp*Timepoint+Genotype,
             data = metadata_30)

z_30 <- cbind(" "=rownames(z_30), z_30)
write.xlsx(z_30,"output directory",colNames = TRUE)


rare_RNA_T1 <- subset_samples(rare_RNA, Timepoint == "T3")
rare_RNA_T1_LT25 <- subset_samples(rare_RNA, LT_temp == "30")
relative_rare_T1 <- microbiome::transform(rare_RNA_T1_LT25, 'compositional')
ps_dist_matrix_T1 <- phyloseq::distance(relative_rare_T1, method ="bray")
#dim(ps_dist_matrix_T1)
#dim(relative_rare_T1)
pair_T1_27 <- pairwise.adonis(ps_dist_matrix_T1, sample_data(relative_rare_T1)$CBASS_temp)
write.xlsx(pair_T1_27,"output directory",colNames = TRUE)





RNA <- subset_samples(rare_phy, type == "RNA")
metadata <- as(sample_data(RNA), "data.frame")
adonis2(distance(RNA, method="bray") ~ LT_temp*Timepoint+CBASS_temp+Genotype,
       data = metadata)
RNA.25 <- subset_samples(RNA, LT_temp == 25)
distance.RNA <- phyloseq::distance(RNA.25, method = "bray")
pairwise.adonis(distance.RNA, sample_data(RNA.25)$LT_temp)



########################################################################
### prepare for Picrust
########################################################################

ASVs <- rownames(tax_table(phy))
all_ASVs <- rownames(asv)


# Neuen Vektor erstellen mit den gewünschten Änderungen
rm(Vektor2)
Vektor2 <- sub("ASV", "ASV_", ASVs)
Vektor2 <- sub("ASV_(0+)", "ASV_", Vektor2)
all_ASVs_without0 <- Vektor2
print(Vektor2)[c(17000:17060)]

Vektor3 <- sub("ASV", "ASV_", ASVs)
Vektor3 <- sub("ASV_(0+)", "ASV_", Vektor3)
ASVs_without0 <- Vektor3
print(Vektor3)

writeLines(ASVs_without0, "output directory")

otu <- t(as(otu_table(rare_phy), "matrix"))
rownames(otu) <- Vektor2
ASV_biom <- make_biom(data = otu)
write_biom(ASV_biom, "output directory")


metadata_rare <- as.data.frame(sample_data(rare_phy))
write.table(metadata_rare, "output directory", sep = "\t", row.names = T)
write_delim(metadata_rare, "output directory", delim = "\t")




