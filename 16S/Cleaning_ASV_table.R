setwd("your working directory")
library(openxlsx)
library(stringr)
library(dplyr)


##################################################################
#####Identifying and removing contaminant ASVs normalized data#####
###################################################################

asv=read.table("ASV_table_raw.txt", header = TRUE, row.names = 1)[,1:205] # A dataframe with ASVs as rows and treatment as cloumns showing the ASV abundance
# I didn't included row 206!!! Those are the sums
tax=read.table("ASV_table_raw.txt", header = TRUE, row.names = 1)[,207:213] # A dataframe with the taxonomy of the ASVs
#map=read.table("inputs/Microbiome_16S_tutorial_metadata.txt", header = T, row.names = 1, sep = "\t")
#match("Endozoicomonadaceae",tax$Family)
# identify  samples with low sequencing depth
hist(colSums(asv),  breaks = 50, labels = F)#first plot the distribution of reads per sample
colSums(asv) # see the number of reads per sample name
#sum(asv$T1S4_11wkLT30_._30_16SDNA_a113)
asv.o=asv[, colSums(asv) > 5000]#to remove samples with less than 5000 reads
dim(asv.o)#number of observations (ASVs) and columns (samples)
message(ncol(asv.o)," samples with > 5000 reads were retained out of ", ncol(asv), " total samples")
print(colnames(asv[, colSums(asv) < 5000]))# Samples that were removed 
# subset(asv, select = c("T1S4_11wkLT30_._30_16SDNA_a113", "T2S3_24wkLT27_._39_16SRNA_R212")) |> 
#sum(asv$sum)

#Identify and removing contaminant ASVs raw data
asv.r=as.data.frame(sweep(asv.o,2,colSums(asv.o),`/`)) # calculate relative abundances on each sample
asv.r$Mean=rowMeans(asv.r[,1:ncol(asv.r)]) # create new column with the average
names(asv.r)
asv.r$NegMean=rowMeans(asv.r[,1:2]) # choose columns with the control samples
asv.r$contaFactor=(asv.r$NegMean/asv.r$Mean)*100
Conta=subset(asv.r, asv.r$contaFactor > 100)
Conta$Family=tax$Family[match(rownames(Conta), rownames(tax))]
message("Number of total ASVs: ", nrow(asv)) #17060 total ASVs
message("Number of identified contaminant ASVs removed from the analysis: ", 
        length(rownames(Conta)), "\n", Conta$Family[1],"\n", 
        Conta$Family[2],"\n", Conta$Family[3],"\n", 
        Conta$Family[4],"\n", Conta$Family[5])

#remove any chloroplast or mitochobdria and the control samples
unwant_tax=tax %>% filter_all(any_vars(str_detect(., 'Mitochondria|Chloroplast')))
message("Number of Mitochondria or Chloroplasts: ", length(rownames(unwant_tax)))
colnames(asv.o)
asv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(Conta) & !rownames(asv.o) %in% rownames(unwant_tax))[, -grep("PCR|Spcr_neg", colnames(asv.o))]
#asv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(Conta) & !rownames(asv.o) %in% rownames(unwant_tax))[, -grep("PCR", colnames(asv.o))][, -grep("Spcr_neg", colnames(asv.o))]

colnames(asv.noConta)
dim(asv.noConta)

#remove ASVs with only zeros
asv.noRare=asv.noConta[rowSums(asv.noConta[])>0,]
message("Number of ASVs used in the analysis: ", length(rownames(asv.noRare)))

# Export QC filtered ASV tables
asv.noConta.f=merge(asv.noRare, tax, by="row.names")
#write.table(asv.noConta.f, "output directory",  quote = FALSE, row.names=T, sep = "\t") 
#write.table(Conta, "output directory",  quote = FALSE, row.names=F, sep = "\t")


############
# subset asv.noConta.f for ASVs under 500 ####
############

asv_500 <- asv.noConta.f |> filter(rowSums(asv.noConta.f[,2:204]) > 500)
#abc <- colSums(asv_500[,c(2:204)])
#sort(abc, decreasing = T)
#asv_500 <- asv_500 |> select(-"T1S4_11wkLT30_._30_16SDNA_a113")
#abc <- colSums(asv_500[,c(2:203)])
#sort(abc, decreasing = T)

write.table(asv_500, "output directory",  quote = FALSE, row.names=F, sep = "\t")

