# I used the script from https://benjjneb.github.io/dada2/tutorial.html and modified it for my purpose.

################################
# setwd
# load needed packages
################################

setwd("your working directory")
library(dada2)

################################
# Create a path(), the path is the location of the fastq files
# printing the file names
################################

path <- "your output location path"
list.files(path) # prints file names

################################
# Read the files in
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fq.gz and SAMPLENAME_2.fq.gz
################################

fnFs <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_1\\."), `[`, 1) # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#cat("Processing",length(sample.names),"samples:", sample.names)

#Inspect read quality profiles
#plotQualityProfile(fnFs[1:3])
#plotQualityProfile(fnRs[1:3])

################################
# 1. Filter and trim
################################

cat("Filtering and trimming")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
### set parameters  maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out) # shows first 10 lines of output


################################
# 2. Dereplicate 
# Dereplication combines all identical sequencing reads into into 
# “unique sequences” with a corresponding “abundance” equal to the number of 
# reads with that unique sequence. 
# Dereplication substantially reduces computation time by eliminating redundant comparisons.
################################

derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)


################################
# 3. Learn error rates
################################

cat("Learning error rates")
errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepR1, multithread=TRUE)


################################
# 4. Sample inference
# apply the core sample inference algorithm to the filtered, trimmed and
# dereplicated sequence data.
################################ 

dadaFs <- dada(derepF1, err=errF, multithread=TRUE,pool = F)
dadaRs <- dada(derepR1, err=errR, multithread=TRUE,pool = F)

#Inspecting the returned dada-class object:
#dadaFs[[1]]


################################
# 5. Merge paired reads
################################
mergers <- mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose=TRUE)


################################ 
# step 6. Construct sequence table
################################ 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:265]
dim(seqtab2)



################################ 
# 7. Remove chimeras
################################ 
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # how many merged sequence variants are chimeric 

# calculate the frequency of chimeric sequences. Out of all merged sequence reads how many percent are chimeras
sum(seqtab.nochim)/sum(seqtab2)

################################ 
# Track reads through the pipeline
################################ 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

################################ 
#  8. Assign taxonomy
################################ 

taxa <- assignTaxonomy(seqtab.nochim, "your path for taxonomy", multithread=TRUE) # I used silva_nr_v138_train_set.fa
taxa <- addSpecies(taxa, "your path for taxonomy") # I used silva_species_assignment_v138.fa
#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


################################ 
# 9. Export table
################################ 

samples.out <- rownames(seqtab.nochim)

asv.1=t(seqtab.nochim)
asv.1=cbind(asv.1, "sum"=rowSums(asv.1)) 
asv.2=merge(asv.1, taxa, by="row.names")
colnames(asv.2)[1]="Sequence"
asv.3=asv.2[,c(2:ncol(asv.2),1)]
asv.final=asv.3[order(-asv.3$sum),]
rownames(asv.final) = sprintf("ASV%04d", 1:nrow(asv.final))
write.table(asv.final, file ="output directory",  quote = FALSE)
write.table(track, file ="output directory", quote = FALSE)
write.table(seqtab.nochim, file ="output directory")

################################ 
# 10. Preparation for Picrust
################################ 

#### Change seq headers to asv number (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#### Create and write out a fasta of our final ASV seqs

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/home/maischt/2023.01_spis_OWCBASS_16S/ASVs.fa")

#### ASV count table

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "output directory", sep="\t", quote=F, col.names=NA)

#### Tax table

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "output directory", sep = "\t", quote=F, col.names=NA)

#### To merge asv abundance and taxonomy into one file

ASV_TAX_table <- merge(asv_tab, asv_tax, by=0)
write.table(ASV_TAX_table, "output directory", sep = "\t", quote=F, col.names=NA)


