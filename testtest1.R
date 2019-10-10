#script for DADA2 analysis in R

#install package
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.12")

#load package
library("dada2")

#load ur data
path<- "~/qbb2/dada2_analysis/reads/" #path to directory containing all the FASTA/Q files
list.files(path) #returns list of all files
fnFs<-sort(list.files(path, pattern="*_R1_001.fastq", full.names = TRUE)) #list all forward
fnRs<-sort(list.files(path,pattern="*_R2_001.fastq", full.names = TRUE)) #list all reverse reads
length(fnRs) #just checking I have all my 200 samples
length(fnFs)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #extract sample names
fnFs[1:2]


#view quality profiles
plotQualityProfile(fnFs[1:2])

#filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
View(filtFs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,260),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose = FALSE) # On Windows set multithread=FALSE
head(out)

#learn error r8
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#denoise stats
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


##extra analysis: bar plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")

##find where's the file with info on clusters
dada-class$clustering (or birth_subs)

