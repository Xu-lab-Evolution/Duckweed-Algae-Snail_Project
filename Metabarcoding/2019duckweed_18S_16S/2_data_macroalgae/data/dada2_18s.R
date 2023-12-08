library("tidyverse")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")
library("ape")
library("seqinr")

cat("Getting ready ...", format(Sys.time(), "%c"), "\n")

n_cores <- 8
indir <- 'data/18s'
outdir <- 'output/18s'
filter_dir <- 'output/18s/filt_data'
taxdir <- 'PR2'
tax_key <- 'PR2/pr2.fa.gz'
metadata_file <- 'metadata.txt'
fastaname <- 'asv.fasta'

dir.create('PR2', showWarnings = FALSE, recursive = TRUE)
download.file("https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_dada2.fasta.gz", file.path('PR2', 'pr2.fa.gz'))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)

metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE, row.names = 1)

cat("Beginning analyses", indir, "...", format(Sys.time(), "%c"), "\n")

fnFs <- sort(list.files(indir, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(indir, pattern="_R2_001.fastq.gz", full.names = TRUE))
sampleIDs <- sub("(.*?)-.*", "\\1", basename(fnFs))

filtFs <- file.path(filter_dir, paste0(sampleIDs, "_F_filt.fastq.gz"))
filtRs <- file.path(filter_dir, paste0(sampleIDs, "_R_filt.fastq.gz"))

names(filtFs) <- sampleIDs
names(filtRs) <- sampleIDs

cat("Filtering and trimming ...", format(Sys.time(), "%c"), "\n")

filter_results <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 200,
                                rm.phix = FALSE,
                                multithread = n_cores, 
                                compress = FALSE, verbose = TRUE)

exists <- file.exists(fnFs) & file.exists(filtFs)
filtFs <- filtFs[exists]
exists <- file.exists(fnRs) & file.exists(filtRs)
filtRs <- filtRs[exists]


cat("Learning errors ...", format(Sys.time(), "%c"), "\n")

errF <- learnErrors(filtFs, multithread = n_cores, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = n_cores, verbose = TRUE)


cat("Dereplicating ...", format(Sys.time(), "%c"), "\n")

dadaFs <- dada(filtFs, err=errF, pool = FALSE, multithread = n_cores)
dadaRs <- dada(filtRs, err=errR, pool = FALSE, multithread = n_cores)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


cat("Creating ASV table ...", format(Sys.time(), "%c"), "\n")

seqtab_all <- makeSequenceTable(mergers)


cat("Removing chimeras ...", format(Sys.time(), "%c"), "\n")

seqtab <- removeBimeraDenovo(seqtab_all,
                             method = 'consensus',
                             multithread = n_cores,
                             verbose = TRUE)


cat("Assigning taxonomy ...", format(Sys.time(), "%c"), "\n")

taxa <- assignTaxonomy(seqtab, tax_key, multithread = n_cores)

colnames(taxa) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', "Species")

otutable <- seqtab
colnames(otutable) <- paste('ASV', 1:ncol(seqtab), sep = '_')
row.names(taxa) <- paste('ASV', 1:ncol(seqtab), sep = '_')

asv_seqs <- colnames(seqtab)
asv_headers <- paste('ASV', 1:ncol(seqtab), sep = '_')
asv_fasta <- c(rbind(asv_headers, asv_seqs))

cat("Aligning sequences ...", format(Sys.time(), "%c"), "\n")

seqs <- getSequences(seqtab)
names(seqs) <- asv_headers
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = n_cores)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
 
cat("Calculating tree ...", format(Sys.time(), "%c"), "\n")

dna_dist <- dist.ml(phang.align, model="JC69")
asv_UPGMA <- upgma(dna_dist)
fit <- pml(asv_UPGMA, phang.align)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
tree <- bootstrap.pml(fitJC, bs=1, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))

cat("Saving ASV sequences ...", format(Sys.time(), "%c"), "\n")

seqs <- getSequences(seqtab)
write.fasta(as.list(seqs), asv_headers, file.path(outdir, fastaname))

cat("Creating phyloseq object ...", format(Sys.time(), "%c"), "\n")
  
ps <- phyloseq(otu_table(otutable, taxa_are_rows = FALSE),
                 sample_data(metadata_df),
                 tax_table(taxa))
  
getN <- function(x) sum(getUniques(x))
track <- cbind(filter_results, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleIDs

rm(list=setdiff(ls(), c("ps", "track")))
save.image(file = 'output/18s/ASV_table_final.rds')
