#!/usr/bin/Rscript
library(dada2); packageVersion("dada2")
path <- "/home/pt40963/misc_python/THAPBI/THAPBI-pycits/data/reads"
filtpath <- file.path(path, "filtered")
fns <- list.files(path, pattern="primers_trimmed.fastq.gz")
filterAndTrim(file.path(path,fns),  file.path(filtpath,fns),  maxEE=2, truncQ=11, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)
filtpath <- "/home/pt40963/misc_python/THAPBI/THAPBI-pycits/data/reads/filtered"
filts <- list.files(filtpath, pattern="primers_trimmed.fastq.gz", full.names=TRUE)
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
names(filts) <- sample.names
set.seed(100)
err <- learnErrors(filts, multithread=TRUE, randomize=TRUE)
plotErrors(err, nominalQ=TRUE)
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
	cat("Processing:", sam, "\n")
	derep <- derepFastq(filts[[sam]])
	dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
seqtab <- makeSequenceTable(dds[sam])
saveRDS(seqtab, "/mnt/shared/users/pt40963/misc_python/THAPBI/THAPBI-pycits/bin/dada2/seqtab.rds")
st1 <- readRDS("/mnt/shared/users/pt40963/misc_python/THAPBI/THAPBI-pycits/bin/dada2/seqtab.rds")
seqtab <- removeBimeraDenovo(st1, method="pooled", multithread=TRUE)
write.table(seqtab, "sequence_table.txt")
