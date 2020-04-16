# Title     : TODO
# Objective : TODO
# Created by: jshleap
# Created on: 2019-08-22
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tools))

print(paste('running on DECIPHER version', packageVersion("DECIPHER")))
packages <- c('DECIPHER', 'parallel', 'foreach', 'doParallel', 'Biostrings')

### Start of the script
option_list = list(
    make_option(c("-f", "--fasta"), type="character", default=NULL,
    help="Name of training fasta file [default= %default]"),
    make_option(c("-p", "--prefix"), type="character", default='Train',
    help="Prefix for outputs [default= %default]"),
    make_option(c("-c", "--cpus"), type="double", default=1,
    help="Number of processors to use [default= %default]"),
    make_option(c("-r", "--toprank"), type="character", default='Eukaryota',
    help="Top rank in the sequences names [default= %default]"),
    make_option(c("-m", "--maxiter"), type="double", default=3,
    help="max iterations to remove problem sequences [default= %default]")

)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cl <- makeCluster(opt$cpus)
registerDoParallel(cl)
prefix <- file_path_sans_ext(basename(opt$fasta))
if (!file.exists(paste0(prefix,'seqs.rds'))){
    seqs_path <- opt$fasta
    system(paste('seqkit', '-j', opt$cpus, 'split2', '-p', opt$cpus, opt$fasta,
    '-f'))
    l <- foreach(i=Sys.glob(paste0(seqs_path, ".split/*.part_*")),
    .packages=packages) %dopar% readDNAStringSet(i)
    seqs <- do.call(c, l)
    saveRDS(seqs, paste0(prefix,'seqs.rds'))
} else{
    seqs <- readRDS(paste0(prefix,'seqs.rds'))
}

groups <- names(seqs)
groups <- gsub(paste0("(.*)(", opt$toprank,";)"), "\\2", groups)
groups <- gsub(" ", "_", groups)
if (opt$toprank != 'Root'){groups <- paste0('Root;', groups)}
groupCounts <- table(groups)
u_groups <- names(groupCounts)
maxGroupSize <- ceiling(mean(groupCounts))
print(paste("There are", length(u_groups), "unique taxa in training file"))
print(paste("Keeping a maximum of", maxGroupSize, "sequences per taxa"))
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
    index <- which(groups==u_groups[i])
    keep <- sample(length(index),
    maxGroupSize)
    remove[index[-keep]] <- TRUE
}
print(paste(sum(remove), "sequences removed"))
maxIterations <- opt$maxiter
allowGroupRemoval <- FALSE
probSeqsPrev <- integer()
for (i in seq_len(maxIterations)) {
    cat("Training iteration: ", i, "\n", sep="")
    # train the classifier
    trainingSet <- LearnTaxa(seqs[!remove],
    groups[!remove])
    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
        cat("No problem sequences remaining.\n")
        break
        } else if (length(probSeqs)==length(probSeqsPrev) &&
        all(probSeqsPrev==probSeqs)) {
        cat("Iterations converged.\n")
        break}
    if (i==maxIterations){break}
    probSeqsPrev <- probSeqs
    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    if (!allowGroupRemoval) {
        # replace any removed groups
        missing <- !(u_groups %in% groups[!remove])
        missing <- u_groups[missing]
        if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
        }
    }
}
print(paste(sum(remove), "sequences removed during training"))
print(paste(length(probSeqs), "problem sequences remaining"))
saveRDS(trainingSet, paste0(prefix, ".rds"))