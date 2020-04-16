# Title     : train_idtaxa
# Objective : Train the idtaxa classifier
# Created by: jshleap from decipher documentation
# Created on: 2019-06-04

suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

print(paste('running on DECIPHER version', packageVersion("DECIPHER")))


### Start of the script
option_list = list(
    make_option(c("-f", "--fasta"), type="character", default=NULL,
    help="Name of training fasta file [default= %default]"),
    make_option(c("-p", "--prefix"), type="character", default='Train',
    help="Prefix for outputs [default= %default]"),
    make_option(c("-i", "--maxIterations"), type='integer', default=3,
    help="Number of iteration during training for problem seqeunce removal")
)

remove_nokmers <- function(train){
    K <- floor(log(100*quantile(width(train), 0.99), 4))
    kmers <- .Call("enumerateSequence", train, K, PACKAGE="DECIPHER")
    kmers <- lapply(kmers, function(x) sort(unique(x + 1L), na.last=NA))
    bad <- lengths(kmers)==0
    return(bad)
}

joining <- function (x){
    good <- paste(strsplit(x," ")[[1]][2:3], collapse=' ')
    if (!grepl("Root;", good)){
    good <- paste0("Root;",good, collapse='')
    }
    return(good)
}

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

taxid <- NULL
if (!file.exists(paste0(opt$prefix, '_seqs.rda'))){
    fas <- opt$fasta
    seqs <- readDNAStringSet(fas)
    seqs_clean <- RemoveGaps(seqs, processors = detectCores())
    save(seqs_clean, file=paste0(opt$prefix, '_seqs.rda'))
} else{
    load(paste0(opt$prefix, '_seqs.rda'))
}
groups <- names(seqs_clean)
groups <- gsub(";$", "", groups)
groups <- gsub("(.*)(Root;)", "\\2", groups)
groups <- unlist(lapply(groups, joining))
groupCounts <- table(groups)
u_groups <- names(groupCounts)
print(paste("There are", length(u_groups), "unique taxa in training file"))
maxGroupSize <- 10
print(paste("Keeping a maximum of", maxGroupSize, "sequences per taxa"))
remove <- logical(length(seqs_clean))
for (i in which(groupCounts > maxGroupSize)) {
index <- which(groups==u_groups[i])
keep <- sample(length(index),
maxGroupSize)
remove[index[-keep]] <- TRUE
}
print(paste(sum(remove), "sequences removed"))
train <- seqs_clean[!remove]
taxonomy <- groups[!remove]
bad <- remove_nokmers(train)
train <- train[!bad]
taxonomy <- taxonomy[!bad]
maxIterations <- opt$maxIterations
probSeqsPrev <- integer()
remove <- logical(length(train))
for (i in seq_len(maxIterations)) {
    cat("Training iteration: ", i, "\n", sep="")
    # train the classifier
    trainingSet <- LearnTaxa(train, taxonomy)
    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs)==0) {
        cat("No problem sequences remaining.\n")
        break
    } else if (length(probSeqs)==length(probSeqsPrev) &&
    all(probSeqsPrev==probSeqs)) {
        cat("Iterations converged.\n")
        break
    }
    if (i==maxIterations) {
        break} else {
    probSeqsPrev <- probSeqs
    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE # remove all problem sequences
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
    }
    print(paste(sum(remove), "problem sequences removed (cumulative)"))
    }
}
print(paste(sum(remove), "problem sequences removed (cumulative)"))
# trainingSet <- LearnTaxa(train, taxonomy)
# trainingSet <- LearnTaxa(seqs, groups)
pdf(paste0(opt$prefix,"_rplot.pdf"))
plot(trainingSet)
dev.off()
save(trainingSet, file=paste0(opt$prefix, '_idtaxa.rda'))