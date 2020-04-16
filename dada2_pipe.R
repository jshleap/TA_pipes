#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dada2));
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(phylotools))
print(paste('running on DADA2 version', packageVersion("dada2"),
'and DECIPHER version', packageVersion("DECIPHER")))


# I am assuming we have a demultiplex samples if not go to bcl2fastq.
# This will also only run with pairend (at least for now) remove adapter,
# linker etc...
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

process_df <- function(df, minqual, wdw, tr=Null){
    if (is.null(tr)){
        l <- max(df$Cycle)
        df$windows <- cut(df$Cycle, l/wdw)
        df$a <- df$Count * df$Score
         av <- ddply(df,~windows,summarise, Median=median(Score),
                    Mean=sum(a)/sum(Count), Mean2=mean(Score), q25=quantile(Score)[1],
                    Min=min(Score), mode=getmode(Score), mCycle=max(Cycle))
        save(av, df, file = 'df.rda')
        # print(av)
        trimm <- min(av[av$Mean  < minqual,]$mCycle)}
    else{trimm <- tr}
  return(c(trimm, l))
}
### Some functions
quality_plot <- function (fl, n = 5e+05, aggregate = FALSE, minqual = 20, wdw = 10, tr = NULL) {
  # modification of the dada2 plot to save the plot to file and output the dataframe as well
    outpref <- strsplit(basename(fl), ".fastq")[[1]][1]
    statdf <- data.frame(Cycle = integer(0), Mean = numeric(0), 
                         Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0), 
                         file = character(0))
    anndf <- data.frame(minScore = numeric(0), label = character(0), 
                        rclabel = character(0), rc = numeric(0), file = character(0))
    FIRST <- TRUE
    for (f in fl[!is.na(fl)]) {
      srqa <- qa(f, n = n)
      df <- srqa[["perCycle"]]$quality
      rc <- srqa[["readCounts"]]$read
      if (rc >= n) {
        rclabel <- paste("Reads >= ", n)
      }
      else {
        rclabel <- paste("Reads: ", rc)
      }
      means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                            df$Cycle)
      get_quant <- function(xx, yy, q) {
        xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
      }
      q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                       foo$Count, 0.25), simplify = TRUE)
      q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                       foo$Count, 0.5), simplify = TRUE)
      q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                       foo$Count, 0.75), simplify = TRUE)
      if (!all(sapply(list(names(q25s), names(q50s), names(q75s)), 
                      identical, rownames(means)))) {
        stop("Calculated quantiles/means weren't compatible.")
      }
      if (FIRST) {filterAndTrim
        plotdf <- cbind(df, file = f)
        FIRST <- FALSE
      }
      else {
        plotdf <- rbind(plotdf, cbind(df, file = f))
      }
      statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)), 
                                         Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                         Q75 = as.vector(q75s), file = f))
      anndf <- rbind(anndf, data.frame(minScore = min(df$Score), 
                                       label = basename(f), rclabel = rclabel, rc = rc, 
                                       file = f))
    }
    anndf$minScore <- min(anndf$minScore)
    if (aggregate) {
      plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, 
                                  sum)
      means <- rowsum(plotdf.summary$Score * plotdf.summary$Count, 
                      plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, 
                                                   plotdf.summary$Cycle)
      q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                               foo$Count, 0.25), simplify = TRUE)
      q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                               foo$Count, 0.5), simplify = TRUE)
      q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, 
                                                                               foo$Count, 0.75), simplify = TRUE)
      statdf.summary <- data.frame(Cycle = as.integer(rownames(means)), 
                                   Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                   Q75 = as.vector(q75s))
      #save(plotdf, "plotdf.rda")
      trim <- process_df(plotdf.summary, minqual, wdw, tr=tr)
      trimm <- trim[1]
      read_lenght <- trim[2]
      ggplot(data = plotdf.summary, aes(x = Cycle, y = Score)) + 
        geom_tile(aes(fill = Count)) + scale_fill_gradient(low = "#F5F5F5", high = "black") + 
        geom_line(data = statdf.summary, aes(y = Mean), color = "#66C2A5") + 
        geom_line(data = statdf.summary, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        geom_line(data = statdf.summary, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
        geom_line(data = statdf.summary, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        ylab("Quality Score") + xlab("Cycle") + 
        annotate("text", x = 0, y = min(anndf$minScore) + 2, label = sprintf("%d files (aggregated)", nrow(anndf)), 
                 hjust = 0) + 
        annotate("text", x = 0, y = min(anndf$minScore) + 1, label = sprintf("Total reads: %d", sum(anndf$rc)), 
                 hjust = 0) + 
        theme_bw() + theme(panel.grid = element_blank()) + 
        guides(fill = FALSE) + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
        geom_vline(xintercept = trimm, linetype="dashed", color = "red")
      ggsave(file.path('.', "Trimming", paste0(outpref, ".pdf")))
      plotdf <- plotdf.summary
    }
    else {
        #save(plotdf, "plotdf.rda")
        trim <- process_df(plotdf, minqual, wdw, tr = tr)
        trimm <- trim[1]
        read_lenght <- trim[2]
        ggplot(data = plotdf, aes(x = Cycle, y = Score)) + geom_tile(aes(fill = Count)) +
        scale_fill_gradient(low = "#F5F5F5", high = "black") +
        geom_line(data = statdf, aes(y = Mean), color = "#66C2A5") + 
        geom_line(data = statdf, aes(y = Q25), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        geom_line(data = statdf, aes(y = Q50), color = "#FC8D62", size = 0.25) + 
        geom_line(data = statdf, aes(y = Q75), color = "#FC8D62", size = 0.25, linetype = "dashed") + 
        ylab("Quality Score") + xlab("Cycle") + theme_bw() + 
        theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
        geom_vline(xintercept = trimm, linetype="dashed", color = "red") +
        geom_hline(yintercept = minqual, linetype="dashed", color = "blue") +
        geom_text(data = anndf, aes(x = 0, label = label, y = minScore + 2), hjust = 0, vjust = 0) + 
        geom_text(data = anndf, aes(x = 0, label = rclabel, y = minScore + 2), hjust = 0, vjust = 2) + 
        facet_wrap(~file) + theme(strip.background = element_blank(), strip.text.x = element_blank())
        ggsave(file.path('.', "Trimming", paste0(outpref, ".pdf")))
    }
  return(list(df = plotdf, trimm = trimm, read_lenght = read_lenght))}

parse_log <- function(logfile, name){
    alllines <- readLines(logfile)
    bl <- strsplit(alllines[8:12], ' ')
    input <- as.numeric(gsub(',', '', tail(bl[[1]], n=1)))
    cutadapt <- as.numeric(gsub(',', '', head(tail(bl[[5]], n=2),n=1)))
    v <- cbind(input, cutadapt)
    row.names(v) <- name
    return(v)
}

execute_cutdadapt <- function (forward, reverse, primerF, primerR, len){#}, cutting){len <- len - max(nchar(primerF), nchar(primerR))
    name <- sapply(strsplit(basename(forward), "_R"), `[`, 1)
    logfile = paste0('./Trimming/', name, '.log')
    Adapter2rc = paste(reverseComplement(DNAString(primerR)),  collapse='')
    Adapter1rc = paste(reverseComplement(DNAString(primerF)),  collapse='')
    command = paste('cutadapt -m', len, '-g', primerF, '-G', primerR, '-a',
    Adapter2rc, '-A', Adapter1rc, '-o', paste0('./Trimming/', 'trimmed_', name,
    '_R1.fastq'), '-p', paste0('./Trimming/trimmed_', name, '_R2.fastq'),
    '--match-read-wildcards', '--trim-n',  '-n 2', '--untrimmed-output',
    './Trimming/untrimmed_U.fastq', '--untrimmed-paired-output',
    './Trimming/untrimmed_P.fastq', forward, reverse, '>', logfile)
    if (!file.exists(paste0('./Trimming/', 'trimmed_', name, '_R1.fastq'))){
      print(paste('Executing cutadapt command:', command))
      system(command)}
      write(command, file=paste0(name, '_cutadapt.command.txt'))
      alllines <- readLines(logfile)
      bl <- strsplit(alllines[8:12], ' ')
      input <- as.numeric(gsub(',', '', tail(bl[[1]], n=1)))
      cutadapt <- as.numeric(gsub(',', '', head(tail(bl[[5]], n=2),n=1)))
      v <- cbind(input, cutadapt)
      row.names(v) <- name
      return(v)
}

getN <- function(x) {
  sum(getUniques(x))
}

seqkitRC <- function(filelist){
    newl <- c()
    for (i in seq(length(filelist))){
        name <- filelist[[i]]
        dirn <- dirname(name)
        basn <- basename(name)
        j <- file.path(dirn, paste0('RC_', basn))
        comm <- paste('seqkit', 'seq', '-r', '-p', name, '|', 'gzip', '-c',
        '>', j)
        system(comm)
        newl <- c(newl, j)
    }
return(newl)}


run_cutadapt <- function(fnFs, fnRs, primerF, primerR, min_read_length,
use_parallel, packages){
    exports=c('execute_cutdadapt', 'parse_log')
    if (!use_parallel){
    caout <- vector()
    for (i in seq(length(fnFs))){
        m <- execute_cutdadapt(fnFs[i], fnRs[i], primerF, primerR, min_read_length)
        caout <- rbind(caout, m)}}else{
        caout <- foreach(i=seq(length(fnFs)), .combine = rbind, .export=exports,
        .packages=packages) %dopar% execute_cutdadapt(fnFs[i], fnRs[i], primerF, primerR,
        min_read_length)}
    #save(caout, file = 'cutadapt.rda')
    return(caout)
}

### Start of the script
option_list = list(
  make_option(c("-p", "--path"), type="character", default=getwd(), 
              help=paste("Path where code should be executed (where the files",
              "are) [default= %default]")),
  make_option(c("-P", "--prefix"), type="character", default=NULL, 
              help="prefix of the fastq files [default= %default]"),
  make_option(c("-m", "--min_qual"), type="numeric", default=15, 
              help=paste("Minimum average (within window) allowed before",
              "trimming [default= %default]")),
  make_option(c("-w", "--window"), type="numeric", default=1,
              help="Window size to compute average quality [default= %default]"),
  make_option(c("-f", "--fwd_primer"), type="character", default=NULL, 
              help="Forward primer for demultiplexing and trimming [default= %default]"),
  make_option(c("-r", "--rev_primer"), type="character", default=NULL, 
              help="Reverse primer for demultiplexing and trimming [default= %default]"),
  make_option(c("-T", "--trim_type"), type="character", default='auto',
              help=paste("Type of trimming (auto, fixed and Null), if fixed",
              "pass comma separated values [default= %default]")),
  make_option(c("-d", "--db"), type="character", default=NULL,
              help="Path to database (dada formated)  [default= %default]"),
  make_option("--predict", action="store_true", help="Path to database (dada
  formated)  [default= %default]"),
  make_option(c("--reverse_complement"), action="store_true", default=FALSE,
  help=paste0("After filter and trimming convert the reads to reverse ",
    "complement. This is necesary for RNA [default= %default]")),
  make_option(c("--min_overlap"), action="store", default=20, type="numeric",
  help=paste0("Minimum overlap for merging [default= %default]")),
  make_option(c("--min_amplicon_length"), action="store", default=200, type="numeric",
  help=paste0("Minimum amplicon lenght expected [default= %default]")),
    make_option(c("--max_amplicon_length"), action="store", default=0, type="numeric",
  help=paste0("Minimum amplicon lenght expected [default= %default]")),
  make_option(c("--cpus"), action="store", default=detectCores(
    all.tests = FALSE, logical = TRUE), type="numeric", 
    help=paste0("Number of cpus to use [default= %default]")),
  make_option(c("--maxEE_fwd"), action="store", default=5, type="numeric",
              help=paste("Maximum number of 'expected errors' allowed in",
                         "a forward read [default= %default]")),
   make_option(c("--maxEE_rev"), action="store", default=5, type="numeric",
              help=paste("Maximum number of 'expected errors' allowed in",
                         "a reverse read [default= %default]")),
  make_option(c("--min_asv_size"), action="store", default=8, type="numeric",
              help=paste("Minumum size of an ASV  [default= %default]")),
  make_option(c("--priors"), action='store', default=character(0),
              help="Prior sequences. File with one headerless sequence per line")
  # make_option(c("--max_cpus"), action="store", default=-1, type="numeric",
  #             help=paste("Maximum number of cpus to use  [default= %default]"))
  )
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("THIS SCRIPT ASSUMES THAT YOUR DATA HAS BEEN MULTIPLEXED PER SAMPLE")
print(opt, sep = "\n")

cpus <- opt$cpus
use_parallel <- !(cpus == 0 | cpus == 1)
if (!use_parallel){
    cpus <- FALSE
} else if (cpus == -1){
    cpus <- TRUE
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
} else{
    cl <- makeCluster(cpus)
    registerDoParallel(cl)
}
if(length(opt$priors) == 0){priors<-opt$priors}else{
priors <- scan(file = opt$priors, what=character())
}
path <- opt$path
prefix <- opt$prefix
min_qual <- opt$min_qual
wdw <- opt$window
primerF <- opt$fwd_primer
primerR <- opt$rev_primer
trimm_type <- opt$trim_type
min_asv_size <- opt$min_asv_size
if (grepl(',', trimm_type)){
  trimm_type = strsplit(trimm_type,',')[[1]]
}
db <- opt$db
reverse <- opt$reverse_complement
min_overlap <- opt$min_overlap
min_amplicon <- opt$min_amplicon_length
max_amplicon <- opt$max_amplicon_length
min_read_length <- round( (min_amplicon/2) + min_overlap )
maxEE=c(opt$maxEE_fwd,opt$maxEE_rev) # this is relaxed and differs from tutorial

# Forward and reverse fastq filenames have format: SAMPLENAME<PATTERN>
if (is.null(prefix))
{
  fnFs <- sort(list.files(path, pattern='_R1.fastq.gz', full.names = TRUE))
  fnRs <- sort(list.files(path, pattern='_R2.fastq.gz', full.names = TRUE))
} else {
  fnFs <- sort(list.files(path, pattern=paste0("^", prefix, "_R1.fastq.gz"),
  full.names = TRUE))
  fnRs <- sort(list.files(path, pattern=paste0("^", prefix, "_R2.fastq.gz"),
  full.names = TRUE))
}

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# plot qualities before trimming and trim
tr1=vector()
tr2=vector()
dir.create('Trimming')
# Remove adapters (I am lazy so will use cutadapt)
packages <- c('tools', 'plyr', 'ShortRead', 'dada2', 'Biostrings')
if (file.exists('cutadapt.rda')) {
  load(file = 'cutadapt.rda')} else {
    caout <- run_cutadapt(fnFs, fnRs, primerF, primerR, min_read_length,
    use_parallel, packages)
    save(caout, file = 'cutadapt.rda')
    print('Head of output of after cutadapt:')
    print(head(caout))
  }
# Get the new adaptors-free filenames 
fnFs <- sort(Sys.glob("./Trimming/trimmed*R1.fastq"))
fnRs <- sort(Sys.glob("./Trimming/trimmed*R2.fastq"))


if ('fixed' %in% trimm_type){
  trimm = as.numeric(trimm_type[2:3])}else{trimm=c(NULL, NULL)}
for (i in seq(length(fnFs))){
  o1 <- quality_plot(fnFs[[i]], minqual=min_qual, wdw=wdw, tr=trimm[1])
  o2 <- quality_plot(fnRs[[i]], minqual=min_qual, wdw=wdw, tr=trimm[2])
  tr1 <- c(tr1, o1$trimm)
  tr2 <- c(tr2, o2$trimm)
  read_lenght <- o1$read_lenght
  }

if (is.null(trimm_type)) {
  print('Trimming disabled')
  trimm = 0
} else{
  if ('fixed' %in% trimm_type){
    trimm = as.numeric(trimm_type[2:3])} else if ('minoverlap' %in% trimm_type) {
      trimm = c(min_read_length, min_read_length)
    }
  else{
    trimm = c(max(min_read_length, min(tr1)), max(min_read_length, min(tr2)))
  }
  print(paste('Trimming forwards at', trimm[1], 'and reverse at', trimm[2]))
}


# Place filtered files in filtered/ subdirectory
filtFs <- file.path("./Filtering", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path("./Filtering", paste0(sample.names, "_R2_filt.fastq.gz"))



if (file.exists('out.rda')) {
  load(file = 'out.rda')} else {
out <- filterAndTrim(fnFs, filtFs, rev=fnRs, filt.rev=filtRs, truncLen=trimm,
                     maxN = 0, maxEE = maxEE, truncQ = 2, rm.phix = TRUE, minQ = 15,
                     compress = TRUE, multithread = cpus, minLen=min_read_length,
                     matchIDs=TRUE)
row.names(out) <- sample.names
save(out, file = 'out.rda')
}
print('Head of output of after trimming:')
head(out)
write.table(out, file='out.tsv', sep='\t')
if (reverse){
    filtFs <- seqkitRC(filtFs)
    filtRs <- seqkitRC(filtRs)
}
print('Learning the Error Rates')
if (file.exists('errors.rda')) { 
  load(file = 'errors.rda')} else {
errF <- learnErrors(filtFs, multithread=cpus,  randomize = TRUE)
errR <- learnErrors(filtRs, multithread=cpus,  randomize = TRUE)
save( errF, errR , file = 'errors.rda')
  }
plotErrors(errF, nominalQ=TRUE)
ggsave(filename="Errors_Forward.pdf")
plotErrors(errR, nominalQ=TRUE)
ggsave(filename="Errors_Reverse.pdf")
print('Dereplicating...')
if (file.exists('dereps.rda')) { 
  load(file = 'dereps.rda')} else {
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
save(derepFs, derepRs, file = 'dereps.rda')}
# Name the derep-class objects by the sample names if multiple
if (is.null(prefix)){
names(derepFs) <- sample.names
names(derepRs) <- sample.names}
print('Sample inference')
if (file.exists('dadas.rda')) { 
  load(file = 'dadas.rda')} else {
dadaFs <- dada(derepFs, err=errF, multithread=cpus, pool="pseudo", priors=priors)
dadaRs <- dada(derepRs, err=errR, multithread=cpus, pool="pseudo", priors=priors)
save(dadaFs, dadaRs, file = 'dadas.rda')}
print('Merging reads ...')
if (file.exists("mergers.rda")) { 
  load(file = 'mergers.rda')} else {
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,
                        minOverlap=min_overlap )
save(mergers, file = "mergers.rda")}
# Inspect the merger data.frame from the first sample
print('Info on the first sample after merge')
head(mergers[[1]])
#write.table(mergers, file='mergers.tsv', sep='\t')
print('Making the sequence table ...')
if (file.exists("seqtab.rda")) { 
  load(file = "seqtab.rda")} else {
seqtab <- makeSequenceTable(mergers)
# remove size smaller than min_asv_size
col.sums <- apply(seqtab, 2, sum)
seqtab <- seqtab[,col.sums >= min_asv_size]

if (is.null(prefix)){
    fnc <- colnames
    if (max_amplicon != 0){
    #  “cutting a band” in-silico
    seqtab <- seqtab[,nchar(fnc(seqtab)) %in% min_amplicon:max_amplicon]
}
}else{
    fnc <- names
    if (max_amplicon != 0){
    #  “cutting a band” in-silico
    seqtab <- seqtab[nchar(fnc(seqtab)) %in% min_amplicon:max_amplicon]
}
}
save(seqtab, file = "seqtab.rda")}
write.table(seqtab, file='seqtab.tsv', sep='\t')
print('Distribution of sequence lengths')
table(nchar(getSequences(seqtab)))
print('Removing Chimeras')
if (file.exists("nochimaeras.rda")) { 
  load(file = "nochimaeras.rda")} else {
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=cpus, verbose=TRUE)
all <- getUniques(seqtab.nochim)
uniquesToFasta(all, 'all_samples.fasta')
seqnames <- get.fasta.name('all_samples.fasta', clean_name = FALSE)
otutab <- seqtab.nochim
if (is.null(prefix)){
    colnames(otutab) <- seqnames
}else{names(otutab) <- seqnames}
save(otutab, seqtab.nochim, file = "nochimaeras.rda")}
print(paste('Proportion of non-chimaeric sequences', sum(seqtab.nochim)/sum(seqtab)))
write.table(otutab, file='ASV_table.tsv', sep='\t')
if (length(filtFs) > 1){
    track <- cbind(caout, out[, 2], sapply(dadaFs, getN), sapply(dadaRs, getN),
    sapply(mergers, getN), rowSums(otutab))
}else{
    track <- cbind(caout, out[, 2], getN(dadaFs), getN(dadaRs), getN(mergers),
    sum(otutab))
}
colnames(track) <- c('raw', 'cutadapt', "filtered", "denoisedF", "denoisedR",
"merged", "nonchim")
rownames(track) <- sample.names
print('Head of the summary of the reads through the pipeline')
head(track)
write.table(track, file = "allsamples_track.tsv", sep='\t')
dir.create('ASVs')
for(sample in row.names(seqtab.nochim)){
    st <- seqtab.nochim[sample,]
    boolean <- st > 0
    st <- st[boolean]
    seqs <- DNAStringSet(getSequences(st))
    names(seqs) <- seqnames[boolean]
    outfn <- file.path("./ASVs", paste0(sample,'_ASV.fasta'))
    writeXStringSet(seqs, outfn)
}


# #with dada
# taxa <- assignTaxonomy(seqtab.nochim, db, multithread=TRUE)
# #with decipher
# if(grepl('Rdata', db)){load(db)}else{
# if(grepl('Rdata', db)){load(db)}else{
# train <- readDNAStringSet(db)
# s <- strsplit(names(train), ";") # get names (has to be formated properly)
# domain <- sapply(s, '[', 1)
# phylum <- sapply(s, '[', 2)
# class <- sapply(s, '[', 3)
# order <- sapply(s, '[', 4)
# family <- sapply(s, '[', 5)
# genus <- sapply(s, '[', 6)
# species <- sapply(s, '[', 7)
# taxonomy <- paste("Root", domain, phylum, class, order, family, genus, species,
# sep="; ")
# trainingSet <- LearnTaxa(train, taxonomy)
# plot(trainingSet)}
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#         m <- match(ranks, x$rank)
#         taxa <- x$taxon[m]
#         taxa[startsWith(taxa, "unclassified_")] <- NA
#         taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
# write.table(taxa, file='naive_bayes.tsv', sep='\t')
# write.table(taxid, file='decipher.tsv', sep='\t')
#

