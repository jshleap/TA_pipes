# Title     : TODO
# Objective : TODO
# Created by: jshleap
# Created on: 2019-08-21

suppressPackageStartupMessages(library('optparse'))
suppressPackageStartupMessages(library('parallel'))
suppressPackageStartupMessages(library('foreach'))
suppressPackageStartupMessages(library('doParallel'))
suppressPackageStartupMessages(library('Biostrings'))
suppressPackageStartupMessages(library('DECIPHER'))

fnc <- function (id) {paste(id$taxon,sep=";", collapse=";")}

### Start of the script
option_list = list(
    make_option(c("-t", "--trainingset"), type="character", default=NULL,
    help="Name of training set rds file [default= %default]"),
    make_option(c("-p", "--processors"), type="double", default=1,
    help="Number of cpus to use [default= %default]"),
    make_option(c("-q", "--queryfile"), type="character", default=NULL,
    help="Fasta file with query sequeces [default= %default]"),
    make_option(c("-T", "--threshold"), type="double", default=0.7,
    help="Confidence at which to truncate the output taxonomic classifications
     [default= %default]"),
    make_option(c("-b", "--bootstraps"), type="double", default=100, help="
    number of bootstrap replicates to perform for each sequence"),
    make_option(c("-m", "--minDescend"), type="double", default=0.95,
    help="Minimum fraction of bootstraps required to descend the tree")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#db <- readRDS(opt$trainingset)
load(opt$trainingset)
db <- trainingSet
dna <- readDNAStringSet(opt$queryfile)
ids <- IdTaxa(dna, db, type="extended", processors=opt$processors,
threshold=opt$threshold, bootstraps=opt$bootstraps, minDescend=opt$minDescend)
out <- data.frame(seq=names(ids), taxa=gsub("Root;", "",  sapply(ids, fnc)))
fn <- file.path(dirname(opt$queryfile), 'idtaxa.tsv')
#paste0(tools::file_path_sans_ext(basename(opt$trainingset)),'.tsv')
write.table(out, file=fn,  sep='\t',  row.names =FALSE)

