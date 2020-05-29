# TA_pipes
Scripts and pipelines for the manuscript "Assessment of current taxonomic
 assignment strategies for metabarcoding eukaryotes: Insights from mock 
 communities". If you use any of this code in your reserch, please cite
 this repository:
> Hleap, 2020. TA_pipes. https://github.com/jshleap/TA_pipes.git

## Obtaining all the code
If you have git installed, you can get all the scripts by:
```bash
git clone https://github.com/jshleap/TA_pipes.git
```
This will generate a folder called TA_pipes where all the scripts are

### Installing dependencies
To install all dependencies required to run any part of the benchmarking:

1. Move to the `TA_pipes` directory by:
```bash
cd TA_pipes
```
2. Execute the `install_all_dependencies.sh`:
```bash
bash install_all_dependencies.sh
```

That is all! All R and python dependencies will be installed. The
`install_all_dependencies.sh` script assumes you have R and python3 already 
installed.

## dada2_pipe.R
dada2_pipe.R is the main R script for the pre-processing of the mock 
communities, and first pass for quality control. It has an automatic check
on when the quality of the reads fall below your desired median quality.
It makes use of `cutadapt` to remove adapters and to pre-trim/filter the\
 sequences.
 
### Usage of dada2_pipe.R
To execute this script you need to use the `Rscript` command. To check the
help of `dada2_pipe.R` you can do:

```bash
Rscript dada2_pipe.R -h
```
You should see something like this:

```
[1] "running on DADA2 version 1.14.1 and DECIPHER version 2.14.0"
Usage: dada2_pipe.R [options]


Options:
	-p PATH, --path=PATH
		Path where code should be executed (where the files are) [default= /home/jshleap/my_gits/TA_pipes]

	-P PREFIX, --prefix=PREFIX
		prefix of the fastq files [default= NULL]

	-m MIN_QUAL, --min_qual=MIN_QUAL
		Minimum average (within window) allowed before trimming [default= 15]

	-w WINDOW, --window=WINDOW
		Window size to compute average quality [default= 1]

	-f FWD_PRIMER, --fwd_primer=FWD_PRIMER
		Forward primer for demultiplexing and trimming [default= NULL]

	-r REV_PRIMER, --rev_primer=REV_PRIMER
		Reverse primer for demultiplexing and trimming [default= NULL]

	-T TRIM_TYPE, --trim_type=TRIM_TYPE
		Type of trimming (auto, fixed and Null), if fixed pass comma separated values [default= auto]

	--reverse_complement
		After filter and trimming convert the reads to reverse complement. This is necesary for RNA [default= FALSE]

	--min_overlap=MIN_OVERLAP
		Minimum overlap for merging [default= 20]

	--min_amplicon_length=MIN_AMPLICON_LENGTH
		Minimum amplicon lenght expected [default= 200]

	--max_amplicon_length=MAX_AMPLICON_LENGTH
		Minimum amplicon lenght expected [default= 0]

	--cpus=CPUS
		Number of cpus to use [default= 12]

	--maxEE_fwd=MAXEE_FWD
		Maximum number of 'expected errors' allowed in a forward read [default= 5]

	--maxEE_rev=MAXEE_REV
		Maximum number of 'expected errors' allowed in a reverse read [default= 5]

	--min_asv_size=MIN_ASV_SIZE
		Minumum size of an ASV  [default= 8]

	--priors=PRIORS
		Prior sequences. File with one headerless sequence per line

	-h, --help
		Show this help message and exit

```

Let's imagine you have a set of amplicon sequences of the Leray primers ( COI; 
Fwd: GGWACWGGWTGAACWGTWTAYCCYCC; Rev: TAAACTTCAGGGTGACCAAAAAATCA), from 10 
samples (and hence 10 files). Say that the  files are in the following path
`/home/user/Lake_samples`. Also, lets assume you want to retain reads that have
 a minimum quality of 18 in a window of 1 (as soon as it falls below 18, the 
 read is trimmed), and that this is done automatically. We also have some prior 
 knowledge of our amplicon, and we know that it ranges between 200bp and 400bp,
  and we want to avoid shorter or larger fragments. We also want to run it only
   on one cpu (sometimes the  multiprocessor fails). You want to run the dada2
   denoising with maximum of 4 errors in the forward reads and 5 in the reverse.
 You also want to consider only fragments that can be merged with an overlap of
 25bp. 
 

To execute the above, you can type:

```bash
Rscript dadatest.R \
    --path=/home/user/Lake_samples
    --fwd_primer=GGWACWGGWTGAACWGTWTAYCCYCC \
    --rev_primer=TAAACTTCAGGGTGACCAAAAAATCA \ 
    --min_qual=18 \
    --window=1 \    
    --max_amplicon_length=400 \
    --min_amplicon_length=200 \
    --cpus=1 \
    --maxEE_fwd=4 \
    --maxEE_rev=5
```

## format_all_references.sh
This script assumes that you have access to the reference sets of the
manuscript (DRYAD DOI will be included soon). It also assumes that you
have access to the PROTAX script provided by Panu Somervuo for conserved
sequences. You will need to write Panu Directly as we do not have permission
to share his code.
Now, assuming you have all the references in the path `/home/user/references`,
that the Protax scripts are in the path `/home/user/PROTAX/scripts`, and
that this git can be found at `/home/user/TA_pipes` you can execute the 
formatting of the reference datasets by:

```bash
cd /home/user/references 
source /home/user/TA_pipes/format_all_references.sh /home/user/PROTAX/scripts
```

## process_mock.py
This script is intended to do the first pass (before manual curation) on
the mock communities. It assumes that the mock communities have been 
preprocessed and that ASV/ZOTUS are reported. It also requires a file with
the list of intended species in the mock community to guide the curation.
It makes use of blast and blast-formatted databases and tree reconstruction.
It also uses the NCBI's Entrez service, and therefore internet connection
is required.

You can check the help by:

```bash
python3 process_mock.py -h
```

You should see something like this:

```bash
usage: process_mock.py [-h] [-d DB] [-b DBNAME] [-e EVALUE] [-g GAPS]
                       [-c CPUS] [-n NTOP] [-m MARKER] [-s SIG_LEVEL]
                       [-p PREFINQUERY]
                       email query intended

positional arguments:
  email                 email for Entrez queries
  query                 Fasta file with query sequences
  intended              file with intended list of species

optional arguments:
  -h, --help            show this help message and exit
  -d DB, --db DB        Path to reference database to use
  -b DBNAME, --dbname DBNAME
                        name of reference database
  -e EVALUE, --evalue EVALUE
                        Evalue threshold for Blast
  -g GAPS, --gaps GAPS  Requested fraction of ungapped columns in alignment
  -c CPUS, --cpus CPUS  Number of cpus to use
  -n NTOP, --ntop NTOP  Number of top hits to retain
  -m MARKER, --marker MARKER
                        NCBI-compliant gene name
  -s SIG_LEVEL, --sig_level SIG_LEVEL
                        Significance level on outlier test
  -p PREFINQUERY, --prefinquery PREFINQUERY
                        prefix of all sequences of interest in query file
```

## benchmark_run.py
This script is the executable of the benchmark. It does not
require any input. However, you **HAVE** to modify the `databases.py`
file to represent your system, especially:
1. base: This represents the path to your working directory, including where
the mock communities and references are.
2. path2basta_exe: This variable points to where the basta binary is
3. path2protax_scripts: Path to the conserved version of protax scripts

You will have to adjust all database names to match to your own formatted
databases, and those databases should be found in the path defined in the`base`
variable. For the  methods that require global alignments, we used 
A2G<sup>2</sup> (https://github.com/jshleap/A2G), and therefore the gobal
and local gene references need to be provided in the databases.py. By default
it will use the ones provided in this git.