#!/usr/bin/env bash
# train_protax.sh --- Train and test protax taxonomic classifier
# Author: Jose Sergio Hleap <jose.hleaplozano@mcgill.ca>
# Created: 2019-08-28

# Usage: bash train_protax.sh taxfile seqid2taxfile alnfile pathtoprotaxscripts priorstr mintax
set -e


taxfile=$1
s2tfile=$2
alnfile=$3
PROTAX=$4
priors=$5
prefix=$(basename ${s2tfile%%_protax.seqid2tax})
d=$(echo ${prefix}| tr '.' '_')
out=${d}/${prefix}_protax_
mkdir -p ${d}
#cd ${prefix}
train_n=$(grep -c '>' "${alnfile}")
# 1) set variable NUM_TAXLEVELS based on your taxonomy,
NUM_TAXLEVELS=$(cut -f 2 "${s2tfile}"| uniq)



# Functions
execute_step5(){
Rscript - << EOF
source("${1}/amcmc.rcode_noweight.txt")
num_tax <- ${2}
readapt <- function(pp1, stepstr){
initstate=initialize.adaptation(pp1[['params']][2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title=paste(stepstr, "readapted"))
return(pp1)}

run_adaptmcmc <- function(pp1, train, step, ind, num.params){
print(paste("Processing", train))
l <- paste0("L", step)
dat <- read.xdata(train)
pp1 <- adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title=l)
# if problems with mixing, usually one or two re-adaptations are enough, run
#pp1 <- readapt(pp1, L1)
k <- which.max(pp1[['postli']][ind])
write.postparams(pp1, paste0("${out}", "mcmc", step), ind[k])
return(pp1)
}

library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)
num.params=1+4
ind=1001:2000
fn <- file.path(dirname("${out}"), basename("${out}"))
dev.new(width=10, height=5, file=paste0(fn, ".pdf"))
# get scxdat files
files <- Sys.glob("${out}*.scxfdat")
pp1 <- NULL
for (i in seq(length(files))){pp1 <- run_adaptmcmc(pp1, files[[i]], i, ind, num.params)}
EOF
}


#TODO: solve levels
bap(){
Rscript - << EOF
source("${1}/amcmc.rcode_noweight.txt")
nimi <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
num_tax <- $2
nimi <- rev(rev(nimi)[1:num_tax])
pdf(file="${out}accuracy.pdf")
par(mfrow=c(2,2))
for (i in 1:num_tax) {
 file=sprintf("${out}query%d.cor",i)
 a=read.table(file,header=F)
 accuracy.plot(a[,3],a[,4],name=sprintf("Model performance: %s",nimi[i]))
}
dev.off()
EOF
}


# 2) add priors for unknown taxa in taxonomy, one value for each level
#    this relates to how much you think there are taxa not included in your
#    taxonomy, larger values add more uncertainty to all predictions unk prior
#    is level-specific where units correspond to leaf nodes of the taxonomy
#    (value * prior(known_species)) NOTE: the number of priors in the ,,, list
#    must equal to $NUM_TAXLEVELS
echo "Setting priors and thinning taxonomy"
if [[ ! -s ${out}taxonomy.priors ]]; then
perl "${PROTAX}"/taxonomy_priors.pl "${priors}" "${taxfile}" > ${out}taxonomy.priors
fi

for ((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
do
  if [[ ! -s  ${out}tax"${LEVEL}" ]]; then
    perl "${PROTAX}"/thintaxonomy.pl "${LEVEL}" ${out}taxonomy.priors > ${out}tax"${LEVEL}"
  fi
done

# 3) generate training data for each level, here 10000 training samples per
# level
echo "Generating training data"
for ((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL ${LEVEL}"
 if [[ ! -s ${out}train"${LEVEL}".numeric ]]; then
   perl "${PROTAX}"/seqid2taxlevel.pl "${LEVEL}" "${s2tfile}" > ${out}ref.tax"${LEVEL}"
   perl "${PROTAX}"/get_all_reference_sequences.pl "${LEVEL}" ${out}tax"${LEVEL}" \
   ${out}ref.tax"${LEVEL}" ${out}rseqs"${LEVEL}"
   perl "${PROTAX}"/taxrseq2numeric.pl "${LEVEL}" ${out}tax"${LEVEL}" \
   "${alnfile}" > ${out}rseqs"${LEVEL}".numeric
   perl "${PROTAX}"/generate_training_data.pl "${LEVEL}" ${out}tax"${LEVEL}" \
   ${out}ref.tax"${LEVEL}" ${out}rseqs"${LEVEL}" "${train_n}" 1 no ${out}train"${LEVEL}"
   perl "${PROTAX}"/traindat2numeric.pl "${alnfile}" \
   ${out}train"${LEVEL}" > ${out}train"${LEVEL}".numeric
   # to check what kind of training data there is for each level:
   cut -f4 -d" " ${out}train"${LEVEL}" | cut -f1 -d"," | sort | uniq -c >${out}_info.${LEVEL}
 fi
done

# 4) calculate xdat file (sequence similarity predictors), scale the values and
# save the scaling parameters for later use ...this can take a while if large
# training data...
echo "Generating sequence similarity predictors"
for ((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL ${LEVEL}"
 if [[ ! -s ${out}train${LEVEL}.xdat ]]; then
   "${PROTAX}"/create_xdata_best2 ${out}tax"${LEVEL}" "${alnfile}" \
   ${out}rseqs"${LEVEL}".numeric ${out}train"${LEVEL}".numeric > ${out}train"${LEVEL}".xdat
 fi
done

if [[ ! -s  ${out}train"${LEVEL}".scxfdat ]]
then
  for ((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
  do
     perl "${PROTAX}"/scale_xdat.pl ${out}sc"${LEVEL}" ${out}train"${LEVEL}".xdat > ${out}train"${LEVEL}".scxfdat
  done
fi


# 5) parameter estimation in R
#    MCMC for parameters separately in each taxonomy level you need to check
#    the convergence and continue iterations or re-initialize adaptive proposal
#    if needed
echo "Estimating parameters"
if [[ ! -s ${out}query0.prob.numeric ]]
then
  execute_step5 "${PROTAX}" "${NUM_TAXLEVELS}"
fi

##################################
# Check model with training data
##################################

# Check classification with training samples with the lowest level
echo "Check classification with training samples"
#train=$(ls train*| sort -V| head -1)
if [[ ! -s ${out}par"${LEVEL}" ]]
then
    train=${out}train"${NUM_TAXLEVELS}"
    perl "${PROTAX}"/init_train_taxprob2numeric.pl "${train}" > ${out}query0.prob.numeric
    for i in ${out}mcmc*
    do
     LEVEL=$(grep -o "[0-9]" <<< ${i})
     cut -f3-6 -d" " "${i}" > ${out}par"${LEVEL}"
    done
fi

# parent probs from previous level classification
for((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
do
 echo "LEVEL ${LEVEL}"
 PREVLEVEL=$((LEVEL-1))
 IFILE=${out}query${PREVLEVEL}.prob.numeric
 OFILE=${out}query${LEVEL}.prob
 # add scaling
 "${PROTAX}"/trainclassify_best2 ${out}tax"${LEVEL}" "${alnfile}" \
 ${out}rseqs"${LEVEL}".numeric ${out}par"${LEVEL}" ${out}sc"${LEVEL}" 0.01 \
 "${train}".numeric "${IFILE}" > "${OFILE}"
 perl "${PROTAX}"/train_taxprob2numeric.pl ${out}tax"${LEVEL}" "${OFILE}" > "${OFILE}".numeric
 perl $PROTAX/trainsample2correct.pl $LEVEL ${taxfile} ${train} > ${out}query${LEVEL}.tax
 perl $PROTAX/trainsample2addcor.pl ${out}query${LEVEL}.prob ${out}query${LEVEL}.tax > ${out}query${LEVEL}.cor
done

##### calculate correctness (note: correct label is the one which training sample mimicked, it can be e.g. unknown species)
#for ((LEVEL=1; LEVEL<=NUM_TAXLEVELS; LEVEL++))
#do
# perl $PROTAX/trainsample2correct.pl $LEVEL ${taxfile} ${train} > ${out}query${LEVEL}.tax
# perl $PROTAX/trainsample2addcor.pl ${out}query${LEVEL}.prob ${out}query${LEVEL}.tax > ${out}query${LEVEL}.cor
#done



echo "Creating the bias accuracy plots"
##### bias accuracy plots
bap "${PROTAX}" "${NUM_TAXLEVELS}"

#### preparing for classification
for ((LEVEL=1; LEVEL<=$NUM_TAXLEVELS; LEVEL++))
do
 cut -f3-6 -d" " ${out}mcmc$LEVEL >> ${out}model.pars
 cat ${out}sc$LEVEL >> ${out}model.scs
 cat ${out}rseqs${LEVEL}.numeric >> ${out}model.rseqs.numeric
done

### cleanup
find ${d} -type f  \( -iname "${prefix}*" ! -iname "*model*" ! \
-iname "*taxonomy.priors" \) -print0 |tar -cjf ${out}intermediateFiles.tar.bz2 \
--null -T -

find ${d} -type f  \( -iname "${prefix}*" ! -iname "*model*" ! \
-iname "*taxonomy.priors" \) -exec rm {} +
