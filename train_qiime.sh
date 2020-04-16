#!/usr/bin/env bash
fasta=$1
prefix=$(echo "${fasta}"| awk 'BEGIN{FS=OFS="."} NF--')
taxfile=$2
#fwd=$3 #GGWACWGGWTGAACWGTWTAYCCYCC
#rvs=$4 #TAAACTTCAGGGTGACCAAAAAATCA
echo "Processing ${prefix} "
seqkit seq -u ${fasta} > tmp

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path tmp \
  --output-path ${prefix}_feature.qza

rm tmp

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ${taxfile} \
  --output-path ${prefix}_ref-taxonomy.qza

#qiime feature-classifier extract-reads \
#  --i-sequences ${prefix}_feature.qza \
#  --p-f-primer ${fwd}  \
#  --p-r-primer ${rvs} \
#  --p-min-length 100 \
#  --p-max-length 600 \
#  --o-reads ${prefix}_ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  ${prefix}_feature.qza \
  --i-reference-taxonomy ${prefix}_ref-taxonomy.qza \
  --o-classifier ${prefix}_classifier.qza