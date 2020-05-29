#!/usr/bin/env bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
PROTAX="${1}"
clean() {
  if ls -U *.shelve 1>/dev/null 2>&1; then
    rm *.shelve
  fi
}

format_it() {
  pyndada=${DIR}
  prefix=${1%%.withoutliers}
  krak=$(${prefix//\./_})
  PROTAX="${2}"
  priors=10,5,1
  # fix most names
  python3 "${pyndada}"/fix_taxonomy.py "${1}" >tmp && mv tmp "${1}"
  # format kraken
  if ls -U "${krak}"/*.k2d 1>/dev/null 2>&1; then
    echo "Kraken already done"
  else
    python3 "${pyndada}"/format_kraken_reference.py ${1}
  fi
  # Format and train hmmufotu
  clean
  if [[ ! -f "${prefix}".tree ]]; then
    fasttree -nt -gtr < "${1}" >"${prefix}".tree
  fi
  if [[ ! -f "${prefix}"_hmmufotu_GTR_dG.csfm ]]; then
    bash ${pyndada}/train_hmmufotu.sh "${prefix}".tree "${1}"
  fi
  # Train idtaxa
  if [[ ! -f "${prefix}"_idtaxa.rda ]]; then
    Rscript ${pyndada}/train_idtaxa.R -f ${1} -p "${prefix}"
  fi
  # Train qiime
  if [[ ! -f "${prefix}"_classifier.qza ]]; then
    python3 ${pyndada}/lineage2qiime.py "${1}" "${prefix}"
    bash /home/jshleap/my_gits/pyndada/train_qiime.sh "${1}" "${prefix}".tax
  fi
  # make blastDB
  if [[ ! -f "${prefix}".nog ]]; then
    makeblastdb -dbtype nucl -in "${1}" -out "${prefix}" -title "${prefix} database" \
      -parse_seqids
  fi
    # format protax and train protax
  clean
  if ls -U ${krak}/*_protax_intermediateFiles.tar.bz2 1>/dev/null 2>&1; then
    echo "Protax already done"
  else
    python3 ${pyndada}/format_protax.py ${1} 4
    bash ${pyndada}/train_protax.sh "${prefix}"_protax.taxonomy \
      "${prefix}"_protax.seqid2tax "${prefix}"_protax.aln ${PROTAX} ${priors}
  fi
}

export -f format_it
# align them all
for mock in FSIS SIS MIS PSS; do
  if [[ ! -f ${mock}_realized.reference_aligned.withoutliers ]]; then
    seqkit grep -v -n -r -p ';;' ${mock}_realized.reference.fasta >tmp
    mv tmp ${mock}_realized.reference.fasta
    python3 ~/my_gits/A2G/A2G/align2consensus.py midori50.consensus \
      Leray_consensus.fst ${mock}_realized.reference.fasta
  else
    echo "Alignment of ${mock} previously done ... skipping"
  fi
done
# process ISIS
if [[ ! -f ${mock}_realized.reference_aligned.withoutliers ]]; then
  python3 ~/my_gits/A2G/A2G/align2consensus.py midori50.consensus \
    Mlep_consensus.fst ISIS_realized.reference.fasta
else
  echo "Alignment of ISIS previously done ... skipping"
fi
# Format the resulting alignments
for aln in *.withoutliers; do
  format_it "${aln}" "${PROTAX}"
done
