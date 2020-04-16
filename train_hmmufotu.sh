#!/usr/bin/env bash
tree=${1}
fasta=${2}
prefix=$(echo "${tree}" | awk 'BEGIN{FS=OFS="."} NF--')

get_root(){
python3 - << EOF
import pandas as pd
df = pd.read_csv('$1', sep='\t', header=None)
tax = df[1].str.split('; ', expand=True)
nun = tax.nunique()
print(tax[max(nun[nun == 1].index)].unique()[0].split('__')[1] )
EOF
}

fix_tree(){
python3 - << EOF
from ete3 import Tree
import sys
from tqdm import tqdm
fn = "${1}"
print('Processing', fn)
prefix = fn[:fn.rfind('.')]
t = Tree(fn)
print('Tree length before:', len(t))
t.resolve_polytomy(default_dist=1E-10)
for n in tqdm(t, total=len(t)):
    n.dist = n.dist + 1E-10
with open('%s.hmmufotu.tree' % prefix, 'w') as out:
    out.write(t.write())
print('Tree length after:', len(t))
EOF
}

grep '>' "${fasta}" | sed -e "s/>//g; s/ /\t/1" -e $'s/\t/\tk__/' \
-e 's/;/; p__/1' -e 's/;/; c__/2' -e 's/;/; o__/3'  -e 's/;/; f__/4' \
-e 's/;/; g__/5' -e 's/;/; s__/6' -e 's/ /_/7g' -e 's/;$//g' > ${prefix}.hmmufotu.tax

fix_tree "${tree}"
root=$(get_root ${prefix}.hmmufotu.tax)
hmmufotu-build "${fasta}" "${prefix}".hmmufotu.tree -r ${root} \
-a ${prefix}.hmmufotu.tax -n ${prefix}_hmmufotu -s GTR -V -v --fmt fasta