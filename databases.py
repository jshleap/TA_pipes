# this file is read by optimize_n_score to set the databases
import os
base = '/home/jshleap/Playground/Taxonomic_assignment/classdb/'
blastdb = base + '%s_realized.reference_aligned'
qiime_class = base + '%s_realized.reference_aligned_classifier.qza'
krakendb = base + '%s_realized_reference_aligned'
path2basta_exe = '/home/jshleap/Programs/BASTA/bin/basta'
hmmufotu_db = base + '%s_realized.reference_aligned_hmmufotu_GTR_dG'
idtaxa_db = base + '%s_realized.reference_aligned_idtaxa.rda'
basta_map = 'benchmark'
path2protax_scripts = '/home/jshleap/Programs/protaxA/scripts'
protaxdb = '%s_realized.reference_aligned_protax_'
protax_ref = base + '%s_realized_reference_aligned'
protax_global = os.path.join(base, 'midori50.consensus')
protax_local = os.path.join(base, 'Leray_consensus.fst')
protax_local_isis = os.path.join(base,'Mlep_consensus.fst')
