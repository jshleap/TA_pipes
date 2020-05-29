# this file is read by optimize_n_score to set the databases
import os
git_path = os.path.dirname(os.path.realpath(__file__))
base = 'put your base path here'
blastdb = base + '%s_realized.reference_aligned'
qiime_class = base + '%s_realized.reference_aligned_classifier.qza'
krakendb = base + '%s_realized_reference_aligned'
path2basta_exe = 'put the path to basta here'
hmmufotu_db = base + '%s_realized.reference_aligned_hmmufotu_GTR_dG'
idtaxa_db = base + '%s_realized.reference_aligned_idtaxa.rda'
basta_map = 'benchmark'
path2protax_scripts = 'path to the conserved version of protax scripts'
protaxdb = '%s_realized.reference_aligned_protax_'
protax_ref = base + '%s_realized_reference_aligned'
protax_global = os.path.join(git_path, 'a2g_ref', 'midori50.consensus')
protax_local = os.path.join(git_path, 'a2g_ref', 'Leray_consensus.fst')
protax_local_isis = os.path.join(git_path, 'a2g_ref', 'Mlep_consensus.fst')
