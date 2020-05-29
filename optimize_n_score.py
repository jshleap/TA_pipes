"""
**Copyright (C) 2019  Jose Sergio Hleap**

This script takes a curated mock community, its taxonomic mapping file to test
taxonomic assignment
"""
import shutil
import tempfile
import time
from collections import namedtuple
from io import BytesIO
from itertools import product
from multiprocessing import cpu_count
from subprocess import check_output
import pandas as pd
import dask
import dill
from dask.diagnostics import ProgressBar, ResourceProfiler
from dask_ml.model_selection import GridSearchCV
from ete3 import NCBITaxa
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import matthews_corrcoef, make_scorer
from sklearn.model_selection import KFold
from sklearn.utils.multiclass import unique_labels
import numpy as np
import argparse
from A2G.align2consensus import main as align
# from align2consensus import *
# from align2consensus import main as align
from databases import *
from subprocess import run, PIPE, CalledProcessError
from tqdm import tqdm

ncbi = NCBITaxa()
myenv = {'PATH': os.environ['PATH'], 'PYTHONPATH': os.environ["PYTHONPATH"]}
taxlevels = dict(
    zip(range(7), ['Kingdom', 'Phylum', 'Class', 'Order', 'Family',
                   'Genus', 'Species']))
stop_words = ['sp.', 'sp', 'spp', 'aff', 'aff.', 'group']
levels = ['superkingdom', 'phylum', 'order', 'class', 'family', 'genus',
              'species']

def get_lineage(name):
    """
    Function to get the six level lineage of a taxonomic name
    :param name: Taxonomic name
    """
    sp = ''
    name = ' '.join([x for x in name.split() if x not in stop_words][:2])
    try:
        taxid = ncbi.get_name_translator([name])[name][0]
    except KeyError:
        if name == 'unk':
            return ';' * len(levels)
        g, s = name.split()
        if len(name.split()) > 1:
            taxid = ncbi.get_name_translator([g])[g][0]
        else:
            taxid = ncbi.get_fuzzy_name_translation(g)[0]
        sp = '%s %s' % (g, s)
    taxid_lineage = {v: k for k, v in ncbi.get_rank(
        ncbi.get_lineage(taxid)).items() if v != 'no rank'}
    lineage = ncbi.translate_to_names(
        [taxid_lineage[x] for x in levels if x in taxid_lineage])
    lineage.append(sp)
    return ';'.join(lineage)


def fasta2array(fasta, truth=None):
    """
    Helper function to transform the fasta into an array to use the gridsearch
    :param fasta:
    :return:
    """
    l = {}
    with open(fasta) as fas:
        for i, line in enumerate(fas):
            if line.startswith('>'):
                if i > 0:
                    l[name] = list(seq)
                name = line.strip()[1:]
                seq = ''
            else:
                seq += line.strip()
    X = pd.DataFrame.from_dict(l, orient='index')
    if truth is None:
        return X
    d = {}
    truth = pd.read_csv(truth, sep='\t', header=None, names=['seq', 'taxa'])
    lin = truth.taxa.str.strip(';').str.split(';', expand=True).rename(
        columns=taxlevels)
    truth = pd.concat([truth, lin], axis=1).copy()
    truth.loc[truth.seq.str.contains('shuffled'), taxlevels.values()
    ] = 'shuffled'
    y = truth[truth.seq.isin(X.index.tolist())]
    assert y.shape[0] == X.shape[0]
    return X, y, truth


def array2fasta(array):
    """
    Helper function to translate a pandas array with the sequences into a
    fasta string
    :param array:
    :return:
    """
    fas = []
    for tup in array.itertuples():
        header = tup[0]
        seq = ''.join(x for x in tup[1:] if x is not None)
        fas.append('>%s\n%s' % (header, seq))
    return '\n'.join(fas)


class Estimator(BaseEstimator, ClassifierMixin):
    """
    Custom estimator/classifier to 'fit' the parameters to be optimized in each
    of the programs

    Parameters
    ----------

    """
    prefix = None
    taxlevel = None
    params = {}

    def __init__(self, program, db, **kwargs):
        self.program = program
        self.db = db
        for key, value in kwargs.items():
            if key not in ['db', 'program']:
                setattr(self, key, value)
                self.params[key] = value

    def fit(self, X, y):
        if X is None or y is None:
            print(X, y, self.__dict__)
        assert X.shape[0] == len(y)
        # print(type_of_target(y))
        try:
            assert y.ndim == 1
        except AssertionError:
            print(y)
            raise
        assert X.shape[0] == y.shape[0]
        self.map = dict(zip(X.index, y))
        self.y_ = y.tolist()
        self.classes_ = unique_labels(self.y_)
        self.X_ = X
        self.is_fitted_ = True
        return self

    def predict(self, X):
        idx = pd.DataFrame(X.index.tolist(), columns=['seq'])
        query = array2fasta(X)
        result = self.program(prefix=self.prefix, query=query, db=self.db,
                              **self.params)
        self.result = result
        taxa = taxlevels[self.taxlevel]
        try:
            # y = [result[result.seq == i][taxa].str.replace('_', ' ').values for
            #      i in idx]
            # y = result[result.seq.isin(idx)][taxa].str.replace('_', ' ')
            try:
                y = idx.merge(result, on='seq', how='left')[taxa].str.replace(
                    '_', ' ')
            except KeyError:
                print(result.columns)
                print(idx.columns)
                raise
            except AttributeError:
                print(idx.merge(result, on='seq', how='left').head())
                raise
            y[pd.isna(y)] = 'Not predicted'
            y = y.values
            try:
                assert y.ndim == 1
            except AssertionError:
                print(y)
                raise
        except KeyError:
            print(result)
            raise
        try:
            assert y.shape[0] == X.shape[0]
        except AssertionError:
            with open('failed.pickle', 'wb') as f:
                dill.dump((taxa, result, X, y), f)
            raise
        return y.tolist()

    def get_params(self, deep=None):
        return self.__dict__

    def set_params(self, **params):
        if not params:
            # Simple optimization to gain speed (inspect is slow)
            return self
        self.params = {}
        for key, value in params.items():
            if key not in ['db', 'program']:
                setattr(self, key, value)
                self.params[key] = value
        return self


mathews_scorer = make_scorer(matthews_corrcoef, greater_is_better=True)


def score(prediction, truth, timed, taxa):
    """
    Function to score a prediction given the true labels
    :param prediction: Vector of predicted labels
    :param truth: Vector of True labels
    :param timed: Time to run the prediction
    :param taxa: Taxonomic level
    """
    np.seterr(all='raise')
    tr = '%s_truth' % taxa
    pr = '%s_predicted' % taxa
    truth.copy().loc[truth.seq.str.contains('shuffled'), taxa] = 'shuffled'
    mer = truth.merge(prediction, on='seq', how='outer', suffixes=(
        '_truth', '_predicted'))
    all_calls = mer.shape[0]
    nshuf = sum(truth[taxa].str.contains('shuffled'))
    shuffled = mer.seq.str.contains('-shuffled')
    assert nshuf == sum(shuffled)
    try:
        assert not set(truth.seq).difference(prediction.seq).difference(
            mer.loc[shuffled, 'seq'])
    except AssertionError:
        with open('set_test.pckl', 'wb') as f:
            dill.dump((prediction, truth, mer), f)
    tp = mer[~shuffled & (mer[tr] == mer[pr])].shape[0]
    fp = mer[(~shuffled & (mer[tr] != mer[pr]) & (~mer[pr].isnull())) | (
            shuffled & ~mer[pr].isnull())].shape[0]
    fn = mer[~shuffled & (mer[tr] != mer[pr]) & (mer[pr].isnull())].shape[0]
    tn = nshuf - mer[shuffled & ~mer[pr].isnull()].shape[0]
    try:
        fdr = fp / (fp + tp)
    except ZeroDivisionError:
        fdr = np.nan
    try:
        tpr = tp / (tp + fn)
    except ZeroDivisionError:
        tpr = np.nan
    try:
        ppv = tp / (tp + fp)
    except ZeroDivisionError:
        ppv = np.nan
    try:
        f1s = 2 * ((tpr * ppv) / (tpr + ppv))
    except ZeroDivisionError:
        f1s = np.nan
    try:
        mcc = ((tp * tn) - (fp * fn)) / np.sqrt((tp + fp) * (tp + fn) *
                                                (fp + tn) * (tn + fn))
    except (ZeroDivisionError, FloatingPointError) as e:
        mcc = -np.inf
    try:
        pr_count = tp + fp + tn + fn
        assert (pr_count == all_calls) or ((pr_count + nshuf) == all_calls)
    except AssertionError:
        print('scoring', taxa)
        print(tp, fp, tn, fn, all_calls, nshuf)
        with open('mer.pckl', 'wb') as f:
            dill.dump((prediction, truth, mer), f)
        raise
    Scores = namedtuple('Scores', ('TP', 'FP', 'TN', 'FN', 'FDR', 'TPR', 'F1S',
                                   'MCC', 'nshuf', 'allcalls', 'time', 'Mem',
                                   'cpu', 'taxlevel'))
    return Scores(tp, fp, tn, fn, fdr, tpr, f1s, mcc, nshuf, all_calls,
                  timed.time, timed.mem, timed.cpu, taxa)


def run_lca(prefix, db, query, evalue, p_id, m_hit, p_hit, n_hit,
            asfile=False, **kwargs):
    """
    Function to run BASTA. It is a wrapper over the BASTA call to optimize its
    parameters
    :param prefix: Prefix of the output
    :param db: Database to use in the inference
    :param query: Query fasta file or fasta as a string 
    :param evalue: Evalue treshold
    :param p_id: Percent identity threshold
    :param m_hit: Same as basta's -m/--minimum option
    :param p_hit: Same as basta's -p/--maj_perc option
    :param n_hit: Same as basta's -n/--number option
    :param asfile: Wether to output the file or not
    :param kwargs: Other Key-word arguments
    """
    if not os.path.isfile('config.txt'):
        with open('config.txt', 'w') as c:
            c.write("query_id\t0\nsubject_id\t1\nalign_length\t6\nevalue\t3\n")
            c.write("pident\t2\n")
    names = ['qaccver', 'saccver', 'pident', 'evalue', 'qcovs', 'qlen',
             'length', 'staxid', 'stitle', 'Kingdom', 'Phylum', 'Class',
             'Order', 'Family', 'Genus', 'Species']
    directory = '%s_LCA_files' % prefix
    src_db = os.path.join(os.environ['HOME'], '.basta', 'taxonomy')
    if not os.path.isdir(directory):
        os.mkdir(directory)
    if asfile:
        outpref = '%s_lca_final' % prefix
    else:
        outpref = os.path.join(directory, 'LCA_ev%.0e_pid%d_m%d_p%d_n%d' % (
            evalue, p_id, m_hit, p_hit, n_hit))
    blast_file = '%s.blast' % outpref
    df = run_blast(db, query, evalue, p_id, 500, tophit=False)
    df.to_csv(blast_file, sep='\t', index=False, header=False)
    args = ['python2', path2basta_exe, 'sequence', blast_file,
            '%s.out' % outpref, basta_map, '-v', '%s.verbose' % outpref, '-e',
            str(evalue), '-i', str(p_id), '-m', str(m_hit), '-d', src_db, '-p',
            str(p_hit), '-n', str(int(n_hit)), '-b', 'T', '-x', 'T', '-c',
            'config.txt']
    st = run(args, stdout=PIPE, stderr=PIPE, env=myenv)
    try:
        st.check_returncode()
    except CalledProcessError:
        print(st.stdout, st.stderr)
        raise
    with open('LOG', 'a') as log:
        log.write('\n'.join([' '.join(args), st.stderr.decode('utf-8'),
                             st.stdout.decode('utf-8')]))
    basta = pd.read_csv('%s.out' % outpref, sep='\t', header=None,
                        names=['seq', 'Lineage', 'best_hit'])
    lin = basta.Lineage.str.strip(';').str.split(';', expand=True).rename(
        columns=taxlevels)
    lin = lin.replace('Unknown', None).reindex(columns=taxlevels.values(),
                                               fill_value='')
    basta = pd.concat([basta, lin], axis=1)
    # print('basta.columns', basta.columns)
    basta = basta.replace('', None)
    basta['Species'] = basta.Species.str.replace('_', ' ')
    return basta


def run_blast(db, query, evalue, p_id, max_target_seqs, tophit=True,
              asfile=False, **kwargs):
    q = query if asfile else '-'
    outfmt = '6 qaccver saccver pident evalue qcovs qlen length staxid stitle'
    args = ['blastn', '-db', db, '-query', q, '-evalue', str(evalue),
            '-perc_identity', str(int(p_id)), '-outfmt', outfmt,
            '-max_target_seqs', str(int(max_target_seqs)), '-num_threads',
            str(cpu_count())]
    st = execute(args, input=query.encode('utf-8'))
    df = pd.read_csv(BytesIO(st.stdout), sep='\t', header=None,
                     names=outfmt.split()[1:])
    if not tophit:
        return df
    taxa = df.stitle.str.strip(';').str.split(';', expand=True).rename(
        columns=taxlevels)
    df = pd.concat([df, taxa], axis=1).replace('', None)
    df.rename(columns={'qaccver': 'seq'}, inplace=True)
    top_hits = df.sort_values(by='evalue', ascending=True).groupby(
        'seq', as_index=False).first()
    return top_hits


def execute(args, input=None, env=dict()):
    st = run(args, input=input, stdout=PIPE, stderr=PIPE, env=env)
    try:
        st.check_returncode()
    except CalledProcessError:
        print(st.stderr.decode('utf-8'), st.stdout.decode('utf-8'))
        raise
    return st


def run_qiime(query, db, conf, **kwargs):
    env = os.environ
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
        try:
            tmpname = tmpdirname.name
        except AttributeError:
            tmpname = tmpdirname
        tempfile.tempdir = tmpname
        os.environ['TMPDIR'] = tmpname
        quer = os.path.join(tmpname, 'query')
        if query.startswith('>'):
            with open(quer, 'w') as q:
                q.write(query.replace('\n\n', '\n'))
        else:
            shutil.copy(query, quer)
        seq = os.path.join(tmpname, 'query_tmp.qza')
        cls = os.path.join(tmpname, 'qiime_res.qza')
        export = ['qiime', 'tools', 'import', '--type',
                  'FeatureData[Sequence]', '--input-path', quer,
                  '--output-path', seq]
        _ = execute(export, env=env)
        classify = ['qiime', 'feature-classifier', 'classify-sklearn',
                    '--i-classifier', db, '--i-reads', seq,
                    '--o-classification', cls, '--p-confidence',
                    str(conf)]  # ,
        # '--p-n-jobs', '28']
        try:
            _ = execute(classify, env=env)
        except:
            print(' '.join(classify))
            shutil.move(seq, 'query.qza')
            _ = execute(classify, env=env)
            raise
        get_taxonomy = ['qiime', 'tools', 'export', '--input-path', cls,
                        '--output-path', tmpname]
        _ = execute(get_taxonomy, env=env)
        df = pd.read_csv(os.path.join(tmpname, 'taxonomy.tsv'), sep='\t')
        try:
            lineage = df.Taxon.str.strip(';').str.split('; ', expand=True
                                                        ).apply(
                lambda x: x.str.split('__').str[1])
            try:
                lineage[6] = lineage[5] + ' ' + lineage[6]
            except KeyError:
                pass
            lineage = lineage.rename(columns=taxlevels).reindex(
                columns=taxlevels.values(), fill_value=None)
            df = pd.concat([df, lineage], axis=1)
        except:
            with open('df.pckl', 'wb') as pckl:
                dill.dump(df, pckl)
            raise
    return df.rename(columns={'Feature ID': 'seq', 'Taxon': 'taxa'})


def run_protax(query, db, conf, **kwargs):
    global_consensus = protax_global
    local_consensus = protax_local_isis if 'ISIS' in db else protax_local
    reference = protax_ref % db.split('_')[0]

    def tosingle(sequence):
        bl = sequence.split('\n')
        return '%s;\n%s' % (bl[0].strip(';'), ''.join(bl[1:]).upper())
    shuff = [tosingle(x) for x in query if 'shuffled' in x]
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
        try:
            tmpname = tmpdirname.name
        except AttributeError:
            tmpname = tmpdirname
        tempfile.tempdir = tmpname
        os.environ['TMPDIR'] = tmpname
        quer = os.path.join(tmpname, 'query')
        aln = os.path.join(tmpname, 'aln')
        with open(quer, 'w') as q:
            q.write(query.replace('\n\n', '\n'))
        seqs, subset = align(global_consensus, local_consensus, quer,
                             no_write=True)
        good = [tosingle(x) for x in seqs if 'shuffled' not in x] + shuff
        fas = '%s.fas' % aln
        with open(fas, 'w') as out:
            out.write('\n'.join(good))
            out.write('\n')
        args = [os.path.join(path2protax_scripts, 'classify_best2'),
                os.path.join(reference, '%staxonomy.priors' % db),
                os.path.join(base, '%s.aln' % db.strip('_')),
                os.path.join(reference, '%smodel.rseqs.numeric' % db),
                os.path.join(reference, '%smodel.pars' % db),
                os.path.join(reference, '%smodel.scs' % db), str(conf), fas]
        st = execute(args)
        result = []
        for line in st.stdout.strip().split(b'\n'):
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            bl = line.strip().split()
            name = bl.pop(0)
            b = pd.Series(bl)
            taxa = b[list(range(0, b.shape[0], 2))].reset_index(drop=True)
            prob = b[list(range(1, b.shape[0], 2))].astype(float).reset_index(
                drop=True)
            levels = taxa.str.split(',').apply(lambda x: len(x))
            min_levels = levels[levels == levels.max()].index
            p = prob[min_levels]
            d = {items[0].split(',')[-1]: items[1] for i, items in
                 enumerate(zip(taxa, prob))}
            idx = p[p == p.max()].index[0]
            best = taxa[idx].split(',')
            lineage = get_lineage(best[-1].replace('_', ' '))
            expand = dict(zip(taxlevels.values(), lineage.strip().split(';')))
            probs = ';'.join('%.2f' % d[x] for x in best)
            d = {'seq': name, 'taxa': lineage, 'probs': probs}
            d.update(expand)
            result.append(d)
    return pd.DataFrame(result)


def lineage_from_taxid(tid):
    if ncbi.get_lineage(tid) is None:
        return ';;;;;;'
    d = {v: k for k, v in ncbi.get_rank(ncbi.get_lineage(tid)).items()}
    lineage = [ncbi.translate_to_names([d[x]])[0] if x in d else '' for x in
               levels]
    return ';'.join(lineage)


def run_kraken(query, db, conf, **kwargs):
    cpu_str = str(cpu_count())
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
        try:
            tmpname = tmpdirname.name
        except AttributeError:
            tmpname = tmpdirname
        tempfile.tempdir = tmpname
        os.environ['TMPDIR'] = tmpname
        quer = os.path.join(tmpname, 'query')
        with open(quer, 'w') as q:
            q.write(query.replace('\n\n', '\n'))
        args = ['kraken2', '--db', db, quer, '--threads', cpu_str,
                '--confidence', str(conf)]
        # print(' '. join(args))
        st = execute(args, env=myenv)
    # tkcmd = "taxonkit lineage -i 3 | taxonkit reformat -i 6"
    # st = run(tkcmd, input=st.stdout, shell=True, stdout=PIPE)
    # st.check_returncode()
    df = pd.read_csv(BytesIO(st.stdout), sep='\t', header=None,
                     names=['C/U', 'seq', 'krakenid', 'len', 'kraken_lin',
                            'basta_lin'])
    df['reformat'] = df.krakenid.apply(lineage_from_taxid)
    taxa = df.reformat.str.strip(';').str.split(';', expand=True, n=7).reindex(
        columns=range(7)).rename(columns=taxlevels)
    df = pd.concat([df, taxa], axis=1).replace('', None)
    return df


def run_hmmufotu(query, db, seed, max_pd, err, method, prior, **kwargs):
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
        try:
            tmpname = tmpdirname.name
        except AttributeError:
            tmpname = tmpdirname
        tempfile.tempdir = tmpname
        os.environ['TMPDIR'] = tmpname
        quer = os.path.join(tmpname, 'query')
        if query.startswith('>'):
            with open(quer, 'w') as q:
                q.write(query.replace('\n\n', '\n'))
        else:
            shutil.copy(query, quer)
        args = ['hmmufotu', db, quer, '--single', '-C', '-v', '-p', '28', '-N',
                str(seed), '--fmt', 'fasta', '-d', str(max_pd), '-e', str(err),
                '-m', method, '--prior', prior, '-s', '1']
        # print(' '.join(args))
        st = execute(args, env=myenv)
    df = pd.read_csv(BytesIO(st.stdout), sep='\t', comment='#')
    u = df.taxon_anno[~df.taxon_anno.str.startswith('k__')].unique()
    comm = "taxonkit -j 28 name2taxid| taxonkit -j 28 lineage -i 2| " \
           "taxonkit -j 28 reformat -i3| cut -f 4| sed -e 's/^/k__/' -e " \
           "'s/;/; p__/1' -e 's/;/; c__/2' -e 's/;/; o__/3' " \
           "-e 's/;/; f__/4' -e 's/;/; g__/5' -e 's/;/; s__/6'"
    d = {i: check_output(comm, shell=True,
                         input=i.encode('utf-8')).strip().decode('utf-8') for i
         in u}
    df.replace(d, inplace=True)
    lineage = df.taxon_anno.str.strip(';').str.replace(' ', '').str.split(
        ';', expand=True).apply(lambda x: x.str.split('__').str[1])
    lineage[6] = lineage[6].str.split('_').str.join(' ')
    lineage = lineage.rename(columns=taxlevels).reindex(
        columns=taxlevels.values(), fill_value=None)
    cols = ['id', 'taxon_id', 'taxon_anno', 'anno_dist', 'Q_placement',
            'Q_taxon'] + lineage.columns.tolist()
    df = pd.concat([df, lineage], axis=1).reindex(columns=cols)
    return df.rename(columns={'id': 'seq', 'taxon_anno': 'taxa'})


def run_idtaxa(query, db, thresh, bootst, mind, **kwargs):
    prefix = kwargs['prefix']
    code_path =os.path.dirname(os.path.realpath(__file__))
    exe = os.path.join(code_path, 'idtaxa.R')
    args = 'Rscript %s -t %s -q %s -T %.2f -b %d -m %.2f'
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdirname:
        try:
            tmpname = tmpdirname.name
        except AttributeError:
            tmpname = tmpdirname
        tempfile.tempdir = tmpname
        os.environ['TMPDIR'] = tmpname
        quer = os.path.join(tmpname, 'query')
        if query.startswith('>'):
            with open(quer, 'w') as q:
                q.write(query.replace('\n\n', '\n'))
        else:
            shutil.copy(query, quer)
        command = args % (exe, db, quer, thresh, bootst, mind)
        print(command)
        _ = execute(command.split(), env=myenv)
        df = pd.read_csv(os.path.join(tmpname, 'idtaxa.tsv'), sep='\t')
        df = df.replace('unclassified_Root', 'None')
        lineage = df.taxa.str.strip(';').str.split(';', expand=True)
        lineage = lineage.rename(columns=taxlevels).reindex(
            columns=taxlevels.values(), fill_value=None)
    return pd.concat([df, lineage], axis=1)


def custom_cv(X, y, folds=3):
    kf = KFold(n_splits=folds)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X.loc[train_index, :], X.loc[test_index, :]
        y_train, y_test = y[train_index], y[test_index]
        yield X_train, X_test, y_train, y_test


def split_train(X, y, folds=3):
    def once(X, y, folds):
        X = X.sample(frac=1)
        arr = X.index.values
        spl_arr = np.array_split(arr, folds)
        test_idx = np.concatenate([a for a in spl_arr[:-1]])
        train_idx = spl_arr[-1]
        X_train, X_test = X.loc[train_idx, :], X.loc[test_idx, :]
        y_train, y_test = y[y.seq.isin(train_idx)], y[y.seq.isin(test_idx)]
        return X_train, X_test, y_train, y_test
    X_train, X_test, y_train, y_test = once(X, y, folds)
    while not (y_train.seq.str.contains('shuffled').any() and
               y_test.seq.str.contains('shuffled').any()):
        X_train, X_test, y_train, y_test = once(X, y, folds)
    assert X_train.index.isin(y_train.seq).all()
    assert X_test.index.isin(y_test.seq).all()
    assert (y_train.seq.str.contains('shuffled').any() and
            y_test.seq.str.contains('shuffled').any())
    return X_train, X_test, y_train, y_test


def compute_simple_grid(function, i, pnames, f_param, y_train, tax):
    """
    This function computes a grid search, that is, a search over a grid or set of parameters
    :param function: Function to be optimized. 
    :param i: Parameters values to optimize
    :param pnames: name of the parameter values to optimize
    :param f_param: Other parameters to the function (these would not be optimized)
    :param y_train: Training fold labels
    :param tax: Taxonomic level being evaluated
    """
    # Generate an empty namedtupled named timed
    timed = namedtuple('timed', ('time', 'mem', 'cpu'))  
    # filled with dummy variables. This is to make it compatible 
    # with the score function
    ti = timed(0, 0, 0)  
    # Generate a dictionary mapping parameter values with their name
    d = dict(zip(pnames, i))
    # Merge the parameters to be optimized with the static parameter
    z = {**d, **f_param}
    # Get a prediction with the function and the parameters
    prediction = function(**z)
    # Remove empty lines in the prediction
    prediction = prediction.replace(r'^\s*$', np.nan, regex=True)
    # Score the prediction using the traingset labels
    sc = score(prediction, y_train, ti, tax)
    # Transform the resulting namedtuple into a dataframe
    scd = pd.DataFrame(data=[sc])
    # Add the optimized parameters as a column in the dataframe
    scd = pd.concat([scd, pd.DataFrame(data=[d])], axis=1)
    # Add names to the index referring to the function name
    scd = scd.rename(index={0: function.__name__.split('_')[1]})
    return scd


def optimize_it(query, truthfile, db, function, taxlevel=6, n_jobs=-1,
                grid='simple', **kwargs):
    taxon = taxlevels[taxlevel]
    prefix = os.path.basename(query)
    prefix = os.path.splitext(prefix)[0]
    mock = prefix.split('_')[0]
    name = function.__name__.split('_')[1]
    pickle_pref = '%s_%s_%s' % (prefix, name, taxon)
    print('Optimizing %s in %s' % (name, taxon))
    if os.path.isfile('%s.pckl' % pickle_pref):
        print('Loading previous run')
        with open('%s.pckl' % pickle_pref, 'rb') as p:
            d = dill.load(p)
    else:
        X, y, truth = fasta2array(query, truthfile)
        params = product(*[kv for kv in kwargs.values()])
        pnames = [x[0] for x in product(kwargs)]
        if os.path.isfile('%s_split.pckl' % prefix):
            with open('%s_split.pckl' % prefix, 'rb') as f:
                X_train, X_test, y_train, y_test = dill.load(f)
        else:
            X_train, X_test, y_train, y_test = split_train(X.copy(deep=True),
                                                           y.copy(deep=True))
            with open('%s_split.pckl' % prefix, 'wb') as f:
                dill.dump((X_train, X_test, y_train, y_test), f)
        query_train = array2fasta(X_train)
        query_test = array2fasta(X_test)
        prefix = query[query.rfind('.fa')]
        f_param = dict(db=db, taxlevel=taxlevel, asfile=False,
                       prefix=prefix, query=query_train)
        if grid == 'simple':
            asfile = False
            if os.path.isfile('%s_training.tsv' % pickle_pref):
                scores = pd.read_csv('%s_training.tsv' % pickle_pref, sep='\t')
            else:
                tax = taxlevels[taxlevel]
                delayed_results = [dask.delayed(
                    compute_simple_grid(function, i, pnames, f_param, y_train,
                                        tax)) for i in tqdm(params)]
                with ProgressBar():
                    scores = dask.compute(*delayed_results)
                scores = pd.concat(scores)
                scores.to_csv('%s_training.tsv' % pickle_pref, sep='\t')
            highcols = ['MCC']
            if "hmmufotu" in name:
                highcols = ['F1S']
            if 'blast' in name:
                highcols += ['evalue', 'max_target_seqs']
                lowcols = ['p_id']
            elif 'lca' in name:
                highcols += ['evalue']
                lowcols = ['p_id', 'm_hit', 'p_hit']
            else:
                lowcols = None
            best = scores.nlargest(1, columns=highcols, keep='all')
            if lowcols is not None:
                best = best.nsmallest(1, columns=lowcols)
            best = best.reset_index(drop=True)
            try:
                best = best.loc[0, pnames].to_dict()
            except TypeError:
                print(best)
                raise
        else:
            print('  Performing Grid Search CV')
            asfile = True
            with ProgressBar():
                tuned_parameters = [kwargs]
                p = dict(program=function, db=db, taxlevel=taxlevel,
                         asfile=False, prefix=prefix)
                c = GridSearchCV(Estimator(**p), tuned_parameters, cv=3,
                                 n_jobs=n_jobs, scoring=mathews_scorer).fit(X,
                                                                            y)
            with open('%s_trainingres.pckl' % prefix, 'wb') as f:
                dill.dump(c.cv_results_, f)
            best = c.best_params_
        print('Processing best parameters in %s' % name)
        print(best)
        with ResourceProfiler(dt=0.01) as prof:
            start = time.time()
            results = function(prefix=prefix, db=db, query=query_test,
                               asfile=asfile, taxlevel=taxlevel, **best)
            elapsed = time.time() - start
            results = results.replace(r'^\s*$', np.nan, regex=True)
        with open('%s_cvresults.pckl' % pickle_pref, 'wb') as q:
            dill.dump((results, y), q)
        resources = pd.DataFrame(data=prof.results).mean()
        resources.rename(columns={'time': 'timestamp'})
        resources['time'] = elapsed
        sc = [score(results, y_test, resources, taxa) for taxa in
              taxlevels.values()]
        d = pd.DataFrame(data=sc)
        d['Method'] = name
        d['Mock'] = mock
        print(d)
        with open('%s.pckl' % pickle_pref, 'wb') as p:
            dill.dump(d, p)
    return d


def plot_scores(d, out='test.pdf'):
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    f, ax = plt.subplots()
    d.reindex(columns=['FDR', 'TPR', 'F1S', 'MCC']).plot.barh(
        rot=0, legend=False, ax=ax)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
              ncol=4, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(out)


def main(query, truthfile, taxlevel=6):
    query_name = os.path.splitext(os.path.basename(query))[0].split('_')[0]
    functions = {
        run_blast: (blastdb % query_name, dict(
            evalue=[1E-50, 1E-25, 1E-15, 1E-7, 1E-4, 1E-2],
            p_id=list(range(70, 105, 5)),
            max_target_seqs=[1, 25, 50, 75, 100, 500]), -1),
        run_lca: (blastdb % query_name, dict(
            evalue=[1E-50, 1E-25, 1E-15, 1E-7, 1E-4, 1E-2],
            p_id=list(range(70, 105, 5)), m_hit=list(range(1, 11)),
            p_hit=[60, 70, 80, 90, 99], n_hit=list(range(11))), -1),
        run_qiime: (qiime_class % query_name, {'conf': np.linspace(0.5, 1, 5)},
                    -1),
        run_kraken: (krakendb % query_name, {'conf': np.linspace(0, 1, 10)},
                     -1),
        run_hmmufotu: (hmmufotu_db % query_name, dict(
            seed=[1, 25, 50, 75, 100], max_pd=[0, 0.1, 1, 10, 'inf'],
            err=[1, 10, 20, 30, 40], method=['unweighted', 'weighted'],
            prior=['uniform', 'height']
        ), -1),
        run_idtaxa: (idtaxa_db % query_name, dict(
            thresh=[10, 20, 40, 60, 80, 100], bootst=[50, 100, 500, 1000],
            mind=np.arange(0.8, 1.05, 0.05)), -1),
        run_protax: (protaxdb % query_name, {'conf':[0.1, 0.05, 0.01, 0.001]}, -1)
    }
    d = pd.concat(
        [optimize_it(query, truthfile, tune[0], func, taxlevel=taxlevel,
                     n_jobs=tune[2], **tune[1]) for func, tune in
         functions.items()])
    print(d)
    return d


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('query', help='Fasta file with query sequences')
    parser.add_argument('truthfile', help='Tab delimited file with the maping '
                                          'of each of the sequences in query '
                                          'to species')
    parser.add_argument('-l', '--tax_level', default=6, type=int,
                        help='Taxonomic level to evaluate. 6 = species, '
                             '5 = genus ...')
    args = parser.parse_args()
    main(args.query, args.truthfile)
