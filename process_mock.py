"""
**Copyright (C) 2019  Jose Sergio Hleap**

This script takes a blast table and assigns as much as it can with greater than
99% id, and make a tree with the rest
"""
import argparse
import os
import shelve
import warnings
from collections import Counter
from io import BytesIO, StringIO
from itertools import combinations
from subprocess import run, PIPE, CalledProcessError
from time import sleep

import PyQt5
import dill
import numpy as np
import pandas as pd
import skbio as sk
from Bio import Entrez, AlignIO
from ete3 import PhyloTree, NCBITaxa
from joblib import Parallel, delayed
from scipy import stats
from tqdm import tqdm

ncbi = NCBITaxa()


def parse_fasta(filename, fn2=None, dictionary=False):
    """
    Parse a single fasta file and put it in a shelve DB

    :param filename: fasta file
    :param fn2: name of the shelve
    :param rename: boolean of whether to append prefix to each sequence
    :return: shelve db name
    """
    if isinstance(filename, shelve.DbfilenameShelf):
        Warning('THIS IS ALREADY A SHELF!!!')
        return filename
    elif isinstance(filename, bytes):
        fnc = BytesIO
        fn = 'shelve' if fn2 is None else fn2
    else:
        fnc = open
        fn = filename[:filename.find('.fa')]
    if dictionary:
        db = dic = {}
    else:
        db = '%s.shelve' % fn
        dic = shelve.open(db)
    name = None
    seq = ''
    with fnc(filename) as F:
        for line in F:
            line = line.decode('utf-8') if isinstance(line, bytes) else line
            if line.startswith('>'):
                if name is not None:
                    dic[name] = seq
                seq = ''
                name = line.strip()
            else:
                seq += line
        if name not in dic:
            dic[name] = seq
    if isinstance(db, shelve.DbfilenameShelf):
        db.close()
    return db


def stdin_run(args, inpt, return_error=False, **kwargs):
    file_name = None
    if isinstance(inpt, tuple):
        file_name, inpt = inpt
    inpt = inpt.encode('utf-8') if isinstance(inpt, str) else inpt
    exe = run(args, input=inpt, stderr=PIPE, stdout=PIPE, **kwargs)
    try:
        exe.check_returncode()
    except CalledProcessError:
        raise Exception(exe.stdout, exe.stderr)
    if return_error:
        return exe.stdout, exe.stderr, file_name
    else:
        return exe.stdout, file_name


def iter_fasta(fn, done=None):
    if not done:
        done = []
    with shelve.open(fn) as fastas:
        for header, sequence in fastas.items():
            if header not in done:
                done.append(header)
                yield '%s\n%s' % (header, sequence)


def single_blast(args, inpt):
    e = True
    while e:
        o, e, f = stdin_run(args, inpt, return_error=True)
        if e:
            print(e)
        if not o:
            o, e, f = stdin_run(args, inpt, return_error=True)
            if e:
                print(e)
    return o, f


def blast(db, query, evalue, max_target_seqs=50, cpus=-1, out='hits.hits'):
    outfmt_str = 'qaccver saccver pident evalue qcovs length staxid stitle ' \
                 'sstart send'
    fasta = parse_fasta(query)
    if os.path.isfile(out):
        blasts = pd.read_csv(out, sep='\t', header=None,
                             names=outfmt_str.split())
    else:
        args = ['blastn', '-db', db, '-query', '-', '-evalue', str(evalue),
                '-outfmt', '6 %s' % outfmt_str, '-max_target_seqs',
                str(max_target_seqs), '-perc_identity', '90']
        blasts = Parallel(n_jobs=cpus)(delayed(single_blast)(args, inp)
                                       for inp in iter_fasta(fasta))
        blasts, _ = zip(*blasts)
        blasts = [pd.read_csv(BytesIO(x), sep='\t', header=None,
                              names=outfmt_str.split()) for x in blasts]
        blasts = pd.concat(blasts)
        blasts.to_csv(out, sep='\t', index=False, header=False)
    return blasts, fasta


def process_othernames(obj):
    a_list = []
    for k, v in obj.items():
        if k != 'Name':
            a_list += v if isinstance(v, list) else [v]
    return a_list


def get_sp_synonyms(email, tax_id_series):
    mapping = {}
    tax2sp = {}
    Entrez.email = email
    taxids = ','.join(str(x) for x in tax_id_series.unique())
    with Entrez.efetch(db="Taxonomy", id=taxids, retmode="xml") as handle:
        rec = Entrez.read(handle)
    for r in rec:
        sp = r['ScientificName']
        tax2sp[int(r['TaxId'])] = sp
        try:
            mapping[sp] = process_othernames(r['OtherNames'])
        except KeyError:
            mapping[sp] = []
    return mapping, tax2sp


def per_group(name, d):
    nd = pd.DataFrame()
    sps = d.sp.unique()
    nd.loc[0, 'Name'] = name
    nd.loc[0, 'Accessions'] = ';'.join(x for x in d.saccver.unique())
    return pd.concat((nd, pd.DataFrame(sps).T), axis=1)


def get_sequences(ids, db):
    args = ['blastdbcmd', '-db', db, '-entry_batch', '-']
    try:
        fas, _ = stdin_run(args, '\n'.join(ids))
    except:
        ids = [i.split('|')[0] for i in ids]
        fas, _ = stdin_run(args, '\n'.join(ids))
    fas = parse_fasta(fas, dictionary=True)
    fas = {k.split()[0]: v for k, v in fas.items()}
    return fas


def trimaln(aln, target_ids, gaps=0.9):
    # Read the alignemnt into biopython structure
    aln = AlignIO.read(StringIO(aln), 'fasta')
    dfaln = pd.DataFrame(aln)
    dfaln.index = [x.id for x in aln]
    # get only the targets
    aln = AlignIO.MultipleSeqAlignment([x for x in aln if x.id in target_ids])
    nseqs = len(aln)
    c = pd.DataFrame(aln).apply(lambda x: sum(x == '-') / nseqs, axis=0)
    dfaln = dfaln.loc[:, c[c < gaps].index]
    e = dfaln.apply(lambda x: sum(x == '-') / dfaln.shape[1], axis=1)
    dfaln = dfaln[e < gaps].drop_duplicates()
    l = ['>%s\n%s\n' % (x[0], ''.join(x[1]).upper().strip()) for x in
         dfaln.iterrows()]
    return '\n'.join(l).replace('\n\n', '\n')


def do_alnntree(pref, ndf, fasta, refs, congen, targetids, gaps=0.9, cpus=-1):
    # TODO: add checkpoint to avoid repeating
    to_phy = congen
    for name, data in ndf.groupby('saccver'):
        # mi = data.sstart.min()
        # ma = data.send.max()
        tx = data.staxid.iloc[0]
        try:
            seq = refs['>%s' % name].replace('\n', '').strip()  # [mi-1:ma+1]
        except KeyError:
            name = name.split('|')[0]
            seq = refs['>%s' % name].replace('\n', '').strip()
        to_phy += '>%d.%s\n%s\n' % (tx, name, seq)
    with shelve.open(fasta) as dic:
        for h, s in dic.items():
            if h.strip()[1:] in targetids:
                print(h)
                to_phy += '%s\n%s\n' % (h, s.strip().replace('\n', ''))
            else:
                print(h, 'not in')
    aln, _ = stdin_run(['mafft', '--thread', str(cpus), '--auto', '-'], to_phy)
    trm = trimaln(aln.decode('utf-8'), targetids, gaps=gaps)
    tre, _ = stdin_run(['fasttreeMP', '-nt', '-gtr', '-gamma'], trm)
    tre = tre.strip()[:-1].replace(b';', b'-').decode('utf-8') + ';'
    t = PhyloTree(tre, sp_naming_function=lambda name: name.split('.')[0])
    with open('%s.aln' % pref, 'w') as al, open('%s.treepickle' % pref, 'wb') \
            as tp:
        al.write(trm)
        t.write(outfile='%s.tree' % pref)
        dill.dump(t, tp)
    tax2 = t.annotate_ncbi_taxa()
    fix_species(t)
    print(t)
    return t, tax2


def subtree2matrix(node):
    m = pd.DataFrame()
    for a in combinations(node.get_leaf_names(), 2):
        dist = node.get_distance(a[0], a[1])
        m.loc[a] = dist
        m.loc[a[::-1]] = dist
    indices = m.sort_index().index
    m = m.fillna(0).loc[indices, indices]
    dm = sk.DistanceMatrix(m)
    groups = [x.split('.')[0] if '.' in x else 'sq' for x in m.index]
    return dm, groups


def subtree2matrix2(node):
    m = pd.DataFrame()
    for a in combinations(node.get_leaf_names(), 2):
        dist = node.get_distance(a[0], a[1])
        m.loc[a] = dist
        m.loc[a[::-1]] = dist
    indices = m.sort_index().index
    m = m.fillna(0).loc[indices, indices]
    return m


def fix_species(node):
    for n in node.get_descendants():
        if not n.is_leaf() and n.support <= 0.7:
            n.delete()
        elif n.is_leaf():
            taxid = n.taxid
            if taxid is not None and n.rank != 'species':
                taxmap = {v: k for k, v in ncbi.get_rank(n.lineage).items()}
                n.sci_name = ncbi.get_taxid_translator([taxmap['species']])[
                    taxmap['species']]
                n.taxid = taxmap['species']
                n.rank = 'species'
                n.species = taxid


def get_species_clade(tree, taxid):
    nodes = tree.search_nodes(taxid=taxid)
    if len(nodes) == 1:
        node = nodes[0]
        if node.up.taxid is not None:
            return None
        else:
            node = node.up
    node = tree.get_common_ancestor(nodes)
    nsp = set(x.taxid for x in node if x.taxid is not None)
    while (len(nsp) == 1) and not node.is_root():
        node = node.up
        nsp = set(x.taxid for x in node if x.taxid is not None)
    if node.support <= 0.7:
        return node, set(x.taxid for x in node if x.taxid is not None)
    else:
        node = node.get_common_ancestor(node.search_nodes(taxid=taxid))
        return node, set(x.taxid for x in node if x.taxid is not None)


def get_target_clade(leafid, tree, done=None, sig_level=0.05):
    if not done:
        done = {()}
    if leafid not in done:
        node = (tree & leafid)
        while True:
            node = node.up
            node_support = node.support
            sps = set(n.taxid for n in node if n.taxid)
            if not [x.name for x in node if x.taxid is not None]:
                continue
            elif node_support < 0.7:
                continue
            elif node.up.support < 0.7:
                c = node.up.up
                if (c.support < 0.7) and not (test_outlier(c, sig_level)):
                    continue
                else:
                    sps = set(x for x in node.up.get_species() if x.isnumeric())
                    return node.up, sps
            elif len([x for x in node.up.get_species() if x.isnumeric()]) != 1:
                return node, sps
            else:
                return node, set(n.taxid for n in node if n.taxid)


def test_outlier(node, sig_level):
    sig_level = stats.norm.ppf(1 - (sig_level / 2))
    m = subtree2matrix2(node)
    return (m.apply(stats.zscore).apply(np.abs) > sig_level).values.any()


def test_id(tree, tax2name, intended, sig_level=0.05, prefinquery='sq'):
    done = set()
    intended = pd.read_csv(intended, header=None, names=['sps'])
    assignment = []
    attributes = ['support', 'taxid', 'sci_name', 'name']
    ref_sps = set(x.taxid for x in tree.iter_leaves() if x.taxid is not None)
    with tqdm(total=len(ref_sps), desc="Processing Targets") as tq:
        for r in ref_sps:
            tq.desc = 'Processing %s' % r
            tq.update()
            tup = get_species_clade(tree, r)
            if tup is None:
                continue
            else:
                c, sps = tup
                queries = [n.name for n in c if n.taxid is None]
                clade = ';'.join([n.name.split('.')[1] for n in c
                                  if n.taxid is not None])
                if len(sps) == 1:
                    sp = int(list(sps)[0])
                else:
                    sp = sorted((c.get_distance(x), int(x.taxid)) for x in c if
                                x.taxid is not None)[0][1]
                if test_outlier(c, sig_level):
                    line = '\n\nAssignment to species %d created outliers with ' \
                           'tree:\n%s\n\n'
                    a_t = c.get_ascii(attributes=attributes)
                    a = [x.sci_name for x in c.iter_leaves() if x.taxid is not
                         None]
                    if (len(set(a)) == 1) and intended.sps.isin(a).any():
                        for query in queries:
                            l = '%s\t%s\t%s' % (query, clade, tax2name[sp])
                            assignment.append(l)
                            done.add(query)
                    else:
                        with open('Outliers_to_check.txt', 'a') as out:
                            out.write(line % (sp, a_t))
                        continue
                if len(queries) == 0:
                    continue
                for query in queries:
                    l = '%s\t%s\t%s' % (query, clade, tax2name[sp])
                    assignment.append(l)
                    done.add(query)
    missing = set([x.name for x in tree if prefinquery in x.name]
                  ).difference(done)
    return assignment, missing


def get_congenerics(pref, ndf, geneofinterest='COI', maxn=3):
    staxids = ndf.staxid.unique()
    common = dict(db="nucleotide", idtype="acc")
    results = ''
    pcklf = '%s_congenerics.pckl' % pref
    print("Downloading Congenerics")
    lins = ncbi.get_lineage_translator(staxids)
    genera = set(k for l in lins.values() for k, v in ncbi.get_rank(l).items()
                 if v == 'genus')
    idlist = set()
    if os.path.isfile(pcklf):
        warnings.warn("Using previous pickle")
        with open(pcklf, 'rb') as pckl:
            fasta = dill.load(pckl)
    else:
        with tqdm(total=len(genera)) as tq:
            for genus in genera:
                tq.desc = "Getting %s" % ncbi.get_taxid_translator([genus])
                tq.update()
                query = 'txid%d[Organism] AND %s[Gene]' % (genus,
                                                           geneofinterest)
                with Entrez.esearch(term=query, **common) as handle:
                    search_results = Entrez.read(handle)
                    # webenv = search_results["WebEnv"]
                    # idlist += search_results["QueryKey"]
                    idlist.update(search_results["IdList"])
                    sleep(1)
            with Entrez.efetch(rettype="fasta", retmode="text",
                               id=list(idlist), **common) as efetch:
                results += efetch.read()
            fasta = results.strip().split('\n\n')
            with open(pcklf, 'wb') as pckl:
                dill.dump(fasta, pckl)
    res = ''
    sps_counts = Counter()
    for fas in fasta:
        spl = fas.split(' ')
        acc = spl[0][1:]
        sp = ' '.join(spl[1:3])
        try:
            taxid = ncbi.get_name_translator([sp])[sp][0]
        except KeyError:
            continue
        seq = ''.join(spl[-1].split('\n')[1:])
        sps_counts.update([sp])
        if not (spl[0][1:] in ndf.saccver.tolist()) and (
                sps_counts[sp] < maxn):
            res += '>%d.%s\n%s\n' % (taxid, acc, seq)
    return res


def main(email, query='all_samples.fasta', intended='intended.txt', db='nt',
         dbname='nt', evalue=1e-3, gaps=0.9, cpus=-1, ntop=3, sig_level=0.05,
         geneofinterest='COI', prefinquery='sq'):
    pref = query[:query.find('.fasta')]
    dbname = os.path.split(db)[-1] if dbname is None else dbname
    outblast = '%s_%s.hits' % (pref, dbname)
    df, fasta = blast(db, query, evalue, out=outblast)
    if 'nt' not in db:
        name2tax = ncbi.get_name_translator(
            df.stitle.str.split(';', expand=True).iloc[:, 6].unique())
        df.staxid = df.stitle.str.split(';', expand=True).iloc[:, 6].apply(
            lambda x: name2tax[x][0])
    hundreds = df[(df.pident >= 99.5) & (df.qcovs == 100)].copy(deep=True)
    Entrez.email = email
    if not hundreds.empty:
        synonyms, tax2sp = get_sp_synonyms(email, hundreds.staxid)
        hundreds['sp'] = hundreds.staxid.map(tax2sp)
        hundreds['other_names'] = hundreds.sp.map(synonyms)
        gr100 = hundreds.groupby('qaccver')
        d = pd.concat((per_group(name, data) for name, data in gr100),
                      sort=False).reindex()
        d.to_csv('%s_hundreds.tsv' % pref, sep='\t', index=False)
        not100 = set(df.qaccver.unique()).difference(d.Name)
    else:
        not100 = set(df.qaccver)
    ndf = df[df.qaccver.isin(not100) & (df.qcovs > 70)]
    ndf = ndf.sort_values(by=['evalue', 'pident', 'qcovs'], ascending=[
        True, False, False]).groupby('qaccver')
    ndf = pd.concat([w.groupby('staxid').head(3) for _, w in ndf]).sort_values(
        by='qcovs', ascending=False).groupby('qaccver').head(ntop)
    congenerics = get_congenerics(pref, ndf, geneofinterest=geneofinterest,
                                  maxn=ntop)
    refs = get_sequences(df[(df.qcovs > 80) & (df.pident > 90) &
                            (df.evalue < 1E-20)].saccver.unique(), db)
    tree, tax2 = do_alnntree(pref, ndf, fasta, refs, congenerics, not100,
                             gaps=gaps, cpus=cpus)
    tax2names, tax2lineages, tax2rank = tax2
    assignment, missing = test_id(tree, tax2names, intended, sig_level,
                                  prefinquery)
    with open('%s_assigned_by_phylogeny.txt' % pref, 'a') as a, open(
            '%s.treepickle' % pref, 'wb') as tp:
        dill.dump((tree, tax2, missing), tp)
        a.write('\n'.join(assignment))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('email', help='email for Entrez queries')
    parser.add_argument('query', help='Fasta file with query sequences')
    parser.add_argument('intended', help='file with intended list of species')
    parser.add_argument('-d', '--db', help='Path to reference database to use',
                        default='nt')
    parser.add_argument('-b', '--dbname', help='name of reference database',
                        default=None)
    parser.add_argument('-e', '--evalue', help='Evalue threshold for Blast',
                        type=float, default=0.003)
    parser.add_argument('-g', '--gaps', help='Requested fraction of ungapped '
                                             'columns in alignment',
                        type=float, default=0.9)
    parser.add_argument('-c', '--cpus', help='Number of cpus to use',
                        type=int, default=-1)
    parser.add_argument('-n', '--ntop', help='Number of top hits to retain',
                        type=int, default=3)
    parser.add_argument('-m', '--marker', help='NCBI-compliant gene name ',
                        default='COI')
    parser.add_argument('-s', '--sig_level', help='Significance level on '
                                                  'outlier test',
                        type=float, default=0.05)
    parser.add_argument('-p', '--prefinquery', help=(
        'prefix of all sequences of interest in query file'), default='sq')

    args = parser.parse_args()
    main(args.email, args.query, args.intended, db=args.db, evalue=args.evalue,
         gaps=args.gaps, cpus=args.cpus, ntop=args.ntop, dbname=args.dbname,
         sig_level=args.sig_level, geneofinterest=args.marker,
         prefinquery='sq')
