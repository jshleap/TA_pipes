import os
import sys
import shelve
import pandas as pd
from tqdm import tqdm
from copy import deepcopy
from collections import defaultdict, Counter


def parse_fasta(filename):
    """
    Parse a fasta file and put it in a shelve DB and rename it with the prefix
    of the filename. It will also write a combined fasta with the renamed
    entries

    :param fil: name of fasta files
    :param fn2: name of the shelve
    :return: shelve db name
    """
    fn2 = '%s.shelve' % filename[:filename.rfind('.')]
    if not os.path.isfile(fn2):
        with shelve.open(fn2) as dic:
            name = None
            seq = ''
            with open(filename) as F:
                for line in tqdm(F, desc="Parsing %s" % filename):
                    if line.startswith('>'):
                        if name is not None:
                            dic[name] = seq
                        seq = ''
                        line = line.strip()
                        l = line[:-1] if line.endswith(';') else line
                        name = '%s' % l
                    else:
                        seq += line
                if name not in dic:
                    dic[name] = seq
    return fn2


class DummyNode(object):
    children = None

    def __init__(self, taxa, level, parent, data):
        self.__dict__.update(dict(taxa=taxa, level=level, nid=0, parent=parent,
                                  data=data, children=[], idx=0))
        self.nodes = self.get_nodes()

    def get_nodes(self):
        return [self]

    def preorder(self):
        yield self
        for child in self.children:
            yield from child.preorder()


class Node(object):
    def __init__(self, taxa, level, parent, data):

        self.level = level
        self.parent = parent
        self.taxa = taxa
        self.data = data
        self.children = []
        self.get_children()
        self.siblings = None
        self.sorted = defaultdict(list)
        self.sort_by_level()
        self.assign_id()

    @property
    def siblings(self):
        return self.__siblings

    @siblings.setter
    def siblings(self,_):
        try:
            self.__siblings = self.parent.children
        except AttributeError:
            self.__siblings = []

    @property
    def taxa(self):
        return self.__taxa

    @taxa.setter
    def taxa(self, taxa):
        self.__taxa = taxa
        try:
            if self.parent.taxa is not None:
                if self.parent.taxa != 'root':
                    self.__taxa = '%s,%s' % (self.parent.taxa, taxa)
        except AttributeError:
            self.__taxa = None

    def get_children(self):
        # counter = deepcopy(self.nid)
        level = deepcopy(self.level) - 1
        try:
            # if not self.parent.children and not self.parent.taxa  == 'unk':
            unk = 'unk'
            if self.parent.taxa and self.taxa != 'root':
                unk = '%s,unk' % self.taxa
            self.children.append(DummyNode(unk, self.level + 1, self,
                                           pd.DataFrame()))
            for name, data in self.data.groupby(level + 1):
                #print(name)
                # counter += 1
                n = Node(name, self.level + 1, self, data)
                self.children.append(n)
            self.children[0].siblings = self.children
        except (KeyError, AttributeError):
            self.children = []
            pass

    def sort_by_level(self):
        for x in self.preorder():
            self.sorted[x.level].append(x)

    def assign_id(self):
        # for i, n in tqdm(enumerate(self.preorder()), desc="Assigning ids in "
        #                                                   "n.taxa"):
        #     n.idx = i
        pbar = tqdm(desc='Assigning ids')
        counter = 0
        for k in self.sorted.keys():
            for n in self.sorted[k]:
                n.idx = counter
                counter += 1
                pbar.update()

    def preorder(self):
        yield self
        for child in self.children:
            yield from child.preorder()

    def __str__(self):
        return '\n'.join('%d\t%d\t%d\t%s' % (
            n.idx, n.parent.idx, n.level, n.taxa) for k in sorted(
            self.sorted.keys()) for n in self.sorted[k])


def seq2tax(dic, df):
    s2t = ''
    aln = ''
    for k, v in tqdm(dic.items(), desc="Formatting seq2tax and aln"):
        try:
            bl = k.split()
            name = bl[0]
            ser = df[df.index == name].iloc[0]
            n = ser.shape[0]
            tax = ','.join(ser).replace(' ', '_')
            seq = v.strip().replace('\n', '').upper()
            s2t += '%s\t%d\t%s\n' % (name[1:], n, tax)
            aln += '>%s %s\n%s\n' % (name, tax, seq) if '>' not in name else \
                '%s %s\n%s\n' % (name, tax, seq)
        except IndexError:
            continue
    return s2t, aln


def retain_only(lineage, mintaxalevel=4):
    """
    Retain only rows which taxonomic levels have more than two entries. Assumes
    a dataframe with only lineages in the column from 0 (kingdom) to species(6)
    :param lineage:
    :return:
    """
    for col in range(mintaxalevel, 7):
        c = Counter(lineage[col])
        p = [k for k, v in c.items() if v >= 2]
        lineage = lineage[lineage[col].isin(p)]
    return lineage


fasta = sys.argv[1]
mintaxalevel = int(sys.argv[2])
prefix = fasta[:fasta.rfind('.')]
fas = parse_fasta(fasta)
tax = '%s_protax.taxonomy' % prefix
s2i = '%s_protax.seqid2tax' % prefix
aln = '%s_protax.aln' % prefix
with shelve.open(fas) as dic, open(tax, 'w') as t, open(s2i, 'w') as s, open(
        aln, 'w') as a:
    keys = list(dic.keys())
    df = pd.DataFrame(keys)[0].str.split(' ', expand=True, n=1)
    lineage = df[1].str.split(';', expand=True)
    lineage.index = df[0]
    print('Original lineage shape:', lineage.shape)
    lineage = retain_only(lineage, mintaxalevel)
    lineage = lineage.loc[:, range(mintaxalevel, 7)]
    lineage = lineage.loc[:, lineage.nunique() > 1]
    print('Filtered lineage shape:', lineage.shape)
    lineage.index.name = 'ID'
    for col in lineage.columns:
        if len(lineage[col].unique()) == 1:
            lineage.drop(col, inplace=True, axis=1)
    lineage.rename(columns=dict(zip(lineage.columns, range(len(lineage.columns)
                                                           ))), inplace=True)
    lineage.replace(' ', '_', regex=True, inplace=True)
    s2, al = seq2tax(dic, lineage)
    s.write(s2.strip() + '\n')
    a.write(al.strip() + '\n')
    dummy = Node(None, 0, None, pd.DataFrame([]))
    root = Node('root', 0, dummy, lineage)
    t.write(str(root))
    print(root)