from subprocess import run, PIPE, CalledProcessError
from io import BytesIO
from tqdm import tqdm
import pandas as pd
import shelve
import sys
import os

myenv = {'PATH':os.environ['PATH'], 'PYTHONPATH': os.environ["PYTHONPATH"]}

def execute(args, input=None, env={}):
    st = run(args, input=input, stdout=PIPE, stderr=PIPE, env=env)
    try:
        st.check_returncode()
    except CalledProcessError:
        print(st.stderr.decode('utf-8'), st.stdout.decode('utf-8'))
        raise
    return st


def parse_fasta(filename):
    """
    Parse a fasta file and put it in a shelve DB and rename it with the prefix
    of the filename. It will also write a combined fasta with the renamed
    entries

    :param fil: name of fasta files
    :return: shelve db name
    """
    fn2 = '%s.shelve' % filename[:filename.rfind('.fasta')]
    with shelve.open(fn2) as dic:
        name = None
        seq = ''
        with open(filename) as F:
            for line in tqdm(F, desc="Parsing %s" % filename):
                if line.startswith('>'):
                    if name is not None:
                        dic[name] = seq
                    seq = ''
                    name = '%s' % (line.strip())
                else:
                    seq += line
            if name not in dic:
                dic[name] = seq
    return fn2


def taxonkit(uniques):
    args = ['taxonkit', '-j', '28', 'name2taxid']
    st = execute(args, input='\n'.join(uniques).encode('utf-8'), env=myenv)
    df = pd.read_csv(BytesIO(st.stdout), sep='\t', header=None, names=[
        'sp', 'taxid'])
    df = df.fillna(32644)
    return df


def build_db(dbname, fasta):
    dbname = dbname.replace('.', '_')
    args0 = ['kraken2-build', '--download-taxonomy', '--db', dbname]
    args1 = ['kraken2-build', '--download-library', 'protozoa',
            '--db', dbname, '--threads',  '28']
    args2 = ['kraken2-build', '--add-to-library', fasta, '--db',
             dbname, '--threads',  '28']
    args3 = ['kraken2-build', '--build', '--db', dbname]
    _ = execute(args0, input=None, env=myenv)
    _ = execute(args1, input=None, env=myenv)
    _ = execute(args2, input=None, env=myenv)
    _ = execute(args3, input=None, env=myenv)


fas = sys.argv[1]
fasta = parse_fasta(fas)
pref = fas[:fas.rfind('.')]
with shelve.open(fasta) as dic, open('%s.kraken.fas' % pref, 'w') as o:
    keys = set(' '.join(x.strip().strip(';').split()[1:]).split(';')[-1]
            for x in dic.keys())
    taxa = taxonkit(keys)
    for k, v in dic.items():
        bl = k.strip().strip(';').split()
        acc = bl[0]
        lineage = ' '.join(bl[1:]).split(';')
        sp = lineage[-1]
        taxid = int(taxa[taxa.sp == sp].taxid.iloc[0])
        o.write('%s|kraken:taxid|%s %s\n%s' % (acc, str(taxid), ';'.join(
            lineage), v))

build_db(pref, '%s.kraken.fas' % pref)

