"""
**tax4dada.py
** Copyright (C) 2018  Jose Sergio Hleap

With a fasta file from either BOLD or NCBI, translate the fasta into
dada compatible training files. That is, include the apropriate lineage
as names. This will generate two files: the lineage file and the species
file. They should be essentially the same, unless incomplete taxonomic
information is available.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jose.hleaplozano@mcgill.ca, dayana.salas@dal.ca

Requirements:
-------------

Python modules:
1. ete3 module
2. requests
3. pandas
4. joblib
5. eutils
6. biopython

System Programs:
1. NCBI eutils on path

Citation
========
If you use this software, please cite:
"""
from ete3 import NCBITaxa
import eutils
import pandas as pd
from multiprocessing import Process, JoinableQueue
from joblib import delayed, Parallel
from Bio import SeqIO
import time
import optparse
from hashlib import md5
import sqlite3
from sqlalchemy import create_engine
import os
import tqdm
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve

# Globals #############################################################
VALID_DBS = {"NCBI", "BOLD"}
FULL_TAX = ['superkingdom', 'kingdom', 'phylum', 'subphylum',
            'superclass', 'class', 'superorder', 'order', 'suborder',
            'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily',
            'genus', 'species']
SIX = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
#######################################################################

# Functions ###########################################################

# make only functions to parallelize??


def for_sqlalchemy(fns=None, path=os.getcwd(), sep='\t', compression='gzip',
                   create=True):
    """
    Create an sql database from sep delimited file.
    :param path: Path to database
    :param create: whether to create or load
    :param compression: Compression type
    :param sep: Delimiter in file
    :param fns: list with Filenames
    :return: sql connexion
    """
    db = 'sqlite:////%s' % os.path.join(path,'accession2taxid.db')
    print("Creating", db, "database")
    engine = create_engine(db, pool_pre_ping=True)  # , echo=True)
    db_params = dict(con=engine, index=False, if_exists='append')
    read_params = dict(sep=sep, compression=compression, chunksize=500000)
    if create or not engine.dialect.has_table(engine, 'taxas'):
        print("Creating accession2taxid.db...")
        open('accession2taxid.db', 'w').close()
        for fn in fns:
            print('Processing', fn, '...')
            reader = pd.read_table(fn, **read_params)
            for row in tqdm(reader):
                row.to_sql('taxas', **db_params)
    return engine


def create_sqlite3(fn, sep='\t', compression='gzip'):
    """
    Create an sql database from sep delimited file. Inpiration and refctoring
    from P2's csv-to-sqlite.py (https://gist.github.com/p2/7797584)
    :param compression: Compression type
    :param sep: Delimiter in file
    :param fn: Filename
    :return: sql connexion
    """
    conn = sqlite3.connect('accession2taxid.db')
    conn.isolation_level = 'EXCLUSIVE'
    cursor = conn.cursor()
    create_sql = 'CREATE TABLE rows '
    insert_sql = 'INSERT INTO rows '
    reader = pd.read_table(fn, sep=sep, compression=compression, chunksize=1)
    for i, row in enumerate(reader):
        if i == 0:
            # fill names
            varchars = ['"%s" VARCHAR' % field for field in row.columns]
            q = ' '.join('?' * row.shape[1])
            create_sql += '(%s)' % ', '.join(varchars)
            cursor.execute(create_sql, ())
            insert_sql += '("%s") VALUES (%s)' % ('", "'.join(row.columns), q)
        params = tuple(row.itertuples(index=False, name=None))[0]
        print(insert_sql)
        cursor.execute(insert_sql, params)
    conn.commit()
    conn.isolation_level = None
    return conn


# TODO: make it update when needed instead of redoing everything. For now, this
# has to do

def update_dumps(path=os.getcwd(), disable=False):
    """
    Updates the dumps for the accession2taxid files
    :param path: Path whefre the dumps and the database should be
    :returns: engine connection
    """
    if not disable:
        def _inner(fn):
            updated = False
            if check_if_up_to_date(fn):
                print('Updating %s' % fn)
                updated = True
                urlretrieve(url % fn, os.path.join(path, fn))
                print('Done. Parsing...')
            return updated

        url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/%s"
        files = ['nucl_est.accession2taxid.gz', 'nucl_gss.accession2taxid.gz',
                 'nucl_wgs.accession2taxid.gz', 'nucl_gb.accession2taxid.gz',
                 'prot.accession2taxid.gz']
        updates = Parallel(n_jobs=-1)(delayed(_inner)(fn) for fn in files)
    engine = for_sqlalchemy(fns=files, create=any(updates))
    return engine


def check_if_up_to_date(fn):
    """
    Check if the local database is up to date or should be updated. Modified
    from ete toolkit ncbiquery
    :param fn: dump filename
    :return: boolean
    """
    ftp = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/%s.md5"
    (md5_filename, _) = urlretrieve(ftp % fn)
    with open(md5_filename, "r") as md5_file:
        md5_check = md5_file.readline().split()[0]
    do_download = False
    try:
        local_md5 = md5(open(fn, "rb").read()).hexdigest()
        if local_md5 != md5_check:
            do_download = True
            print('Downloading %s from NCBI FTP site (via HTTP)...' % fn)
        else:
            print('Local %s seems up-to-date' % fn)
    except FileNotFoundError:
        do_download = True
        print('Downloading %s from NCBI FTP site (via HTTP)...' % fn)
    return do_download


def initNCBI():
    """
    Build the dabase (if not build before), and update its contents
    :return: connection
    """
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    return ncbi


def lineage2dict(tax_id):
    """
    Translate a taxid lineage to names avoiding no ranks
    :param tax_id: NCBItaxID to get the lineage from
    :return: dictionary mapping the names and the ranks
    """
    print('Retrieving lineages for all sequences...')
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(tax_id)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    ranks_list, names_list = [], []
    for i, k in enumerate(ranks.keys()):
        if ranks[k] == 'no rank':
            continue
        else:
            ranks_list.append(ranks[k])
            names_list.append(names[k])
    return dict(zip(ranks_list, names_list))


def get_NCBItaxid(accession, email=None, api_key=None, engine=None):
    """
    Access NCBI to retrieve taxid information from the accession number
    :param engine: connection to local database
    :param api_key: API key for 10 requests per second
    :param accession: NCBI accession number
    :return: taxid
    esearch -db nuccore -query LC315594.1|elink -target taxonomy |efetch
    -format uid
    """
    print('Get taxid with engine', engine)
    if engine is not None:
        query = "SELECT taxid FROM taxas WHERE accession='%s'" % accession
        tax_id = engine.execute(query).first()[0]
    else:
        if not email:
            raise BaseException('Please provide a valid email address')
        try:
            query_service = eutils.QueryService(email=email, api_key=api_key)
            xml = query_service.elink({'dbfrom': 'nuccore', 'db': 'taxonomy',
                                       'id': accession}).decode('utf')
            tax_id = xml[xml.rfind('<Id>') + 4: xml.rfind('</Id>')]
        except eutils.exceptions.EutilsRequestError:
            print(accession)
            raise
    return tax_id


def process_ncbi_per_entry(seq, lineage_fn, species_fn, tax_levels=SIX,
                           email=None, api_key=None, engine=None):
    """
    Process one sequence at a time
    :param lineage_fn: Oufilename of the lineage file
    :param species_fn: Oufilename of the species file
    :param tax_levels: NCBI taxonomic levels to retain
    :param api_key: Key to NCBI API
    :param seq: SeqRecord from a SeqIO read
    :return: Tuple with the two ouputs in string format
    """
    get_params = dict(email=email, api_key=api_key, engine=engine)
    print('Processing %s' % seq.id)
    seq_str = str(seq.seq)
    accession = seq.id[:seq.id.find('.')]
    lineage_dict = lineage2dict(get_NCBItaxid(accession, **get_params))
    lineage_str = ';'.join([lineage_dict[k] for k in tax_levels if k in
                            lineage_dict])
    if 'species' in lineage_dict:
        species_str = ' '.join([seq.id, lineage_dict['species']])
    else:
        species_str = None
    # lock = RLock()
    # with lock, open(lineage_fn, 'a') as lfn, open(species_fn, 'a') as sfn:
    with open(lineage_fn, 'a') as lfn, open(species_fn, 'a') as sfn:
        lfn.write('>%s\n%s\n' % (lineage_str, seq_str))
        if species_str:
            sfn.write('>%s\n%s\n' % (species_str, seq_str))
    if engine is None:
        time.sleep(1)


def process_bold_per_entry(seq, lineage_fn, species_fn, tax_levels=SIX):
    raise NotImplementedError


FUNCS = {"NCBI": process_ncbi_per_entry, "BOLD": process_bold_per_entry}


def main(fasta_fn, lineage_fn, species_fn, tax_levels=SIX, cpus=-1, db='NCBI',
         api_key=None, update_taxdb=False, email=None, dumps_path=os.getcwd(),
         disable=False):
    if db not in VALID_DBS:
        raise ValueError("results: db must be one of %r." % VALID_DBS)
    engine = update_dumps(dumps_path, disable=disable)
    seqs = SeqIO.parse(fasta_fn, 'fasta')
    if db == 'NCBI':
        n_requests = 8 if api_key else 2
        n_batch = 1
        if update_taxdb:
            _ = initNCBI()
    else:
        n_requests = 'all'
        n_batch = 'auto'
    p_params = dict(n_jobs=cpus, pre_dispatch=n_requests, batch_size=n_batch)
    f_params = dict(tax_levels=tax_levels, email=email, api_key=api_key,
                    engine=engine)
    # for seq in seqs:
    #     FUNCS[db](seq, lineage_fn, species_fn, **f_params)
    Parallel(**p_params)(delayed(FUNCS[db])(seq, lineage_fn, species_fn,
                                            **f_params) for seq in seqs)


#######################################################################


if __name__ == '__main__':
    opts = optparse.OptionParser(usage='%prog [options] ')
    opts.add_option('--fasta_fn', '-f', action='store',
                    help='Fasta file to process')
    opts.add_option('--lineage_fn', '-l', action='store',
                    help='Filename of the lineage outputfile')
    opts.add_option('--species_fn', '-s', action='store',
                    help='Filename of the species outputfile')
    opts.add_option('--tax_levels', '-t', action='store', default=SIX,
                    help='List of the levels to retain')
    opts.add_option('--cpus', '-c', action='store', default=-1, type=int,
                    help='Max number of cpus to use')
    opts.add_option('--db', '-d', action='store', default='NCBI', type=str,
                    help='Format of the fasta file (NCBI or BOLD)')
    opts.add_option('--api_key', '-a', action='store', default=None,
                    help='NCBI API key for more requests.')
    opts.add_option('--update_taxdb', '-u', action='store_true', default=False,
                    help='Update the database before running')
    opts.add_option('--email', '-e', action='store', default=None,
                    help='Email for NCBI requests')
    opts.add_option('--dumps_path', '-D', action='store', default=os.getcwd(),
                    help='Path where accession2taxid dumps are')
    opts.add_option('--disable_check', action='store', default=False,
                    help='Disable database check to update')
    opt, arg = opts.parse_args()

    main(opt.fasta_fn, opt.lineage_fn, opt.species_fn,
         tax_levels=opt.tax_levels, cpus=opt.cpus, db=opt.db,
         api_key=opt.api_key, update_taxdb=opt.update_taxdb, email=opt.email,
         dumps_path=opt.dumps_path, disable=opt.disable_check)