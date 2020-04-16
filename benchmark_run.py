"""
**Copyright (C) 2019  Jose Sergio Hleap**

This script runs optimize_n_score in a range of taxonomic levels, and plots the
result
"""
import os
import dill
import pandas as pd
import argparse
from optimize_n_score import main as opt_n_score
from glob import glob


def main():
    fas = glob('*.fa')
    dfs = []
    for query in fas:
        pref = os.path.splitext(query)[0]
        print('Processing', pref)
        truthfile = '%s.txt' % pref
        mock = pref.split('_')[0]
        pcklfile = '%s.pckl' % mock
        if os.path.isfile(pcklfile):
            with open(pcklfile, 'rb') as p:
                ds = dill.load(p)
        else:
            ds = opt_n_score(query, truthfile)
            with open(pcklfile, 'wb') as p:
                dill.dump(ds, p)
        dfs.append(ds)
    dfs = pd.concat(dfs)
    dfs.to_csv('Benchmark_results.tsv', sep='\t')


if __name__ == '__main__':
    main()
