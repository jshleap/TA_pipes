from ete3 import NCBITaxa
import sys
ncbi = NCBITaxa()
levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
          'species']


def name_translator(sp):
    d = ncbi.get_name_translator([sp])
    if not d:
        gn = sp.split()[0]
        d = ncbi.get_rank(ncbi.get_lineage(ncbi.get_name_translator([gn])[gn][0]))
        d = {v:k for k, v in d.items()}
        d['species'] = sp
    else:
        d = ncbi.get_rank(ncbi.get_lineage(ncbi.get_name_translator([sp])[sp][0]))
        d = {v: k for k, v in d.items()}
    taxids = [d[x] if x in d else '' for x in levels]
    return taxids


def get_lineage(sp):
    taxids = name_translator(sp)
    names = []
    for tid in taxids:
        if isinstance(tid, str):
            names += [tid]
        else:
            names += ncbi.translate_to_names([tid])
    return ';'.join(names)


for line in open(sys.argv[1]):
    if line.startswith('>') and ';;' in line:
        line = line.replace(';\n', '')
        bl = line.strip().split(';')
        acc = bl[0].split()[0]
        sp = bl[-1]
        lineage = get_lineage(sp)
        print('%s %s;' % (acc, lineage))
    else:
        print(line.strip())


