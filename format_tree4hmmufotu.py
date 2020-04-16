from ete3 import Tree
import sys
from tqdm import tqdm

fn = sys.argv[1]
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
