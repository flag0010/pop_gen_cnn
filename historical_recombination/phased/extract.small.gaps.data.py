import numpy as np
from common import *
from gzip import GzipFile as gz
from sklearn.neighbors import NearestNeighbors


def sort_min_diff(amat):
    '''this function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted by genetic similarity.
    this problem is NP-hard, so here we use a nearest neighbors approx.  it's not perfect, but it's fast and generally performs ok.
    assumes your input matrix is a numpy array'''
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def convert_01_to_neg1_1(amat):
    '''convert standard binary 0/1 ms SNP matrix to -1/1 SNP matrix. B/c weights & biases are drawn from a distribution with mean=0
    choosing -1/1 (which is also mean=0) tends to help in training. assumes your input matrix is a numpy array'''
    return (amat*-2+1)*-1 

a = (i.strip().split('asdfasdfd') for i in gz('gap.LD.sims.txt.gz'))

idx = 0
xvals, yvals = {},{}
seg_sites = {}
pos = {}
n = []
ctr = 0
true_idx = 0

for i in a:
    if 'segsites' in i[0]:
        seg_sites[true_idx] = int(i[0].split()[-1])
    if 'positions' in i[0]:
        pos[true_idx] = [float(jj) for jj in i[0].split()[1:]]
    if './ms' in i[0]:
        true_idx+=1
        idx = 0
        theta_line = i[0].split()
        ms_vals = map(float, [theta_line[4], theta_line[6]])
        yvals[true_idx] = np.array(ms_vals)
    if 5 < idx < 56:
        n.append(np.array(map(int, i[0]), dtype='int8'))
    if len(n) == 50:
        n = np.array(n)
        n = sort_min_diff(n)
        n = convert_01_to_neg1_1(n)
        xvals[true_idx] =  n.T 
        idx = 0
        n = []
        if not true_idx % 100: print true_idx, len(xvals), len(yvals), min(xvals), max(xvals), min(yvals), max(yvals) #xvals[max(xvals)].shape, seg_sites[max(seg_sites)]
        #if len(xvals) != len(yvals):
        #    print true_idx, idx, len(xvals), len(yvals), xvals[-2].shape, yvals[-2], '*********'
    #print idx
    idx+=1
    #print xvals, yvals, seg_sites, pos


newx, newy, newpos = [],[],[]
for i in sorted(seg_sites):
    if i in xvals and i in yvals and i in pos:
        s = seg_sites[i]
        if s > 0:
            newpos.append(np.array(pos[i]))
            newx.append(xvals[i])
            newy.append(yvals[i])
max_size = 0
for i in newx:
    if i.shape[0] > max_size: max_size = i.shape[0]

print "max size", max_size
np.savez_compressed('gap.ld.data.npz', xtest=newx, ytest = newy, postest=newpos)
