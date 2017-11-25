##!/bin/bash -l
##PBS -l nodes=1:ppn=2,mem=100gb,walltime=1:30:00
##PBS -e error.txt
##PBS -o out.txt

#cd /home/mcgaughs/flag0010/FILET/

#python extract.big.data.set.py

import numpy as np
from random import shuffle
import gzip

def get_gz_file(filename, splitchar = 'NA', buffered = False):
    print(filename)
    if not buffered:
        if splitchar == 'NA':
            return [i.strip().split() for i in gzip.open(filename, 'rt')]
        else: return [i.strip().split(splitchar) for i in gzip.open(filename, 'rt')]
    else:
        if splitchar == 'NA':
            return (i.strip().split() for i in gzip.open(filename, 'rt'))
        else: return (i.strip().split(splitchar) for i in gzip.open(filename, 'rt'))

def load_data(xfile, max_len):
    g = list(get_gz_file(xfile))#'FILET-master/trainingSims/mig21.msOut'
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and i[0].startswith('//')]
    #print(k)
    f = []
    for i in k:
        L = g[i+3:i+37]
        q = []
        for i in L:
            i = [int(j)*2-1 for j in list(i[0])]
            missing = max_len - len(i)
            for z in range(missing): i.append(0)
            i = np.array(i, dtype=np.int8)
            q.append(i)
        q = np.array(q)
        q = q.astype("int8")
        f.append(np.array(q).T)
        #print(q.shape)
        if not len(f) % 100: print len(f)
    return f

def max_len_only(xfile):
    g = list(get_gz_file(xfile))#'FILET-master/trainingSims/mig21.msOut'
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and str(i[0]).startswith('//')]
    #print(k)
    ml = 0
    for i in k:
        L = g[i+3:i+37]
        q = []
        for i in L:
            if len(i[0]) > ml: ml = len(i[0])
    return ml

max_len = max([max_len_only('noMig_100k.msOut.gz'), max_len_only('mig21_100k.msOut.gz')]) #1201

big_mig = load_data('noMig_100k.msOut.gz', max_len)
big_mig.extend(load_data('mig21_100k.msOut.gz', max_len))
big_mig.extend(load_data('trainingSims/noMig.msOut.gz', max_len))
big_mig.extend(load_data('trainingSims/mig12.msOut.gz', max_len))
big_mig.extend(load_data('trainingSims/mig21.msOut.gz', max_len))
x = np.array(big_mig)
del(big_mig)
y = [0 for i in xrange(100000)]
y.extend([2 for i in xrange(100000)])
y.extend([0 for i in xrange(10000)])
y.extend([1 for i in xrange(10000)])
y.extend([2 for i in xrange(10000)])
y = np.array(y)

print len(x), len(y)
shf = range(len(x))
shuffle(shf)

y = y[shf]
x = x[shf]

xtrain, xtest = x[10000:], x[:10000]
ytrain, ytest = y[10000:], y[:10000]
np.savez_compressed('big_sim.npz', xtest=xtest, xtrain=xtrain, ytest=ytest, ytrain=ytrain)

