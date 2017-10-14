from keras.preprocessing import sequence
from common import *
from gzip import GzipFile as gz
import numpy as np
from random import shuffle
import h5py
from sklearn.neighbors import NearestNeighbors


def resort_min_diff(amat):
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def load_data(xfile, maxlen = 5000, hard_stop = 1e60):
    g = [ii.strip().split() for ii in gz(xfile)] 
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and  i[0].startswith('//')]
    f = []
    lens = []
    all_pos = []
    for idx, i in enumerate(k): 
        l = g[i+3:i+211]
        pos = np.array(map(float, g[i+2][1:]), dtype="float32")
        #print l[0][:10], l[-1][:10]
        q = []
        for i in l:
            i = list(i[0])
            q.append(np.array([int(j) for j in i], dtype='int8'))
        #print len(q)
        #print q[0][:10], q[-1][:10]
        q = resort_min_diff(np.array(q)).T
        q = (q*-2+1)*-1  # this maps 0/1 to -1/1
        if q.shape[0] <= maxlen:
            f.append(q)
            all_pos.append(pos)
        if len(f) > 49:
            if not len(f) % 50: print idx, len(f), len(all_pos)
        lens.append(len(q[0]))
        if len(f) >= hard_stop: break
    print len([i for i in lens if i > maxlen]), len(lens)
    print '*********', len(f), len(all_pos)
    return f, all_pos



testpath = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/testData/JPT/'
trainpath = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/trainingData/JPT/'


hard = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/testData/JPT/psmcHard_5.msOut.gz'
neut = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/testData/JPT/psmcNeut.msOut.gz'
soft = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/testData/JPT/psmcSoft_5.msOut.gz'

all_pos = []
hard, p = load_data(hard)
all_pos.extend(p)
print len(all_pos), len(hard)
#HARD

neut, p = load_data(neut)
all_pos.extend(p)
print len(all_pos)
#NEUT


soft, p = load_data(soft)
all_pos.extend(p)
print len(all_pos)
#soft


x = [i for i in hard]
x.extend(neut)
x.extend(soft)
y = [1 for i in xrange(len(hard))]
y.extend([0 for i in xrange(len(neut))])
y.extend([2 for i in xrange(len(soft))])

hard = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/trainingData/JPT/psmcHard_5.msOut.gz'
neut = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/trainingData/JPT/psmcNeut.msOut.gz'
soft = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/trainingData/JPT/psmcSoft_5.msOut.gz'

hard, p = load_data(hard)
all_pos.extend(p)
print len(all_pos), len(hard)
#HARD

neut, p = load_data(neut)
all_pos.extend(p)
print len(all_pos)
#NEUT


soft, p = load_data(soft)
all_pos.extend(p)
print len(all_pos)
#soft


x.extend(hard)
x.extend(neut)
x.extend(soft)
y.extend([1 for i in xrange(len(hard))])
y.extend([0 for i in xrange(len(neut))])
y.extend([2 for i in xrange(len(soft))])

x = np.array(x)
y = np.array(y)
all_pos = np.array(all_pos)

print len(x), len(y), len(all_pos)

shf = range(len(x))
shuffle(shf)

x = x[shf]
y = y[shf]

maxsize = 0
for i in x:
    if i.shape[0] > maxsize: maxsize = i.shape[0]

x = sequence.pad_sequences(x, maxlen=maxsize, padding='post')

all_pos = all_pos[shf] 

class_ct = [0,0,0]
xtrain, xtest, ytrain, ytest, postrain, postest = [],[],[],[],[],[]
for idx, i in enumerate(y):
    if class_ct[i] < 351: #a hit
        xtest.append(x[idx])
        ytest.append(y[idx])
        postest.append(all_pos[idx])
        class_ct[i]+=1
    else:
        xtrain.append(x[idx])
        ytrain.append(y[idx])
        postrain.append(all_pos[idx])

xtrain,	xtest, ytrain, ytest, postrain, postest = map(np.array, [xtrain, xtest, ytrain, ytest, postrain, postest])
#xtrain, xtest = x[1200:], x[:1200]
#ytrain, ytest = y[1200:], y[:1200]
#postrain, postest = all_pos[1000:], all_pos[:1000]

np.savez_compressed('hard.soft.neut.npz', xtrain = xtrain, xtest=xtest, ytrain=ytrain, ytest=ytest, postrain=postrain, postest=postest)

f = h5py.File('training.data.h5', 'w')
dset = f.create_dataset('xtrain', xtrain.shape, dtype='i')
dset[...] = xtrain
f.close()
