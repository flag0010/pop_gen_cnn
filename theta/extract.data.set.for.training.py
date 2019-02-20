import numpy as np
from random import shuffle
import gzip
import json
from common import choose

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
    g = list(get_gz_file(xfile))
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and i[0].startswith('//')]
    #print(k)
    f = []
    for i in k:
        L = g[i+3:i+43]
        q = []
        for i in L:
            i = [int(j) for j in list(i[0])]
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

max_len = max([int(i.strip().split()[1]) for i in gzip.open('theta.sims.txt.gz') if 'seg' in i])
print max_len
x = load_data('theta.sims.txt.gz', max_len)
y = json.load(open('true.thetas.json'))

xtrain, xtest = x[1000:], x[:1000]

ytrain, ytest = y[1000:], y[:1000]
np.savez_compressed('theta_sim.npz', xtest=xtest, xtrain=xtrain, ytest=ytest, ytrain=ytrain)#, fin_xtest=fin_xtest, fin_ytest=fin_ytest)

#print x[1].shape

#def pi(x):
#    v = []
#    d = choose(x[0].shape[1], 2) ##
#    for i in x:
#        s = i.sum(axis=1)
#        k = []
#        for j in s:
#            n = choose(j, 2)+choose(i.shape[1]-j, 2)
#            k.append((d-n)*d**-1)
#        v.append(sum(k))
#    return np.array(v)

#def theta(x):
#    v = []
#    a = sum([i**-1 for i in range(1, x[0].shape[1])])
#    print a
#    for i in x:
#        S = i.shape[0]
#        print S
#        v.append(S/a)
#    return np.array(v)

#y = np.array(y)
#pix = pi(x)
#thetax = theta(x)
#print np.mean((pix-y)**2)**.5
#print np.mean((thetax-y)**2)**.5
#ensemble = (pix+thetax)/2.
#print np.mean((ensemble-y)**2)**.5


#from matplotlib import pyplot as plt
#plt.subplot(1,2,1)
#plt.scatter(y, pix, alpha=.4)
#plt.subplot(1,2,2)
#plt.scatter(y, thetax, alpha=.4)
#plt.show()
