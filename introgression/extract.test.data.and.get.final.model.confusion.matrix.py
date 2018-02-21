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

max_len = max([max_len_only('mig12.test.msOut.gz'), max_len_only('mig21.test.msOut.gz'), max_len_only('noMig.test.msOut.gz'), 1201]) #1201
print max_len

mig = []
y = []

for i,j in zip(('mig12.test.msOut.gz', 'mig21.test.msOut.gz', 'noMig.test.msOut.gz'), (1, 2, 0)):
    l = load_data(i, max_len)
    mig.extend(l)
    y.extend([j for iii in xrange(len(l))])
x = np.array(mig)

print len(x), len(y)


from keras.models import load_model
from sklearn.metrics import confusion_matrix

model = load_model('big.data.89.2.acc.mod')
pred = model.predict(x)
pred_cat = [i.argmax() for i in pred]
print confusion_matrix(y, pred_cat)
print
print confusion_matrix(y, pred_cat) / float(len(y))

k = []
for idx,i in enumerate(pred):
    #k.append(np.exp(i)/sum(np.exp(i)))
    k.append(i/sum(i))

n = []
for i,j in zip(k,y):
    if i.argmax() == j: val = 1
    else: val=0
    prob = i[i.argmax()]
    n.append((prob, val))

d = {}
s,e = 0,.1
for i in range(10):
    if e > .999: e=1
    d[(s,e)] = [0., 0.]
    s+=.1
    e+=.1

for i,o in n:
    for j1,j2 in d:
        if j1 < i <= j2: 
            d[(j1,j2)][o]+=1

for i in sorted(d):
    if sum(d[i]): print i,'\t', d[i][1] / sum(d[i]), '\t', sum(d[i])
