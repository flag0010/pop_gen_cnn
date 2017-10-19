import keras
model = keras.models.load_model('big.data.89.2.acc.mod')
import numpy as np

u = np.load('mig.runs/big_sim.npz')
pred = model.predict(u['xtest'])

k = []
for idx,i in enumerate(pred):
    k.append(i/sum(i))

n = []
for i,j in zip(k,t):
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
