import cPickle
from matplotlib import pyplot as plt
from math import log
import numpy as np
from collections import deque
from statsmodels.nonparametric.smoothers_lowess import lowess
from common import *
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

p = get_file('ldhat.data/pairs.txt', ',')
X = {}
for i in p:
    z = i[0].split('_')
    theta = float(z[0].replace('sites/', ''))
    rho = float(z[1])
    X[rho] = theta

rmse = lambda x,y: np.mean([(float(iii)-float(jjj))**2. for iii,jjj in zip(x,y)])**0.5
r2   = lambda x,y: np.corrcoef(map(float, x),map(float,y))[0][1]**2

a = cPickle.load(open('realrho.ldhatrho.nnrho.pickle'))

nn, ld, real = [],[],[]

for i in a:                                                                           
    real.extend([j[0]/20000. for j in a[i]])                              
    nn.extend([j[2]/20000. for j in a[i]])
    ld.extend([j[1]/20000. for j in a[i]])

missing = []
for i in X:
    if i not in real:
        print i, X[i]
        missing.append(X[i])
h = []
for i in a:
    for j in a[i]: h.append(i)
#plt.hist(h)
#plt.hist(missing)
#plt.show()
print count_all(h)
print count_all(missing)

print "range", max(real), min(real)
print "NN rsquare", r2(real, nn), "rmse",  rmse(real,  nn)
print "LDHAT rsquare", r2(real, ld), 'rmse', rmse(real, ld)
print rmse(real, nn), rmse(real, ld)
print map(len, [real, nn, ld])
print "LDHAT vs NN",  r2(ld, nn), rmse(ld, nn)

zz = np.array([(i+j)/2. for i,j in zip(nn, ld)])
print 'ensemble', r2(real, zz), rmse(real,zz)

nn_better = []

for i,j,k in zip( real, ld, nn):             
    nnm = (i-k)**2
    ldm = (i-j)**2                 
    if nnm == ldm: print 'tie!'
    if nnm < ldm: nn_better.append(1.)
    else: nn_better.append(0.)

print 'proportion nn better', np.mean(nn_better)
b = zip(real, nn_better)
b.sort()

x = deque(maxlen=100)
meanx = []

for i,j in b:                                                                                            
    x.append(j)                                                                                          
    meanx.append(sum(x)*len(x)**-1)

plt.subplot(1,3,1)
plt.plot([0, max(real)], [0, max(real)], color='r')
plt.scatter(real, ld, alpha=.1, color='k', marker='.')
plt.ylabel('LDhat rho estimate')
plt.xlabel('real rho value')

plt.subplot(1,3,2)
plt.plot([0, max(real)], [0, max(real)], color='r')
plt.scatter(real, nn, alpha=.1, color='k', marker='.')
plt.ylabel('CNN rho estimate')
plt.xlabel('real rho value')

plt.subplot(1,3,3)
#plt.scatter([i[0] for i in b], [i[1] for i in b], alpha=.1, color='k', marker='.')
plt.axhline(y=1/2., color='r')
low = lowess(endog=[i[1] for i in b], exog=[i[0] for i in b], frac=.15)
plt.plot(low[:,0], low[:,1], color='k')
plt.ylabel('probability CNN better than LDhat')
plt.xlabel('real rho value')
plt.ylim(0,1)
plt.show()
#print low