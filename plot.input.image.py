import numpy as np
from matplotlib import pyplot as plt
from random import *
def transpose_shape(x):
    n = []
    for i in x: n.append(i.T)
    return np.array(n)

def flip_bits(x):  
    n = []
    for idx, i in enumerate(x):
        m = []
        for j in i:
            ch = choice([0,1])
            if ch and j.sum(): #avoid flipping bits on the all zeros padded section (which has sum = 0)
                #print j
                j = j*-1+1
                #print j, '\n'#works
            m.append(j)
        n.append(m)
        if not idx % 250: print len(x) - idx, 'left to flip'
    return np.array(n)

def shuffle_rows(x):
    n = []
    s = set(range(x.shape[1]))
    for kdx, i in enumerate(x):
        out = set([idx for idx, j in enumerate(i) if not j.sum()])
        shuffleable = list(s-out)
        #print sorted(out)
        #print i.shape
        shuffle(shuffleable)
        shuffleable.extend(sorted(out))
        n.append(i[shuffleable])
        if not kdx % 250: print len(x) - kdx, 'left to shuffle'
    return np.array(n)



u = np.load('sim.data.npz')

x = transpose_shape(u['x_slim'])
y = u['y_slim']


plt.subplot(1,2,1)
plt.imshow(x[443], aspect=.10, cmap='gray')

plt.subplot(1,2,2)
plt.imshow(shuffle_rows(x)[443], aspect=.10, cmap='gray')
plt.show()



plt.subplot(1,2,1)
plt.imshow(x[143], aspect=.10, cmap='gray')
plt.title('Originial')

plt.subplot(1,2,2)
plt.title("SNP order shuffled")
plt.imshow(flip_bits(x)[143], aspect=.10, cmap='gray')
plt.show()


idx = 0
s = [0, 0]
for i in x:
    g = i.sum()
    if y[idx] > 0:
        if g > s[1]:  s = [idx, g]
    idx+=1

print y[s[0]]

plt.imshow(x[s[0]], aspect=.10, cmap='gray')
plt.show()


idx = 0
s = [0, 0]
for i in x:
    g = i.sum()
    if y[idx] == 0:
        if g > s[1]:  s = [idx, g]
    idx+=1

print y[s[0]]

plt.imshow(x[s[0]], aspect=.10, cmap='gray')
plt.show()

idx = 0
s = [0, 0]
for i in x:
    g = i.sum()
    if y[idx] == 2:
        if g > s[1]:  s = [idx, g]
    idx+=1

print y[s[0]]

plt.imshow(x[s[0]], aspect=.10, cmap='gray')
plt.show()

