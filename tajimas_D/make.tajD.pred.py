import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from generic_tajD import calc_S_and_k_from_seqs, count_all, tajD
from random import choice
from matplotlib import pyplot as plt
#from keras_diagram import ascii
from sklearn.neighbors import NearestNeighbors

def sort_min_diff(amat):
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def transpose_shape(x):
    n = []
    for i in x: n.append(i.T)
    return np.array(n)

x = []
for i in xrange(5000):
    q = []
    v = np.random.randint(2, size=40)*-2+1 
    for j in range(20):
        if choice([0,1]): q.append(v)
        else: q.append(np.random.randint(2, size=40)*2-1)
    x.append(sort_min_diff(np.array(q)))
    
x = np.array(x)
print x.shape

y = []
for i in x:
    S,k =  calc_S_and_k_from_seqs(i)
    #print S,k
    td = tajD(20, S, k)
    #print td
    y.append(td)

ytest, ytrain = y[:1000], y[1000:]
xtest, xtrain = transpose_shape(x[:1000]), transpose_shape(x[1000:])

model = Sequential()
model.add(Conv1D(128, kernel_size=2,
                 activation='relu',
                 input_shape=(40, 20)))
model.add(Conv1D(128, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Dropout(0.1))
model.add(Flatten())
model.add(Dense(128, activation='relu', kernel_initializer='normal'))
model.add(Dropout(0.25))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss='mean_squared_error', optimizer='adam')
model.summary()
#print(ascii(model))
model.fit(xtrain, ytrain, batch_size=64,
          epochs=7, verbose=1,
          validation_data=(xtest, ytest))

pred = model.predict(xtest)
plt.scatter(ytest, pred, alpha=.5)
plt.plot([-3, 3], [-3, 3], color='k')
plt.xlabel("True Tajima's D")
plt.ylabel("Predicted Tajima's D")
plt.show()
