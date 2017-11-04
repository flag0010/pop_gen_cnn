import numpy as np
import keras
#from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Merge
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from generic_tajD import calc_S_and_k_from_seqs, count_all, tajD
from random import choice
from matplotlib import pyplot as plt
#from keras_diagram import ascii
from sklearn.neighbors import NearestNeighbors
from pprint import pprint

nreps, nepoch = 40, 10

def resort_min_diff(amat):
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
    #print i
    q = []
    v = np.random.randint(2, size=60)#*-2+1 
    #print v
    for j in range(35):
        if choice([0,1]): q.append(v)
        else: q.append(np.random.randint(2, size=60))
    x.append(np.array(q))
    
x = np.array(x)

#print x[1]
#print x.shape

y = []
for i in x:
    S,k =  calc_S_and_k_from_seqs(i)
    #print S,k
    td = tajD(35, S, k)
    #print td
    y.append(td)

ytest, ytrain = y[:1000], y[1000:]
xtest, xtrain = x[:1000], x[1000:]
#print xtest.shape
#print xtrain.shape

np.savez_compressed('tajd.npz', ytest=ytest, xtest=xtest, xtrain=xtrain, ytrain=ytrain)

all_out = {'not_transposed':[], 'binary':[], 'neg1_1': [], 'resort':[], 'resort_and_neg1_1':[]}
for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(128, kernel_size=2,
                     activation='relu',
                     input_shape=(35, 60)))
    #model.add(Dropout(0.25))#not helpful
    model.add(Conv1D(128, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.1))
    model.add(Flatten())
    model.add(Dense(128, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.25))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    pred = model.predict(xtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
    f = [rmse]
    for j in range(nepoch):    
        model.fit(xtrain, ytrain, batch_size=64,
                  epochs=1, verbose=0, validation_data=(xtest, ytest))
        pred = model.predict(xtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['not_transposed'].append(f)
    pprint( all_out )

xtrain, xtest = transpose_shape(xtrain), transpose_shape(xtest)

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(128, kernel_size=2,
                     activation='relu',
                     input_shape=(60, 35)))
    #model.add(Dropout(0.25))#not helpful
    model.add(Conv1D(128, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.1))
    model.add(Flatten())
    model.add(Dense(128, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.25))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    pred = model.predict(xtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
    f = [rmse]
    for j in range(nepoch):    
        model.fit(xtrain, ytrain, batch_size=64,
                  epochs=1, verbose=0, validation_data=(xtest, ytest))
        pred = model.predict(xtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['binary'].append(f)
    pprint( all_out )

mtrain, mtest = (xtrain*-2+1)*-1, (xtest*-2+1)*-1

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(128, kernel_size=2,
                     activation='relu',
                     input_shape=(60, 35)))
    #model.add(Dropout(0.25))#not helpful
    model.add(Conv1D(128, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.1))
    model.add(Flatten())
    model.add(Dense(128, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.25))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    pred = model.predict(mtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
    f = [rmse]
    for j in range(nepoch):    
        model.fit(mtrain, ytrain, batch_size=64,
                  epochs=1, verbose=0, validation_data=(mtest, ytest))
        pred = model.predict(mtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['neg1_1'].append(f)
    pprint( all_out )

rtrain = np.array([resort_min_diff(i.T).T for i in xtrain])
rtest = np.array([resort_min_diff(i.T).T for i in xtest])

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(128, kernel_size=2,
                     activation='relu',
                     input_shape=(60, 35)))
    #model.add(Dropout(0.25))#not helpful
    model.add(Conv1D(128, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.1))
    model.add(Flatten())
    model.add(Dense(128, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.25))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    pred = model.predict(rtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
    f = [rmse]
    for j in range(nepoch):    
        model.fit(rtrain, ytrain, batch_size=64,
                  epochs=1, verbose=0, validation_data=(rtest, ytest))
        pred = model.predict(rtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['resort'].append(f)
    pprint( all_out )

rtrain = (rtrain*-2+1)*-1
rtest = (rtest*-2+1)*-1

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(128, kernel_size=2,
                     activation='relu',
                     input_shape=(60, 35)))
    #model.add(Dropout(0.25))#not helpful
    model.add(Conv1D(128, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.1))
    model.add(Flatten())
    model.add(Dense(128, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.25))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    pred = model.predict(rtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
    f = [rmse]
    for j in range(nepoch):    
        model.fit(rtrain, ytrain, batch_size=64,
                  epochs=1, verbose=0, validation_data=(rtest, ytest))
        pred = model.predict(rtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(ytest, pred)])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['resort_and_neg1_1'].append(f)
    pprint( all_out )

def rmean(x):
    k = [sum(i)*len(i)**-1 for i in zip(*x)]
    return k

for i,color in zip(['not_transposed', 'binary', 'neg1_1', 'resort', 'resort_and_neg1_1'], ['r', 'b', 'g', 'k', 'magenta']):
    #for j in all_out[i]:
        #plt.plot(range(11), j, color=color, alpha=.05)
        #plt.scatter(range(11), j, color=color, alpha=.3)
    plt.plot(range(11), rmean(all_out[i]), color=color)

from pickle import dump
dump(all_out, open('all.fitted.nets.pickle', 'w'))

plt.xlabel('Training Epoch')
plt.ylabel("Hold-out Data RMSE")
plt.semilogy()
plt.show()
