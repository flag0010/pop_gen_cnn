import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from random import choice
from matplotlib import pyplot as plt
from sklearn.neighbors import NearestNeighbors
from pprint import pprint

nreps, nepoch = 10, 10

def sort_min_diff(amat):
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def transpose_shape(x):
    n = []
    for i in x: n.append(i.T)
    return np.array(n)

def standardize(y):
    return y.mean(), y.std(), (y-y.mean())/y.std()

def unstandardize(y, mean, std):
    return y*std+mean

a = np.load('theta_sim.npz')
xtrain, xtest, ytest, ytrain = [a[i] for i in ['xtrain', 'xtest', 'ytest', 'ytrain']]
n_seg_sites, n_indv = xtrain[1].shape

xtest_untransposed, xtrain_untransposed = map(transpose_shape, [xtest, xtrain])
ytest_mean, ytest_std, ytest = standardize(ytest)
ytrain_mean, ytrain_std, ytrain = standardize(ytrain)
print(ytest_mean, ytest_std)

all_out = {'not_transposed':[], 'transpose':[], 'sort_and_transpose':[]}
 
for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(64, kernel_size=2,
                     activation='relu',
                     input_shape=(n_indv, n_seg_sites)))
    model.add(Conv1D(64, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    #model.summary()
    pred = model.predict(xtest_untransposed)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
    f = [rmse]
    print(i, f[0])
    for j in range(nepoch):    
        model.fit(xtrain_untransposed, ytrain, batch_size=32,
                  epochs=1, verbose=0, validation_data=(xtest_untransposed, ytest))
        pred = model.predict(xtest_untransposed)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['not_transposed'].append(f)
pprint( all_out )

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(64, kernel_size=2,
                     activation='relu',
                     input_shape=(n_seg_sites, n_indv)))
    model.add(Conv1D(64, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    #model.summary()
    pred = model.predict(xtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
    f = [rmse]
    print(i, f[0])
    for j in range(nepoch):    
        model.fit(xtrain, ytrain, batch_size=32,
                  epochs=1, verbose=0, validation_data=(xtest, ytest))
        pred = model.predict(xtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['transpose'].append(f)
pprint( all_out )


rtrain = np.array([sort_min_diff(i.T).T for i in xtrain])
rtest = np.array([sort_min_diff(i.T).T for i in xtest])

for i in range(nreps):
    model = Sequential()
    model.add(Conv1D(64, kernel_size=2,
                     activation='relu',
                     input_shape=(n_seg_sites, n_indv)))
    model.add(Conv1D(64, kernel_size=2, activation='relu'))
    model.add(AveragePooling1D(pool_size=2))
    model.add(Dropout(0.25))
    model.add(Flatten())
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(64, activation='relu', kernel_initializer='normal'))
    model.add(Dropout(0.5))
    model.add(Dense(1, kernel_initializer='normal'))
    model.compile(loss='mean_squared_error', optimizer='adam')
    #model.summary()
    pred = model.predict(rtest)
    rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
    f = [rmse]
    print(i, f[0])
    for j in range(nepoch):    
        model.fit(rtrain, ytrain, batch_size=32,
                  epochs=1, verbose=0, validation_data=(rtest, ytest))
        pred = model.predict(rtest)
        rmse = np.mean([(iii-jjj)**2 for iii,jjj in zip(unstandardize(ytest, ytest_mean, ytest_std), unstandardize(pred, ytest_mean, ytest_std))])**0.5
        print( i,j, rmse)
        f.append(rmse)
    all_out['sort_and_transpose'].append(f)
pprint( all_out )

