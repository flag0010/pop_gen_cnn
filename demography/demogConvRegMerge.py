import sys, os
import numpy as np
import keras
from keras.models import Model
from keras.layers import Input, Dense, Dropout, Flatten
from keras.layers.merge import concatenate
from keras.layers.convolutional import Conv2D, Conv1D
from keras.layers.pooling import MaxPooling2D, AveragePooling1D
from keras import backend as K
from sklearn.neighbors import NearestNeighbors

batch_size = 200
epochs = 10

convDim, convSize, poolSize, useLog, useInt, sortRows, useDropout, lossThreshold, inDir, weightFileName, modFileName, testPredFileName = sys.argv[1:]
convDim = convDim.lower()
convSize, poolSize = int(convSize), int(poolSize)
useLog = True if useLog.lower() in ["true","logy"] else False
useInt = True if useInt.lower() in ["true","intallele"] else False
sortRows = True if sortRows.lower() in ["true","sortrows"] else False
useDropout = True if useDropout.lower() in ["true","dropout"] else False
lossThreshold = float(lossThreshold)

def resort_min_diff(amat):
    ###assumes your snp matrix is indv. on rows, snps on cols
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

X = []
y = []
print("reading data")
for npzFileName in os.listdir(inDir):
    u = np.load(inDir + npzFileName)
    currX, curry = [u[i] for i in  ('X', 'y')]
    ni,nr,nc = currX.shape
    newCurrX = []
    for i in range(ni):
        currCurrX = [currX[i,0]]
        if sortRows:
            currCurrX.extend(resort_min_diff(currX[i,1:]))
        else:
            currCurrX.extend(currX[i,1:])
        currCurrX = np.array(currCurrX)
        newCurrX.append(currCurrX.T)
    currX = np.array(newCurrX)
    assert currX.shape == (ni,nc,nr)
    #indices = [i for i in range(nc) if i % 10 == 0]
    #X.extend(np.take(currX,indices,axis=1))
    X.extend(currX)
    y.extend(curry)
    #if len(y) == 10000:
    #    break

y = np.array(y)
numParams=y.shape[1]
if useLog:
    y[y == 0] = 1e-6#give zeros a very tiny value so they don't break our log scaling
    y = np.log(y)
totInstances = len(X)
testSize=10000
valSize=10000
print("formatting data arrays")
X = np.array(X)
posX=X[:,:,0]
X=X[:,:,1:]
imgRows, imgCols = X.shape[1:]

if useInt:
    X = X.astype('int8') 
    #this ^^^ is a bug. the intent was to converts to 0/255, but as it it converts to 0/-1
    #see: https://github.com/flag0010/pop_gen_cnn/issues/4
    #leaving as is so the code represents what we actually did. Bugs and all
    #if you want to fix to get 0/255, please see the suggestion at the link above
else:
    X = X.astype('float32')/127.5-1
if convDim == "2d":
    X = X.reshape(X.shape[0], imgRows, imgCols, 1).astype('float32')
posX = posX.astype('float32')/127.5-1

assert totInstances > testSize+valSize

testy=y[:testSize]
valy=y[testSize:testSize+valSize]
y=y[testSize+valSize:]
testX=X[:testSize]
testPosX=posX[:testSize]
valX=X[testSize:testSize+valSize]
valPosX=posX[testSize:testSize+valSize]
X=X[testSize+valSize:]
posX=posX[testSize+valSize:]

yMeans=np.mean(y, axis=0)
yStds=np.std(y, axis=0)
y = (y-yMeans)/yStds
testy = (testy-yMeans)/yStds
valy = (valy-yMeans)/yStds

print(len(X), len(y), len(yMeans))
print(yMeans, yStds)
print(X.shape, testX.shape, valX.shape)
print(posX.shape, testPosX.shape, valPosX.shape)
print(y.shape, valy.shape)
print("ready to learn (%d params, %d training examples, %d rows, %d cols)" %(numParams, len(X), imgRows, imgCols))

if convDim == "2d":
    inputShape = (imgRows, imgCols, 1)
    convFunc = Conv2D
    poolFunc = MaxPooling2D
else:
    inputShape = (imgRows, imgCols)
    convFunc = Conv1D
    poolFunc = AveragePooling1D

b1 = Input(shape=inputShape)
conv11 = convFunc(128, kernel_size=convSize, activation='relu')(b1)
pool11 = poolFunc(pool_size=poolSize)(conv11)
if useDropout:
    pool11 = Dropout(0.25)(pool11)
conv12 = convFunc(128, kernel_size=2, activation='relu')(pool11)
pool12 = poolFunc(pool_size=poolSize)(conv12)
if useDropout:
    pool12 = Dropout(0.25)(pool12)
conv13 = convFunc(128, kernel_size=2, activation='relu')(pool12)
pool13 = poolFunc(pool_size=poolSize)(conv13)
if useDropout:
    pool13 = Dropout(0.25)(pool13)
conv14 = convFunc(128, kernel_size=2, activation='relu')(pool13)
pool14 = poolFunc(pool_size=poolSize)(conv14)
if useDropout:
    pool14 = Dropout(0.25)(pool14)
flat11 = Flatten()(pool14)

b2 = Input(shape=(imgRows,))
dense21 = Dense(32, activation='relu')(b2)
if useDropout:
    dense21 = Dropout(0.25)(dense21)

merged = concatenate([flat11, dense21])
denseMerged = Dense(256, activation='relu', kernel_initializer='normal')(merged)
if useDropout:
    denseMerged = Dropout(0.25)(denseMerged)
denseOutput = Dense(numParams)(denseMerged)
model = Model(inputs=[b1, b2], outputs=denseOutput)
print(model.summary())

model.compile(loss='mean_squared_error', optimizer='adam')
earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=0, mode='auto')
checkpoint = keras.callbacks.ModelCheckpoint(weightFileName, monitor='val_loss', verbose=1, save_best_only=True, mode='min')
callbacks = [earlystop, checkpoint]

model.fit([X, posX], y, batch_size=batch_size, epochs=epochs, verbose=1, validation_data=([valX, valPosX], valy), callbacks=callbacks)

#now we load the weights for the best-performing model
model.load_weights(weightFileName)
model.compile(loss='mean_squared_error', optimizer='adam')

#now we get the loss for our best model on the test set and emit predictions
testLoss = model.evaluate([testX, testPosX], testy)
print(testLoss)
preds = model.predict([testX, testPosX])
with open(testPredFileName, "w") as outFile:
    for i in range(len(preds)):
        outStr = []
        for j in range(len(preds[i])):
            outStr.append("%f vs %f" %(testy[i][j], preds[i][j]))
        outFile.write("\t".join(outStr) + "\n")

#if the loss is lower than our threshold we save the model file if desired
if modFileName.lower() != "nomod" and testLoss <= lossThreshold:
    with open(modFileName, "w") as modFile:
        modFile.write(model.to_json())
else:
    os.system("rm %s" %(weightFileName))
