import numpy as np
import keras
from keras.preprocessing import sequence
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Merge
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D

num_classes = 5
u = np.load('/dataset/final.split.up.seln.big.npz', encoding = 'latin1')

xtest, ytest,postest = [u[i][:2000] for i in 'xtest ytest postest'.split()]
ytest = keras.utils.to_categorical(ytest, num_classes)
postest = sequence.pad_sequences(postest, padding='post', value=-1., dtype='float32',  maxlen=5000)
xtest = sequence.pad_sequences(xtest, padding='post',  maxlen=5000)

xtrains, ytrains, postrains = [i for i in u.keys() if "xtrain" in i],[i for i in u.keys() if "ytrain" in i], [i for i in u.keys() if "postrain" in i] 
xtrains.sort(), ytrains.sort(), postrains.sort()

training_data = list(zip(*[xtrains, ytrains, postrains]))

batch_size, epochs = 32, 3
ksize = 2
l2_lambda = 0.0001

b1 = Sequential()
b1.add(Conv1D(128*2, kernel_size=ksize,
                 activation='relu',
                 input_shape=( xtest.shape[1], xtest.shape[2]),
                 kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b1.add(MaxPooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b1.add(MaxPooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b1.add(AveragePooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b1.add(AveragePooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Flatten())
b2 = Sequential()
b2.add(Dense(64, input_shape = (5000,), activation='relu',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
b2.add(Dropout(0.1))
    
model = Sequential()
model.add(Merge([b1, b2], mode = 'concat'))
model.add(Dense(256, activation='relu', kernel_initializer='normal',kernel_regularizer=keras.regularizers.l2(l2_lambda)))
model.add(Dropout(0.25))
model.add(Dense(num_classes, activation='softmax'))
    
model.compile(loss=keras.losses.categorical_crossentropy,
                  optimizer=keras.optimizers.Adam(),
                  metrics=['accuracy'])
print(model.summary())

for e in range(epochs):
    for xtrainx, ytrainx, postrainx in training_data:
        print (e, xtrainx, ytrainx, postrainx)
        xtrain, ytrain,postrain = u[xtrainx], u[ytrainx], u[postrainx]
        ytrain = keras.utils.to_categorical(ytrain, num_classes)
        postrain = sequence.pad_sequences(postrain, padding='post', value=-1., dtype='float32',  maxlen=5000)
        xtrain = sequence.pad_sequences(xtrain, padding='post',  maxlen=5000)
        model.fit([xtrain, postrain], ytrain, batch_size=64,
                  epochs=1,
                  verbose=1,
                  validation_data=([xtest, postest], ytest))


del(xtest)
del(xtrain)
del(xtrainx)

posfinal, xfinal, yfinal = [u[i] for i in ['posfinal', 'xfinal', 'yfinal']]
z = [0,0,0,0,0]
keep = []
idx = 0
while sum(z) < 3000:
    if z[yfinal[idx]] < 600: 
        keep.append(idx)
        z[yfinal[idx]]+=1
    idx+=1

posfinal, xfinal, yfinal = posfinal[keep], xfinal[keep], yfinal[keep]

posfinal = sequence.pad_sequences(posfinal, padding='post', value=-1., dtype='float32',  maxlen=5000)
xfinal = sequence.pad_sequences(xfinal, padding='post',  maxlen=5000)

pred = model.predict([xfinal, posfinal])

pred = [np.argmax(i) for i in pred]
from sklearn import metrics
print(metrics.confusion_matrix(yfinal, pred))

j = model.to_json()
print(j)  #need to save this weith weights file to run on validation
model.save_weights('seln.split.merge.mod.weights')

