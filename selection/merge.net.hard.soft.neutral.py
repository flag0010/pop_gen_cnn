import numpy as np
import keras
from keras.preprocessing import sequence
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Merge
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D

batch_size, epochs, num_classes = 32*2, 20, 3

u = np.load('/dataset/hard.soft.neut.npz')
xtrain, xtest, ytrain, ytest = [u[i] for i in 'xtrain xtest ytrain ytest'.split()]

ytest = keras.utils.to_categorical(ytest, num_classes)
ytrain = keras.utils.to_categorical(ytrain, num_classes)

ksize = 2

u = np.load('/dataset/hard.soft.neut.npz', encoding = 'latin1')
postrain, postest  = [u[i] for i in ['postrain', 'postest']]
postrain = sequence.pad_sequences(postrain, padding='post', value=-1., dtype='float32', maxlen=5000)
postest = sequence.pad_sequences(postest, padding='post', value=-1., dtype='float32',  maxlen=5000)

#xtrain/test are the SNP matrices
#ytrain/test is the categories
#postrain/test are vectors of SNP positions normalized between 0-1

#branch1
#this processes the SNP matrix
b1 = Sequential()
b1.add(Conv1D(128*2, kernel_size=ksize,
                 activation='relu',
                 input_shape=( xtest.shape[1], xtest.shape[2])))
b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
b1.add(MaxPooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
b1.add(MaxPooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
b1.add(AveragePooling1D(pool_size=ksize))
b1.add(Dropout(0.2))

b1.add(Flatten())

#branch2
#this processes the position data
b2 = Sequential()
b2.add(Dense(64, input_shape = (5000,), activation='relu'))
b2.add(Dropout(0.1))

#merging branch 1 and 2
model = Sequential()
model.add(Merge([b1, b2], mode = 'concat'))

model.add(Dense(256, activation='relu', kernel_initializer='normal'))
model.add(Dropout(0.25))

model.add(Dense(num_classes, activation='softmax'))

#compile and fit

model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])
print(model.summary())

#note trick with fitting, you pass branch1/2 data in as list
model.fit([xtrain, postrain], ytrain, batch_size=64*2,
          epochs=epochs,
          verbose=1,
          validation_data=([xtest, postest], ytest)) #same for testing data, make a list of data
          
#confusion matrix
from sklearn.metrics import confusion_matrix
m = [np.argmax(i) for i in model.predict([xtest,postest])]
n = [np.argmax(k) for k in ytest]
confusion_matrix(n,m)/len(n)
