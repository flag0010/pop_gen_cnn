import numpy as np
import keras
from keras.preprocessing import sequence
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Merge
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D

epochs, num_classes = 20, 3

u = np.load('/dataset/hard.soft.neut.npz')
xtrain, xtest, ytrain, ytest = [u[i] for i in 'xtrain xtest ytrain ytest'.split()]

ytest = keras.utils.to_categorical(ytest, num_classes)
ytrain = keras.utils.to_categorical(ytrain, num_classes)

ksize = 2

model = Sequential()
model.add(Conv1D(128*2, kernel_size=ksize,
                 activation='relu',
                 input_shape=( xtest.shape[1], xtest.shape[2])))
model.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
model.add(MaxPooling1D(pool_size=ksize))
model.add(Dropout(0.2))

model.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
model.add(MaxPooling1D(pool_size=ksize))
model.add(Dropout(0.2))

model.add(Conv1D(128*2, kernel_size=ksize, activation='relu'))
model.add(AveragePooling1D(pool_size=ksize))
model.add(Dropout(0.2))

#model.add(Conv1D(128, kernel_size=ksize, activation='relu'))
#model.add(AveragePooling1D(pool_size=ksize))
#model.add(Dropout(0.1))

model.add(Flatten())

model.add(Dense(128*2, activation='relu'))
#model.add(Dropout(0.5))

model.add(Dense(num_classes, activation='softmax'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])
print(model.summary())
model.fit(xtrain, ytrain, batch_size=64,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest))
