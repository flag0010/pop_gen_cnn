from __future__ import print_function
import numpy as np
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras import backend as K

batch_size = 200
epochs = 12
num_classes = 3

u = np.load('sim.data.npz')
x_test, x_train, y_test, y_train = [u[i] for i in  ['x_test', 'x_train', 'y_test', 'y_train']]

img_rows, img_cols = x_test.shape[1:3]

y_test = keras.utils.to_categorical(y_test, num_classes)
y_train = keras.utils.to_categorical(y_train, num_classes)

n = []
for i in x_train: n.append(i.T)
x_train = np.array(n)

n = []
for i in x_test: n.append(i.T)
x_test = np.array(n)

epochs = 60  #this one isn't so bad.  the pooling is 2^3, so blocks of 8 indv among 34 lines.  that compresses things
model = Sequential()
model.add(Conv1D(128, kernel_size=2,
                 activation='relu',
                 input_shape=(1102, 34)))
model.add(Conv1D(64, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Conv1D(64, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Conv1D(64, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dense(num_classes, activation='softmax'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adadelta(),
              metrics=['accuracy'])
model.fit(x_train, y_train, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test))
model.save('theano.ave.60.mod')

print(model.evaluate(x_test, y_test))
