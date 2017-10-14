#run on p2.xlarge ec2 instance.  requires ~40 GB or ram and a GPU
import numpy as np
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras import backend as K
from random import choice, shuffle

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
        if not idx % 1000: print(len(x) - idx, 'left to flip')
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
        if not kdx % 1000: print(len(x) - kdx, 'left to shuffle')
    return np.array(n)

batch_size = 200
epochs = 120
num_classes = 3

u = np.load('sim.data.npz')
x_test, x_train, y_test, y_train = [u[i] for i in  ['x_test', 'x_train', 'y_test', 'y_train']]

img_rows, img_cols = x_test.shape[1:3]

y_test = keras.utils.to_categorical(y_test, num_classes)
y_train = keras.utils.to_categorical(y_train, num_classes)

x_test, x_train = transpose_shape(x_test), transpose_shape(x_train)

s = (shuffle_rows(i) for i in [x_train for ff in range(7)] )
x_train = [i for i in x_train]
for i in s:
    x_train.extend(i)
del(s)
x_train = np.array(x_train)
print(x_train.shape)

print(x_train.shape)
y_train = u["y_train"]
y_train_full = [i for i in y_train]
for i in range(7):
    y_train_full.extend(y_train)
y_train = np.array(y_train_full)
y_train = keras.utils.to_categorical(y_train, num_classes)



[i.shape for i in [x_test, y_test, x_train, y_train]]


#np.savez_compressed('big.sim.data.npz', x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test)

model = Sequential()
model.add(Conv1D(128, kernel_size=2,
                 activation='relu',
                 #input_shape=(34, 1102))) ##tried this, it's worse as expected
                 input_shape=(1102, 34)))
model.add(Conv1D(64, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Dropout(0.25))
model.add(Conv1D(64, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(128, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(num_classes, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])
model.summary()
model.fit(x_train, y_train, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test))



model.save('tf.smaller.87.8.acc.mod')


from sklearn.metrics import confusion_matrix
m = model.predict(x_test)
m = [i for i in map(np.argmax, m)]
n = [i for i in map(np.argmax, y_test)]
print(100*(confusion_matrix(n,m)/len(m)))



acc = 29.516666+33.21666667+25.01666667
