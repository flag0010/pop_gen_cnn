import numpy as np
from random import shuffle

def get_file(filename, splitchar = 'NA', buffered = False):
    print(filename)
    if not buffered:
        if splitchar == 'NA':
            return [i.strip().split() for i in open(filename)]
        else: return [i.strip().split(splitchar) for i in open(filename)]
    else:
        if splitchar == 'NA':
            return (i.strip().split() for i in open(filename))
        else: return (i.strip().split(splitchar) for i in open(filename))

def load_data(xfile, max_len):
    g = get_file(xfile)#'FILET-master/trainingSims/mig21.msOut'
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and  i[0].startswith('//')]
    f = []
    for i in k:
        l = g[i+3:i+37]
        q = []
        for i in l:
            i = [int(j) for j in list(i[0])]
            missing = max_len - len(i)
            for z in range(missing): i.append(0)
            i = np.array(i, dtype=np.float32)
            q.append(i)
        q = np.array(q)
        q = q.astype("float32")
        f.append(np.array(q))
        print q.shape
    return f

def max_len_only(xfile):
    g = get_file(xfile)#'FILET-master/trainingSims/mig21.msOut'
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and  i[0].startswith('//')]
    ml = 0
    for i in k:
        l = g[i+3:i+37]
        q = []
	for i in l:
            if len(i[0]) > ml: ml = len(i[0])
    return ml

max_len = max(map(max_len_only, ['trainingSims/noMig.msOut', 'trainingSims/mig12.msOut', 'trainingSims/mig21.msOut']))
print max_len
nomig = load_data('trainingSims/noMig.msOut', max_len)
mig12= load_data('trainingSims/mig12.msOut', max_len)
mig21= load_data('trainingSims/mig21.msOut', max_len)

y = [0 for i in range(10000)]
y.extend([1 for i in range(10000)])
y.extend([2 for i in range(10000)])
y = np.array(y)

#m = 0
#for i in [nomig,mig21, mig12]:
#    for j in i:
#        if len(i[1]) > m: m=len(i[1])

x = nomig+mig12+mig21
x = np.array(x)

#print x
print x.shape

zidx = range(len(x))
shuffle(zidx)

y = y[zidx]
x = x[zidx]

x_train, y_train = x[:24000], y[:24000]
x_test, y_test = x[24000:], y[24000:]

#from cPickle import dump
#dump([x_train, y_train, x_test, y_test], open('sim.data.pickle', 'w'))

np.savez_compressed('sim.data.npz', x_train=x_train, y_train=y_train, x_test=x_test, 
                                    y_test=y_test, x_slim=x_test[:500], y_slim=y_test[:500])










# def find_max_size(data_list):
#     max_len = 0
#     for i in data_list:
#         for j in i:
#             l = len(j[0])
#             if l > max_len: max_len = l
# max_size = find_max_size([mig12, mig21, nomig])
# max_size
# 
# 
# 
# from __future__ import print_function
# import keras
# from keras.datasets import mnist
# from keras.models import Sequential
# from keras.layers import Dense, Dropout, Flatten
# from keras.layers import Conv2D, MaxPooling2D
# from keras import backend as K
# 
# batch_size = 128
# num_classes = 10
# epochs = 12
# 
# # input image dimensions
# img_rows, img_cols = 28, 28
# 
# # the data, shuffled and split between train and test sets
# (x_train, y_train), (x_test, y_test) = mnist.load_data()
# 
# if K.image_data_format() == 'channels_first':
#     x_train = x_train.reshape(x_train.shape[0], 1, img_rows, img_cols)
#     x_test = x_test.reshape(x_test.shape[0], 1, img_rows, img_cols)
#     input_shape = (1, img_rows, img_cols)
# else:
#     x_train = x_train.reshape(x_train.shape[0], img_rows, img_cols, 1)
#     x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
#     input_shape = (img_rows, img_cols, 1)
# 
# x_train = x_train.astype('float32')
# x_test = x_test.astype('float32')
# x_train /= 255
# x_test /= 255
# print('x_train shape:', x_train.shape)
# print(x_train.shape[0], 'train samples')
# print(x_test.shape[0], 'test samples')
# 
# # convert class vectors to binary class matrices
# y_train = keras.utils.to_categorical(y_train, num_classes)
# y_test = keras.utils.to_categorical(y_test, num_classes)
# 
# model = Sequential()
# model.add(Conv2D(32, kernel_size=(3, 3),
#                  activation='relu',
#                  input_shape=input_shape))
# model.add(Conv2D(64, (3, 3), activation='relu'))
# model.add(MaxPooling2D(pool_size=(2, 2)))
# model.add(Dropout(0.25))
# model.add(Flatten())
# model.add(Dense(128, activation='relu'))
# model.add(Dropout(0.5))
# model.add(Dense(num_classes, activation='softmax'))
# 
# model.compile(loss=keras.losses.categorical_crossentropy,
#               optimizer=keras.optimizers.Adadelta(),
#               metrics=['accuracy'])
# 
# model.fit(x_train, y_train,
#           batch_size=batch_size,
#           epochs=epochs,
#           verbose=1,
#           validation_data=(x_test, y_test))
# score = model.evaluate(x_test, y_test, verbose=0)
# print('Test loss:', score[0])
# print('Test accuracy:', score[1])

