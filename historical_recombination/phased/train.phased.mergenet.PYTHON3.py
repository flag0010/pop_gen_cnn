import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Merge
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras import backend as K
from random import shuffle, choice
from keras.preprocessing import sequence
from matplotlib import pyplot as plt

def rsquare(x,y):
    return np.corrcoef(x,y)[0][1]**2  #r-squared
def rmse(x,y):
    return np.sqrt(np.mean((x-y)**2))

#load and prep data
a = np.load('../dataset/ld.data.npz' , encoding='latin1')

postrain, ytest, ytrain, xtest, xtrain, postest = [a[i] for i in ['postrain', 'ytest', 'ytrain', 'xtest', 'xtrain', 'postest']]

#pad, max len was prev. found to be 418 (when extracting from ms output)
xtrain = sequence.pad_sequences(xtrain, maxlen=418, padding='post')
xtest = sequence.pad_sequences(xtest, maxlen=418, padding='post')

postrain = sequence.pad_sequences(postrain, maxlen=418, padding='post', value=-1., dtype='float32')
postest = sequence.pad_sequences(postest, maxlen=418, padding='post', value=-1., dtype='float32')

#prep y data
ytrain_rho, ytest_rho = np.array([i[1] for i in ytrain]), np.array([i[1] for i in ytest])
#ytrain_rho_over_theta, ytest_rho_over_theta = np.array([i[1]/i[0] for i in ytrain]), np.array([i[1]/i[0] for i in ytest])

mean_test = np.mean(np.log(ytest_rho))
mean_train = np.mean(np.log(ytrain_rho))

ytrain_rho_log_centered = np.log(ytrain_rho) - mean_train
plt.hist(ytrain_rho_log_centered)
plt.show()

ytest_rho_log_centered = np.log(ytest_rho) - mean_test
plt.hist(ytest_rho_log_centered)
plt.show()
print (mean_test, mean_train) #mean train is important, it was 4.78757605959

batch_size = 32
l2_lambda = 0.0001

b1 = Sequential()
b1.add(Conv1D(1250, kernel_size=2,
                 activation='relu',
                 kernel_regularizer=keras.regularizers.l2(l2_lambda),
                 input_shape=(xtest.shape[1], xtest.shape[2])))
b1.add(Conv1D(256, kernel_size=2, kernel_regularizer=keras.regularizers.l2(l2_lambda), activation='relu'))
b1.add(AveragePooling1D(pool_size=2))
b1.add(Dropout(0.25))
b1.add(Conv1D(256, kernel_size=2, kernel_regularizer=keras.regularizers.l2(l2_lambda), activation='relu'))
b1.add(AveragePooling1D(pool_size=2))
b1.add(Dropout(0.25))
b1.add(Flatten())

b2 = Sequential()
b2.add(Dense(64, input_shape = (418,), kernel_regularizer=keras.regularizers.l2(l2_lambda), activation='relu'))
b2.add(Dropout(0.1))

model = Sequential()
#model.add([b1, b2])
#model = keras.layers.concatenate([b1,b2])
model.add(Merge([b1, b2], mode = 'concat'))
model.add(Dense(256, kernel_regularizer=keras.regularizers.l2(l2_lambda), activation='relu'))
model.add(Dense(1, kernel_initializer='normal'))
model.compile(loss='mean_squared_error', optimizer='adam')

for III in range(25):
    if III < 1: model.summary()

    model.fit([xtrain, postrain], ytrain_rho_log_centered, batch_size=batch_size,
              epochs=1,
              verbose=1,
              validation_data=([xtest, postest], ytest_rho_log_centered))
    pred = [i[0] for i in model.predict([xtest, postest])]
    plt.scatter(np.exp(ytest_rho_log_centered+mean_train), np.exp(pred+mean_train), alpha=.4)
    plt.semilogy()
    plt.semilogx()
    plt.show()
    print (III, rsquare(np.exp(ytest_rho_log_centered+mean_train), np.exp(pred)+mean_train), 
           rmse(np.exp(ytest_rho_log_centered+mean_train), np.exp(pred+mean_train)))
    if rmse(np.exp(ytest_rho_log_centered+mean_train), np.exp(pred+mean_train)) < 210: break
 

j = model.to_json()
print(j)  #need to save this weith weights file to run on validation
model.save_weights('regularized.merge.mod.weights')
