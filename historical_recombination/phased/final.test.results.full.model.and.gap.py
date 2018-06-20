import keras
from keras.preprocessing import sequence
import numpy as np
from matplotlib import pyplot as plt
import sys
import json
from scipy import stats
spearman = lambda a,b: stats.spearmanr(a,b)[0]
from sklearn.neighbors import NearestNeighbors


#model magic below. just a json string that specifies model structure
json_str = '{"keras_version": "2.0.6", "class_name": "Sequential", "config": [{"class_name": "Merge", "config": {"mode_type": "raw", "output_mask_type": "raw", "mode": "concat", "dot_axes": -1, "arguments": {}, "name": "merge_2", "output_shape": null, "concat_axis": -1, "output_mask": null, "output_shape_type": "raw", "layers": [{"class_name": "Sequential", "config": [{"class_name": "Conv1D", "config": {"padding": "valid", "use_bias": true, "activity_regularizer": null, "kernel_size": [2], "kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "distribution": "uniform", "seed": null}}, "dtype": "float32", "filters": 1250, "batch_input_shape": [null, 418, 50], "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "dilation_rate": [1], "strides": [1], "trainable": true, "name": "conv1d_4", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "relu", "bias_initializer": {"class_name": "Zeros", "config": {}}}}, {"class_name": "Conv1D", "config": {"padding": "valid", "use_bias": true, "activity_regularizer": null, "kernel_size": [2], "kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "distribution": "uniform", "seed": null}}, "filters": 256, "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "dilation_rate": [1], "strides": [1], "trainable": true, "name": "conv1d_5", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "relu", "bias_initializer": {"class_name": "Zeros", "config": {}}}}, {"class_name": "AveragePooling1D", "config": {"padding": "valid", "strides": [2], "pool_size": [2], "trainable": true, "name": "average_pooling1d_3"}}, {"class_name": "Dropout", "config": {"rate": 0.25, "trainable": true, "name": "dropout_4"}}, {"class_name": "Conv1D", "config": {"padding": "valid", "use_bias": true, "activity_regularizer": null, "kernel_size": [2], "kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "distribution": "uniform", "seed": null}}, "filters": 256, "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "dilation_rate": [1], "strides": [1], "trainable": true, "name": "conv1d_6", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "relu", "bias_initializer": {"class_name": "Zeros", "config": {}}}}, {"class_name": "AveragePooling1D", "config": {"padding": "valid", "strides": [2], "pool_size": [2], "trainable": true, "name": "average_pooling1d_4"}}, {"class_name": "Dropout", "config": {"rate": 0.25, "trainable": true, "name": "dropout_5"}}, {"class_name": "Flatten", "config": {"trainable": true, "name": "flatten_2"}}]}, {"class_name": "Sequential", "config": [{"class_name": "Dense", "config": {"units": 64, "use_bias": true, "activity_regularizer": null, "kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "distribution": "uniform", "seed": null}}, "dtype": "float32", "batch_input_shape": [null, 418], "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "trainable": true, "name": "dense_4", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "relu", "bias_initializer": {"class_name": "Zeros", "config": {}}}}, {"class_name": "Dropout", "config": {"rate": 0.1, "trainable": true, "name": "dropout_6"}}]}]}}, {"class_name": "Dense", "config": {"units": 256, "use_bias": true, "activity_regularizer": null, "kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "distribution": "uniform", "seed": null}}, "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "trainable": true, "name": "dense_5", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "relu", "bias_initializer": {"class_name": "Zeros", "config": {}}}}, {"class_name": "Dense", "config": {"units": 1, "use_bias": true, "activity_regularizer": null, "kernel_initializer": {"class_name": "RandomNormal", "config": {"stddev": 0.05, "seed": null, "mean": 0.0}}, "kernel_regularizer": null, "trainable": true, "name": "dense_6", "kernel_constraint": null, "bias_regularizer": null, "bias_constraint": null, "activation": "linear", "bias_initializer": {"class_name": "Zeros", "config": {}}}}], "backend": "tensorflow"}'


def sort_min_diff(amat):
    '''this function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted by genetic similarity.
    this problem is NP-hard, so here we use a nearest neighbors approx.  it's not perfect, but it's fast and generally performs ok.
    assumes your input matrix is a numpy array'''
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def rsquare(x,y):
    return np.corrcoef(x,y)[0][1]**2  #r-squared

def rmse(x,y):
    return np.sqrt(np.mean((x-y)**2))

s = json.load(open('ldhat.data/validate.data.LD.json'))

xtest = []
postest = []
ytest = []
maxL = 0
idx = 0
for pos, theta_rho, str_mat in s:
    postest.append(np.array(pos))
    nm = []
    for i in str_mat:
        nm.append([float(j) for j in i])
    nm = sort_min_diff(np.array(nm)).T
    #print nm.shape
    xtest.append(nm)
    ytest.append(theta_rho)
    if nm.shape[1] > maxL: maxL = nm.shape[1]
    if not idx % 100: print idx
    idx+=1
print maxL

postest = np.array(postest)
xtest = np.array(xtest)

# xtemp = []
# for i in xtest:
#     v = sequence.pad_sequences(i, maxlen=418, padding='post')
#     xtemp.append(v)
#     
# xtest = np.array(xtemp)
xtest = sequence.pad_sequences(xtest, maxlen=418, padding='post')
postest = sequence.pad_sequences(postest, maxlen=418, padding='post', value=-1., dtype='float32')
ytest_rho = np.array([i[1] for i in ytest])
thetas=[i[0] for i in ytest]

#print postest

mean_test = 4.78757605959  #use training data mean ln(rho)   np.mean(np.log(ytest_rho))
print mean_test

mod = keras.models.model_from_json(json_str)
#mod.load_weights('merge.mod.weights')
mod.load_weights('regularized.merge.mod.weights')

pred = mod.predict([xtest,postest])
#plt.hist(pred, bins=24)
#plt.show()

pred = np.exp([i[0]+mean_test for i in pred])

print map(len, (thetas, pred, ytest_rho))
print 'r-squared', rsquare(pred, ytest_rho), 'rmse', rmse(pred, ytest_rho), 'spearman', spearman(pred, ytest_rho)
q = {}
outfile = open('test.data.results.csv', 'w')
outfile.write('index,theta,num_sites,predict_rho,real_rho\n')
idx = 1
for i,j,k, p in zip(thetas, pred, ytest_rho, postest):
    if i not in q: q[i] = []
    q[i].append((j,k))
    num_sites = len([iiii for iiii in p if iiii >= 0])
    outfile.write(','.join(map(str, [idx, i, num_sites, j, k]))+'\n')
    idx+=1
for i in q: print i, len(q[i])
outfile.close()
# 
# idx=1
# for i in sorted(q):
#     plt.subplot(2,3,idx)
#     pred, real = [n[0] for n in q[i]], [n[1] for n in q[i]]
#     plt.scatter(real, pred, alpha=.3)
#     r = rsquare(np.array(real), np.array(pred))
#     plt.title('theta = '+str(i)+' r2:'+str(round(r,2)))
#     idx+=1
# plt.show()
# 
# 
# d = []
# for i in sorted(q):
#     resid = [p-r for p,r in q[i]]
#     d.append(resid)
# plt.boxplot(d)
# plt.xticks(range(1,6), map(str, sorted(q)))
# plt.show()


#now try same model at interpolation in the big gap in Ne we created

a = np.load('gap.ld.data.npz')
ytest, xtest, postest = [a[i] for i in [ 'ytest', 'xtest', 'postest']]
xtest = sequence.pad_sequences(xtest, maxlen=418, padding='post')
postest = sequence.pad_sequences(postest, maxlen=418, padding='post', value=-1., dtype='float32')
ytest_rho = np.array([i[1] for i in ytest])
thetas=[i[0] for i in ytest]
mean_test = 4.78757605959  #use training data mean ln(rho)  np.mean(np.log(ytest_rho))
print mean_test

pred = mod.predict([xtest,postest])
#plt.hist(pred, bins=24)
#plt.show()

pred = np.exp([i[0]+mean_test for i in pred])

print map(len, (thetas, pred, ytest_rho))
print 'r-squared', rsquare(pred, ytest_rho), 'rmse', rmse(pred, ytest_rho), 'spearman', spearman(pred, ytest_rho)

q = {}
idx = 1
for i,j,k, p in zip(thetas, pred, ytest_rho, postest):
    if i not in q: q[i] = []
    q[i].append((j,k))
    num_sites = len([iiii for iiii in p if iiii >= 0])
    idx+=1

for i in sorted(q): print i, len(q[i])


# idx = 1
# for i in sorted(q):
#     plt.subplot(2, 2, idx)
#     pred, real = [n[0] for n in q[i]], [n[1] for n in q[i]]
#     plt.scatter(real, pred, alpha=.3)
#     r = rsquare(np.array(real), np.array(pred))
#     plt.title('theta = '+str(i)+' r2:'+str(round(r,2)))
#     idx+=1
# 
# plt.show()
