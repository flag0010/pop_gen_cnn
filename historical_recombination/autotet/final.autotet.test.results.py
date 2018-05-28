import keras
from keras.preprocessing import sequence
import numpy as np
from matplotlib import pyplot as plt
import sys
import json
from random import shuffle
from gzip import GzipFile as gz
from sklearn.neighbors import NearestNeighbors


#model magic below. just a json string that specifies model structure
json_str ='{"backend": "tensorflow", "config": [{"class_name": "Merge", "config": {"name": "merge_2", "mode": "concat", "dot_axes": -1, "layers": [{"class_name": "Sequential", "config": [{"config": {"filters": 512, "strides": [1], "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_regularizer": {"config": {"l2": 9.999999747378752e-05, "l1": 0.0}, "class_name": "L1L2"}, "bias_constraint": null, "kernel_initializer": {"config": {"seed": null, "distribution": "uniform", "mode": "fan_avg", "scale": 1.0}, "class_name": "VarianceScaling"}, "kernel_size": [2], "kernel_constraint": null, "trainable": true, "activity_regularizer": null, "batch_input_shape": [null, 460, 12], "use_bias": true, "activation": "relu", "dilation_rate": [1], "bias_regularizer": null, "dtype": "float32", "padding": "valid", "name": "conv1d_4"}, "class_name": "Conv1D"}, {"config": {"activity_regularizer": null, "strides": [1], "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_regularizer": {"config": {"l2": 9.999999747378752e-05, "l1": 0.0}, "class_name": "L1L2"}, "bias_constraint": null, "kernel_initializer": {"config": {"seed": null, "distribution": "uniform", "mode": "fan_avg", "scale": 1.0}, "class_name": "VarianceScaling"}, "kernel_constraint": null, "kernel_size": [2], "trainable": true, "filters": 256, "use_bias": true, "activation": "relu", "dilation_rate": [1], "bias_regularizer": null, "padding": "valid", "name": "conv1d_5"}, "class_name": "Conv1D"}, {"config": {"padding": "valid", "strides": [2], "trainable": true, "pool_size": [2], "name": "average_pooling1d_3"}, "class_name": "AveragePooling1D"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_4"}, "class_name": "Dropout"}, {"config": {"activity_regularizer": null, "strides": [1], "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_regularizer": {"config": {"l2": 9.999999747378752e-05, "l1": 0.0}, "class_name": "L1L2"}, "bias_constraint": null, "kernel_initializer": {"config": {"seed": null, "distribution": "uniform", "mode": "fan_avg", "scale": 1.0}, "class_name": "VarianceScaling"}, "kernel_constraint": null, "kernel_size": [2], "trainable": true, "filters": 256, "use_bias": true, "activation": "relu", "dilation_rate": [1], "bias_regularizer": null, "padding": "valid", "name": "conv1d_6"}, "class_name": "Conv1D"}, {"config": {"padding": "valid", "strides": [2], "trainable": true, "pool_size": [2], "name": "average_pooling1d_4"}, "class_name": "AveragePooling1D"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_5"}, "class_name": "Dropout"}, {"config": {"trainable": true, "name": "flatten_2"}, "class_name": "Flatten"}]}, {"class_name": "Sequential", "config": [{"config": {"kernel_initializer": {"config": {"seed": null, "distribution": "uniform", "mode": "fan_avg", "scale": 1.0}, "class_name": "VarianceScaling"}, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_regularizer": {"config": {"l2": 9.999999747378752e-05, "l1": 0.0}, "class_name": "L1L2"}, "bias_constraint": null, "trainable": true, "kernel_constraint": null, "activity_regularizer": null, "units": 64, "name": "dense_4", "batch_input_shape": [null, 460], "use_bias": true, "activation": "relu", "bias_regularizer": null, "dtype": "float32"}, "class_name": "Dense"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_6"}, "class_name": "Dropout"}]}], "mode_type": "raw", "output_shape": null, "output_mask": null, "output_mask_type": "raw", "arguments": {}, "concat_axis": -1, "output_shape_type": "raw"}}, {"class_name": "Dense", "config": {"kernel_initializer": {"class_name": "VarianceScaling", "config": {"mode": "fan_avg", "scale": 1.0, "seed": null, "distribution": "uniform"}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L1L2", "config": {"l2": 9.999999747378752e-05, "l1": 0.0}}, "bias_constraint": null, "kernel_constraint": null, "activity_regularizer": null, "units": 256, "trainable": true, "use_bias": true, "activation": "relu", "bias_regularizer": null, "name": "dense_5"}}, {"class_name": "Dense", "config": {"kernel_initializer": {"class_name": "RandomNormal", "config": {"mean": 0.0, "stddev": 0.05, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_constraint": null, "kernel_constraint": null, "activity_regularizer": null, "units": 1, "trainable": true, "use_bias": true, "activation": "linear", "bias_regularizer": null, "name": "dense_6"}}], "class_name": "Sequential", "keras_version": "2.0.6"}'

#'{"keras_version": "2.0.6", "config": [{"config": {"mode": "concat", "arguments": {}, "layers": [{"config": [{"config": {"use_bias": true, "kernel_initializer": {"config": {"scale": 1.0, "mode": "fan_avg", "distribution": "uniform", "seed": null}, "class_name": "VarianceScaling"}, "batch_input_shape": [null, 460, 12], "strides": [1], "dilation_rate": [1], "kernel_constraint": null, "padding": "valid", "filters": 512, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_size": [2], "bias_regularizer": null, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "dtype": "float32", "kernel_regularizer": null, "activation": "relu", "name": "conv1d_7"}, "class_name": "Conv1D"}, {"config": {"use_bias": true, "kernel_initializer": {"config": {"scale": 1.0, "mode": "fan_avg", "distribution": "uniform", "seed": null}, "class_name": "VarianceScaling"}, "kernel_regularizer": null, "strides": [1], "dilation_rate": [1], "kernel_constraint": null, "padding": "valid", "filters": 256, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_size": [2], "bias_regularizer": null, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "activation": "relu", "name": "conv1d_8"}, "class_name": "Conv1D"}, {"config": {"padding": "valid", "pool_size": [2], "trainable": true, "strides": [2], "name": "average_pooling1d_5"}, "class_name": "AveragePooling1D"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_7"}, "class_name": "Dropout"}, {"config": {"use_bias": true, "kernel_initializer": {"config": {"scale": 1.0, "mode": "fan_avg", "distribution": "uniform", "seed": null}, "class_name": "VarianceScaling"}, "kernel_regularizer": null, "strides": [1], "dilation_rate": [1], "kernel_constraint": null, "padding": "valid", "filters": 256, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "kernel_size": [2], "bias_regularizer": null, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "activation": "relu", "name": "conv1d_9"}, "class_name": "Conv1D"}, {"config": {"padding": "valid", "pool_size": [2], "trainable": true, "strides": [2], "name": "average_pooling1d_6"}, "class_name": "AveragePooling1D"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_8"}, "class_name": "Dropout"}, {"config": {"trainable": true, "name": "flatten_3"}, "class_name": "Flatten"}], "class_name": "Sequential"}, {"config": [{"config": {"bias_regularizer": null, "kernel_initializer": {"config": {"scale": 1.0, "mode": "fan_avg", "distribution": "uniform", "seed": null}, "class_name": "VarianceScaling"}, "batch_input_shape": [null, 460], "kernel_constraint": null, "units": 64, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "use_bias": true, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "dtype": "float32", "kernel_regularizer": null, "activation": "relu", "name": "dense_7"}, "class_name": "Dense"}, {"config": {"rate": 0.25, "trainable": true, "name": "dropout_9"}, "class_name": "Dropout"}], "class_name": "Sequential"}], "output_shape": null, "output_mask_type": "raw", "output_shape_type": "raw", "concat_axis": -1, "name": "merge_3", "dot_axes": -1, "output_mask": null, "mode_type": "raw"}, "class_name": "Merge"}, {"config": {"use_bias": true, "kernel_initializer": {"config": {"scale": 1.0, "mode": "fan_avg", "distribution": "uniform", "seed": null}, "class_name": "VarianceScaling"}, "kernel_regularizer": null, "kernel_constraint": null, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "units": 256, "bias_regularizer": null, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "activation": "relu", "name": "dense_8"}, "class_name": "Dense"}, {"config": {"use_bias": true, "kernel_initializer": {"config": {"mean": 0.0, "stddev": 0.05, "seed": null}, "class_name": "RandomNormal"}, "kernel_regularizer": null, "kernel_constraint": null, "bias_initializer": {"config": {}, "class_name": "Zeros"}, "units": 1, "bias_regularizer": null, "activity_regularizer": null, "bias_constraint": null, "trainable": true, "activation": "linear", "name": "dense_9"}, "class_name": "Dense"}], "class_name": "Sequential", "backend": "tensorflow"}'

def transform_hap_to_autotet(x, coverage=25):
    idx = [[0+i,1+i,2+i,3+i] for i in range(0,48,4)]
    out = []
    for i in idx:
        cov = np.random.poisson(coverage, x.shape[1]) #draw seq coverage
        p = x[i].mean(axis=0)  #freq of 1 in 4 adjacent indv (i.e. a simulated autotetraploid genotype)
        p_count = np.random.binomial(cov, p) #reads in support of p allele, conditional on freq and cov.
        out.append(p_count * (1./cov)) #sim. freq of p_reads with odd way to do division while coercing to denom. to float
        #print out[-1]
    return np.array(out) * 2 - 1 #between -1 and 1

def sort_min_diff(amat):
    '''this function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted by genetic similarity.
    this problem is NP-hard, so here we use a nearest neighbors approx.  it's not perfect, but it's fast and generally performs ok.
    assumes your input matrix is a numpy array'''
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def convert_01_to_neg1_1(amat):
    '''convert standard binary 0/1 ms SNP matrix to -1/1 SNP matrix. B/c weights & biases are drawn from a distribution with mean=0
    choosing -1/1 (which is also mean=0) tends to help in training. assumes your input matrix is a numpy array'''
    return (amat*-2+1)*-1 

def rsquare(x,y):
    return np.corrcoef(x,y)[0][1]**2  #r-squared

def rmse(x,y):
    return np.sqrt(np.mean((x-y)**2))

a = (i.strip().split('asdfasdfd') for i in gz('autotet.test.data.LD.sims.txt.gz'))

k = range(48)
idx = 0
xvals, yvals = {},{}
seg_sites = {}
pos = {}
n = []
ctr = 0
true_idx = 0
maxL = 0

for i in a:
    if 'segsites' in i[0]:
        seg_sites[true_idx] = int(i[0].split()[-1])
        if seg_sites[true_idx] > maxL: maxL = seg_sites[true_idx]
    if 'positions' in i[0]:
        pos[true_idx] = [float(jj) for jj in i[0].split()[1:]]
    if './ms' in i[0]:
        true_idx+=1
        idx = 0
        theta_line = i[0].split()
        ms_vals = map(float, [theta_line[4], theta_line[6]])
        yvals[true_idx] = np.array(ms_vals)
    if 5 < idx < 54:
        n.append(np.array(map(int, i[0]), dtype='int8'))
    if len(n) == 48:
        n = np.array(n)
#        print 'init', n.shape, n.max(), n.min()
        shuffle(k)
        n = n[k]
#        print 'mix',  n.shape, n.max(), n.min()
        #print n.shape
        #print n
        n = transform_hap_to_autotet(n)
        #print n
#        print 'condense to tet', n.shape, n.max(), n.min()
        n = sort_min_diff(n)
        #print n
#        print 'sort', n.shape, n.max(), n.min(), '\n'
        #n = convert_01_to_neg1_1(n)
        #print 'rescale', n.shape, n.max(), n.min()
        #print n.shape
        #print n
        xvals[true_idx] =  n.T
        #print n.mean()
        idx = 0
        n = []
        if not true_idx % 100:
            print true_idx, len(xvals), len(yvals), min(xvals), max(xvals), min(yvals), max(yvals) #xvals[max(xvals)].shape, seg_sites[max(seg_sites)]
            #print xvals
        #if len(xvals) != len(yvals):
        #    print true_idx, idx, len(xvals), len(yvals), xvals[-2].shape, yvals[-2], '*********'
    #print idx
    idx+=1
    #print xvals, yvals, seg_sites, pos

new_yvals_test = []
new_xvals_test = []
new_pos_test = []

isect = set(xvals.keys()).intersection(yvals.keys()).intersection(pos.keys())
#print 'intersect len', len(isect), min(isect), max(isect)

for i in sorted(seg_sites):
    if i in xvals and i in yvals and i in pos:
        s = seg_sites[i]
        if s > 0 and yvals[i][0] > 3:
            new_pos_test.append(np.array(pos[i]))
            new_xvals_test.append(xvals[i])
            new_yvals_test.append(yvals[i])
#print new_xvals_test
xtest = sequence.pad_sequences(new_xvals_test, maxlen=460, padding='post', dtype='float32')#, value=-1., dtype='float32')
postest = sequence.pad_sequences(new_pos_test, maxlen=460, padding='post', value=-1., dtype='float32')
ytest_rho = np.array([i[1] for i in new_yvals_test])

xtest = (xtest+1)/2.

thetas=[i[0] for i in new_yvals_test]

np.savez_compressed('autotet.test.npz', xtest=xtest)
#print xtest
print "range", min(ytest_rho)/20000., max(ytest_rho)/20000
#for i in xtest:
#    print set([j for j in i.flatten()])
mean_test = 4.78838995038  #mean of training data ln(rho)

mod = keras.models.model_from_json(json_str)
#mod.load_weights('3rd.autotet.mergnet.weights')
mod.load_weights('auto.tet.regularized.merge.mod.weights')

pred = mod.predict([xtest,postest])
#print pred[:10]
#plt.hist(pred, bins=24)
#plt.show()

pred = np.exp([i[0]+mean_test for i in pred])


plt.scatter([ii/20000. for ii in ytest_rho], [pp/20000. for pp in pred], marker='.', color='k', alpha=0.25)
m = max(max([ii/20000. for ii in ytest_rho]), max([pp/20000. for pp in pred]))
plt.plot([0,m], [0,m], color='r')
plt.show()
#print map(len, (thetas, pred, ytest_rho))
print 'rsquare', rsquare(pred, ytest_rho), 'rmse', rmse(pred, ytest_rho)
q = {}
outfile = open('autotet.test.data.results.csv', 'w')
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

