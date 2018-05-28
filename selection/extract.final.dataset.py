from common import *
from gzip import GzipFile as gz
import numpy as np
from random import shuffle
from sklearn.neighbors import NearestNeighbors
from pprint import pprint
import multiprocessing

NCPU = 8


def sort_min_diff(amat):
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]

def load_data(xfile, maxlen = 5000, hard_stop = 1e60):
    g = [ii.strip().split() for ii in gz(xfile)] 
    k = [idx for idx,i in enumerate(g) if len(i) > 0 and  i[0].startswith('//')]
    f = []
    lens = []
    all_pos = []
    for idx, i in enumerate(k):
        l = g[i+3:i+211]
        pos = np.array(map(float, g[i+2][1:]), dtype="float32")
        #print l[0][:10], l[-1][:10]
        #print len(l)
        q = []
        for i in l:
            i = list(i[0])
            q.append(np.array([int(j) for j in i], dtype='int8'))
        #print len(q)
        #print q[0][:10], q[-1][:10]
        q = sort_min_diff(np.array(q)).T
        if q.shape[0] <= maxlen:
            f.append(q)
        #    print q.shape
            all_pos.append(pos)
        if len(f) > 9:
            if not len(f) % 10: print idx, len(f), len(all_pos)
        lens.append(len(q[0]))
        if len(f) >= hard_stop: break
    print len([i for i in lens if i > maxlen]), len(lens)
    print '*********', len(f), len(all_pos)
    return f, all_pos


#neut = 0, hard=1, hard near=2, soft=3, soft near = 4

rawxfiles = {'psmcHard_0.msOut.gz':2, 
          'psmcHard_1.msOut.gz':2,
          'psmcHard_10.msOut.gz':2,
          'psmcHard_2.msOut.gz':2,
          'psmcHard_3.msOut.gz':2,
          'psmcHard_4.msOut.gz':2,
          'psmcHard_5.msOut.gz':1,
          'psmcHard_6.msOut.gz':2,
          'psmcHard_7.msOut.gz':2,
          'psmcHard_8.msOut.gz':2,
          'psmcHard_9.msOut.gz':2,
          'psmcNeut.msOut.gz':0,
          'psmcSoft_0.msOut.gz':4,
          'psmcSoft_1.msOut.gz':4,
          'psmcSoft_10.msOut.gz':4,
          'psmcSoft_2.msOut.gz':4,
          'psmcSoft_3.msOut.gz':4,
          'psmcSoft_4.msOut.gz':4,
          'psmcSoft_5.msOut.gz':3,
          'psmcSoft_6.msOut.gz':4,
          'psmcSoft_7.msOut.gz':4,
          'psmcSoft_8.msOut.gz':4,
          'psmcSoft_9.msOut.gz':4}


testpath = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/testData/JPT/'
trainpath = 'san/data/dan/simulations/discoal_multipopStuff/spatialSVMSims/humanScanSims/trainingData/JPT/'

xfiles = {}
for i in rawxfiles:
    xfiles[testpath+i] = rawxfiles[i]
    xfiles[trainpath+i] = rawxfiles[i]
print xfiles


extra_data = '''psmcHard_0.batch_0.msOut.gz   psmcHard_7.batch_2.msOut.gz   psmcSoft_4.batch_0.msOut.gz   psmcSoft_5.batch_33.msOut.gz
psmcHard_0.batch_1.msOut.gz   psmcHard_7.batch_3.msOut.gz   psmcSoft_4.batch_1.msOut.gz   psmcSoft_5.batch_34.msOut.gz
psmcHard_0.batch_2.msOut.gz   psmcHard_8.batch_0.msOut.gz   psmcSoft_4.batch_2.msOut.gz   psmcSoft_5.batch_35.msOut.gz
psmcHard_0.batch_3.msOut.gz   psmcHard_8.batch_1.msOut.gz   psmcSoft_4.batch_3.msOut.gz   psmcSoft_5.batch_36.msOut.gz
psmcHard_10.batch_0.msOut.gz  psmcHard_8.batch_2.msOut.gz   psmcSoft_5.batch_0.msOut.gz   psmcSoft_5.batch_37.msOut.gz
psmcHard_10.batch_1.msOut.gz  psmcHard_8.batch_3.msOut.gz   psmcSoft_5.batch_10.msOut.gz  psmcSoft_5.batch_38.msOut.gz
psmcHard_10.batch_2.msOut.gz  psmcHard_9.batch_0.msOut.gz   psmcSoft_5.batch_11.msOut.gz  psmcSoft_5.batch_39.msOut.gz
psmcHard_10.batch_3.msOut.gz  psmcHard_9.batch_1.msOut.gz   psmcSoft_5.batch_12.msOut.gz  psmcSoft_5.batch_3.msOut.gz
psmcHard_1.batch_0.msOut.gz   psmcHard_9.batch_2.msOut.gz   psmcSoft_5.batch_13.msOut.gz  psmcSoft_5.batch_4.msOut.gz
psmcHard_1.batch_1.msOut.gz   psmcHard_9.batch_3.msOut.gz   psmcSoft_5.batch_14.msOut.gz  psmcSoft_5.batch_5.msOut.gz
psmcHard_1.batch_2.msOut.gz   psmcSoft_0.batch_0.msOut.gz   psmcSoft_5.batch_15.msOut.gz  psmcSoft_5.batch_6.msOut.gz
psmcHard_1.batch_3.msOut.gz   psmcSoft_0.batch_1.msOut.gz   psmcSoft_5.batch_16.msOut.gz  psmcSoft_5.batch_7.msOut.gz
psmcHard_2.batch_0.msOut.gz   psmcSoft_0.batch_2.msOut.gz   psmcSoft_5.batch_17.msOut.gz  psmcSoft_5.batch_8.msOut.gz
psmcHard_2.batch_1.msOut.gz   psmcSoft_0.batch_3.msOut.gz   psmcSoft_5.batch_18.msOut.gz  psmcSoft_5.batch_9.msOut.gz
psmcHard_2.batch_2.msOut.gz   psmcSoft_10.batch_0.msOut.gz  psmcSoft_5.batch_19.msOut.gz  psmcSoft_6.batch_0.msOut.gz
psmcHard_2.batch_3.msOut.gz   psmcSoft_10.batch_1.msOut.gz  psmcSoft_5.batch_1.msOut.gz   psmcSoft_6.batch_1.msOut.gz
psmcHard_3.batch_0.msOut.gz   psmcSoft_10.batch_2.msOut.gz  psmcSoft_5.batch_20.msOut.gz  psmcSoft_6.batch_2.msOut.gz
psmcHard_3.batch_1.msOut.gz   psmcSoft_10.batch_3.msOut.gz  psmcSoft_5.batch_21.msOut.gz  psmcSoft_6.batch_3.msOut.gz
psmcHard_3.batch_2.msOut.gz   psmcSoft_1.batch_0.msOut.gz   psmcSoft_5.batch_22.msOut.gz  psmcSoft_7.batch_0.msOut.gz
psmcHard_3.batch_3.msOut.gz   psmcSoft_1.batch_1.msOut.gz   psmcSoft_5.batch_23.msOut.gz  psmcSoft_7.batch_1.msOut.gz
psmcHard_4.batch_0.msOut.gz   psmcSoft_1.batch_2.msOut.gz   psmcSoft_5.batch_24.msOut.gz  psmcSoft_7.batch_2.msOut.gz
psmcHard_4.batch_1.msOut.gz   psmcSoft_1.batch_3.msOut.gz   psmcSoft_5.batch_25.msOut.gz  psmcSoft_7.batch_3.msOut.gz
psmcHard_4.batch_2.msOut.gz   psmcSoft_2.batch_0.msOut.gz   psmcSoft_5.batch_26.msOut.gz  psmcSoft_8.batch_0.msOut.gz
psmcHard_4.batch_3.msOut.gz   psmcSoft_2.batch_1.msOut.gz   psmcSoft_5.batch_27.msOut.gz  psmcSoft_8.batch_1.msOut.gz
psmcHard_6.batch_0.msOut.gz   psmcSoft_2.batch_2.msOut.gz   psmcSoft_5.batch_28.msOut.gz  psmcSoft_8.batch_2.msOut.gz
psmcHard_6.batch_1.msOut.gz   psmcSoft_2.batch_3.msOut.gz   psmcSoft_5.batch_29.msOut.gz  psmcSoft_8.batch_3.msOut.gz
psmcHard_6.batch_2.msOut.gz   psmcSoft_3.batch_0.msOut.gz   psmcSoft_5.batch_2.msOut.gz   psmcSoft_9.batch_0.msOut.gz
psmcHard_6.batch_3.msOut.gz   psmcSoft_3.batch_1.msOut.gz   psmcSoft_5.batch_30.msOut.gz  psmcSoft_9.batch_1.msOut.gz
psmcHard_7.batch_0.msOut.gz   psmcSoft_3.batch_2.msOut.gz   psmcSoft_5.batch_31.msOut.gz  psmcSoft_9.batch_2.msOut.gz
psmcHard_7.batch_1.msOut.gz   psmcSoft_3.batch_3.msOut.gz   psmcSoft_5.batch_32.msOut.gz  psmcSoft_9.batch_3.msOut.gz'''.split()


for i in extra_data:
    if "Hard" in i: xfiles[i] = 2
    else:
        if 'Soft_5' in i: xfiles[i] = 3
        else: xfiles[i] = 4


more_hard = '''psmcHard_5.batch_0.msOut.gz   psmcHard_5.batch_17.msOut.gz  psmcHard_5.batch_24.msOut.gz  psmcHard_5.batch_31.msOut.gz  psmcHard_5.batch_39.msOut.gz
psmcHard_5.batch_10.msOut.gz  psmcHard_5.batch_18.msOut.gz  psmcHard_5.batch_25.msOut.gz  psmcHard_5.batch_32.msOut.gz  psmcHard_5.batch_3.msOut.gz
psmcHard_5.batch_11.msOut.gz  psmcHard_5.batch_19.msOut.gz  psmcHard_5.batch_26.msOut.gz  psmcHard_5.batch_33.msOut.gz  psmcHard_5.batch_4.msOut.gz
psmcHard_5.batch_12.msOut.gz  psmcHard_5.batch_1.msOut.gz   psmcHard_5.batch_27.msOut.gz  psmcHard_5.batch_34.msOut.gz  psmcHard_5.batch_5.msOut.gz
psmcHard_5.batch_13.msOut.gz  psmcHard_5.batch_20.msOut.gz  psmcHard_5.batch_28.msOut.gz  psmcHard_5.batch_35.msOut.gz  psmcHard_5.batch_6.msOut.gz
psmcHard_5.batch_14.msOut.gz  psmcHard_5.batch_21.msOut.gz  psmcHard_5.batch_29.msOut.gz  psmcHard_5.batch_36.msOut.gz  psmcHard_5.batch_7.msOut.gz
psmcHard_5.batch_15.msOut.gz  psmcHard_5.batch_22.msOut.gz  psmcHard_5.batch_2.msOut.gz   psmcHard_5.batch_37.msOut.gz  psmcHard_5.batch_8.msOut.gz
psmcHard_5.batch_16.msOut.gz  psmcHard_5.batch_23.msOut.gz  psmcHard_5.batch_30.msOut.gz  psmcHard_5.batch_38.msOut.gz  psmcHard_5.batch_9.msOut.gz'''.split()

for i in more_hard:
    xfiles[i] = 1

more_neut = '''psmcNeut.batch_0.msOut.gz   psmcNeut.batch_17.msOut.gz  psmcNeut.batch_24.msOut.gz  psmcNeut.batch_31.msOut.gz  psmcNeut.batch_39.msOut.gz
psmcNeut.batch_10.msOut.gz  psmcNeut.batch_18.msOut.gz  psmcNeut.batch_25.msOut.gz  psmcNeut.batch_32.msOut.gz  psmcNeut.batch_3.msOut.gz
psmcNeut.batch_11.msOut.gz  psmcNeut.batch_19.msOut.gz  psmcNeut.batch_26.msOut.gz  psmcNeut.batch_33.msOut.gz  psmcNeut.batch_4.msOut.gz
psmcNeut.batch_12.msOut.gz  psmcNeut.batch_1.msOut.gz   psmcNeut.batch_27.msOut.gz  psmcNeut.batch_34.msOut.gz  psmcNeut.batch_5.msOut.gz
psmcNeut.batch_13.msOut.gz  psmcNeut.batch_20.msOut.gz  psmcNeut.batch_28.msOut.gz  psmcNeut.batch_35.msOut.gz  psmcNeut.batch_6.msOut.gz
psmcNeut.batch_14.msOut.gz  psmcNeut.batch_21.msOut.gz  psmcNeut.batch_29.msOut.gz  psmcNeut.batch_36.msOut.gz  psmcNeut.batch_7.msOut.gz
psmcNeut.batch_15.msOut.gz  psmcNeut.batch_22.msOut.gz  psmcNeut.batch_2.msOut.gz   psmcNeut.batch_37.msOut.gz  psmcNeut.batch_8.msOut.gz
psmcNeut.batch_16.msOut.gz  psmcNeut.batch_23.msOut.gz  psmcNeut.batch_30.msOut.gz  psmcNeut.batch_38.msOut.gz  psmcNeut.batch_9.msOut.gz'''.split()

for i in more_neut:
    xfiles[i] = 0

print '****************'
print xfiles

print len(xfiles)

def proc_one_file(i, xfiles=xfiles):
    print i
    code = xfiles[i]
    xdata,posdata = load_data(i)
    ydata = [code for i in range(len(xdata))]
    return xdata, posdata, ydata

if __name__ == '__main__':
    x,y,pos = [],[],[]
    P = multiprocessing.Pool(NCPU) #multicore initialized to NCPUs
    for xdata, posdata, ydata in P.imap(proc_one_file, xfiles.keys()):
        x.extend(xdata)
        pos.extend(posdata)
        y.extend(ydata)
    
    x = np.array(x)
    y = np.array(y)
    all_pos = np.array(pos)
 
    print len(x), len(y), len(all_pos)
 
    shf = range(len(x))
    shuffle(shf)

    x = x[shf]
    y = y[shf]
    all_pos = all_pos[shf]

    final_ct, final_pos = [], []
    idx = 0
    for i in y:
        if final_ct.count(i) < 2000:
            final_ct.append(i)
            final_pos.append(idx)
        idx+=1
        
    rest_pos = sorted(set(range(len(x))).difference(final_pos))
    
    xfinal, yfinal, posfinal = x[final_pos[:5000]], y[final_pos[:5000]], all_pos[final_pos[:5000]]
    xtest, ytest, postest = x[final_pos[5000:]], y[final_pos[5000:]], all_pos[final_pos[5000:]]
    
    #above enforces equal proportions of each of the 5 classes in the final test data, but not training/val
    
    x, y, all_pos = x[rest_pos], y[rest_pos], all_pos[rest_pos]
    
    postrain_0, xtrain_0, ytrain_0 = all_pos[0:3000], x[0:3000], y[0:3000]
    postrain_3000, xtrain_3000, ytrain_3000 = all_pos[3000:6000], x[3000:6000], y[3000:6000]
    postrain_6000, xtrain_6000, ytrain_6000 = all_pos[6000:9000], x[6000:9000], y[6000:9000]
    postrain_9000, xtrain_9000, ytrain_9000 = all_pos[9000:12000], x[9000:12000], y[9000:12000]
    postrain_12000, xtrain_12000, ytrain_12000 = all_pos[12000:15000], x[12000:15000], y[12000:15000]
    postrain_15000, xtrain_15000, ytrain_15000 = all_pos[15000:18000], x[15000:18000], y[15000:18000]
    postrain_18000, xtrain_18000, ytrain_18000 = all_pos[18000:21000], x[18000:21000], y[18000:21000]
    postrain_21000, xtrain_21000, ytrain_21000 = all_pos[21000:24000], x[21000:24000], y[21000:24000]
    postrain_24000, xtrain_24000, ytrain_24000 = all_pos[24000:27000], x[24000:27000], y[24000:27000]
    postrain_27000, xtrain_27000, ytrain_27000 = all_pos[27000:30000], x[27000:30000], y[27000:30000]
    postrain_30000, xtrain_30000, ytrain_30000 = all_pos[30000:33000], x[30000:33000], y[30000:33000]
    postrain_33000, xtrain_33000, ytrain_33000 = all_pos[33000:36000], x[33000:36000], y[33000:36000]
    postrain_36000, xtrain_36000, ytrain_36000 = all_pos[36000:39000], x[36000:39000], y[36000:39000]
    postrain_39000, xtrain_39000, ytrain_39000 = all_pos[39000:42000], x[39000:42000], y[39000:42000]
    postrain_42000, xtrain_42000, ytrain_42000 = all_pos[42000:45000], x[42000:45000], y[42000:45000]
    postrain_45000, xtrain_45000, ytrain_45000 = all_pos[45000:48000], x[45000:48000], y[45000:48000]
    postrain_48000, xtrain_48000, ytrain_48000 = all_pos[48000:51000], x[48000:51000], y[48000:51000]
    postrain_51000, xtrain_51000, ytrain_51000 = all_pos[51000:54000], x[51000:54000], y[51000:54000]
    postrain_54000, xtrain_54000, ytrain_54000 = all_pos[54000:57000], x[54000:57000], y[54000:57000]
    postrain_57000, xtrain_57000, ytrain_57000 = all_pos[57000:60000], x[57000:60000], y[57000:60000]
    postrain_60000, xtrain_60000, ytrain_60000 = all_pos[60000:63000], x[60000:63000], y[60000:63000]
    postrain_63000, xtrain_63000, ytrain_63000 = all_pos[63000:66000], x[63000:66000], y[63000:66000]
    postrain_66000, xtrain_66000, ytrain_66000 = all_pos[66000:69000], x[66000:69000], y[66000:69000]
    postrain_69000, xtrain_69000, ytrain_69000 = all_pos[69000:72000], x[69000:72000], y[69000:72000]
    postrain_72000, xtrain_72000, ytrain_72000 = all_pos[72000:75000], x[72000:75000], y[72000:75000]
    postrain_75000, xtrain_75000, ytrain_75000 = all_pos[75000:78000], x[75000:78000], y[75000:78000]
    postrain_78000, xtrain_78000, ytrain_78000 = all_pos[78000:81000], x[78000:81000], y[78000:81000]
    postrain_81000, xtrain_81000, ytrain_81000 = all_pos[81000:84000], x[81000:84000], y[81000:84000]
    postrain_84000, xtrain_84000, ytrain_84000 = all_pos[84000:87000], x[84000:87000], y[84000:87000]
    postrain_87000, xtrain_87000, ytrain_87000 = all_pos[87000:90000], x[87000:90000], y[87000:90000]
    postrain_90000, xtrain_90000, ytrain_90000 = all_pos[90000:93000], x[90000:93000], y[90000:93000]
    postrain_93000, xtrain_93000, ytrain_93000 = all_pos[93000:96000], x[93000:96000], y[93000:96000]
    postrain_96000, xtrain_96000, ytrain_96000 = all_pos[96000:99000], x[96000:99000], y[96000:99000]
    postrain_99000, xtrain_99000, ytrain_99000 = all_pos[99000:102000], x[99000:102000], y[99000:102000]
    postrain_102000, xtrain_102000, ytrain_102000 = all_pos[102000:105000], x[102000:105000], y[102000:105000]
    postrain_105000, xtrain_105000, ytrain_105000 = all_pos[105000:108000], x[105000:108000], y[105000:108000]
    postrain_108000, xtrain_108000, ytrain_108000 = all_pos[108000:111000], x[108000:111000], y[108000:111000]
    postrain_111000, xtrain_111000, ytrain_111000 = all_pos[111000:114000], x[111000:114000], y[111000:114000]
    postrain_114000, xtrain_114000, ytrain_114000 = all_pos[114000:117000], x[114000:117000], y[114000:117000]
    postrain_117000, xtrain_117000, ytrain_117000 = all_pos[117000:120000], x[117000:120000], y[117000:120000]
    postrain_120000, xtrain_120000, ytrain_120000 = all_pos[120000:123000], x[120000:123000], y[120000:123000]
    postrain_123000, xtrain_123000, ytrain_123000 = all_pos[123000:126000], x[123000:126000], y[123000:126000]
    postrain_126000, xtrain_126000, ytrain_126000 = all_pos[126000:129000], x[126000:129000], y[126000:129000]
    postrain_129000, xtrain_129000, ytrain_129000 = all_pos[129000:132000], x[129000:132000], y[129000:132000]
    postrain_132000, xtrain_132000, ytrain_132000 = all_pos[132000:135000], x[132000:135000], y[132000:135000]
    postrain_135000, xtrain_135000, ytrain_135000 = all_pos[135000:138000], x[135000:138000], y[135000:138000]
    postrain_138000, xtrain_138000, ytrain_138000 = all_pos[138000:141000], x[138000:141000], y[138000:141000]
    postrain_141000, xtrain_141000, ytrain_141000 = all_pos[141000:144000], x[141000:144000], y[141000:144000]
    postrain_144000, xtrain_144000, ytrain_144000 = all_pos[144000:147000], x[144000:147000], y[144000:147000]
    postrain_147000, xtrain_147000, ytrain_147000 = all_pos[147000:150000], x[147000:150000], y[147000:150000]
    postrain_150000, xtrain_150000, ytrain_150000 = all_pos[150000:153000], x[150000:153000], y[150000:153000]
    postrain_153000, xtrain_153000, ytrain_153000 = all_pos[153000:156000], x[153000:156000], y[153000:156000]
    postrain_156000, xtrain_156000, ytrain_156000 = all_pos[156000:159000], x[156000:159000], y[156000:159000]
    postrain_159000, xtrain_159000, ytrain_159000 = all_pos[159000:162000], x[159000:162000], y[159000:162000]
    postrain_162000, xtrain_162000, ytrain_162000 = all_pos[162000:165000], x[162000:165000], y[162000:165000]
    postrain_165000, xtrain_165000, ytrain_165000 = all_pos[165000:168000], x[165000:168000], y[165000:168000]
    postrain_168000, xtrain_168000, ytrain_168000 = all_pos[168000:171000], x[168000:171000], y[168000:171000]
    postrain_171000, xtrain_171000, ytrain_171000 = all_pos[171000:174000], x[171000:174000], y[171000:174000]
    postrain_174000, xtrain_174000, ytrain_174000 = all_pos[174000:177000], x[174000:177000], y[174000:177000]
    postrain_177000, xtrain_177000, ytrain_177000 = all_pos[177000:180000], x[177000:180000], y[177000:180000]
    postrain_180000, xtrain_180000, ytrain_180000 = all_pos[180000:183000], x[180000:183000], y[180000:183000]
    postrain_183000, xtrain_183000, ytrain_183000 = all_pos[183000:186000], x[183000:186000], y[183000:186000]
    postrain_186000, xtrain_186000, ytrain_186000 = all_pos[186000:189000], x[186000:189000], y[186000:189000]
    postrain_189000, xtrain_189000, ytrain_189000 = all_pos[189000:192000], x[189000:192000], y[189000:192000]
    postrain_192000, xtrain_192000, ytrain_192000 = all_pos[192000:195000], x[192000:195000], y[192000:195000]
    postrain_195000, xtrain_195000, ytrain_195000 = all_pos[195000:198000], x[195000:198000], y[195000:198000]
    postrain_198000, xtrain_198000, ytrain_198000 = all_pos[198000:201000], x[198000:201000], y[198000:201000]
    postrain_201000, xtrain_201000, ytrain_201000 = all_pos[201000:204000], x[201000:204000], y[201000:204000]
    postrain_204000, xtrain_204000, ytrain_204000 = all_pos[204000:207000], x[204000:207000], y[204000:207000]
    postrain_207000, xtrain_207000, ytrain_207000 = all_pos[207000:210000], x[207000:210000], y[207000:210000]
    postrain_210000, xtrain_210000, ytrain_210000 = all_pos[210000:213000], x[210000:213000], y[210000:213000]
    postrain_213000, xtrain_213000, ytrain_213000 = all_pos[213000:216000], x[213000:216000], y[213000:216000]
    postrain_216000, xtrain_216000, ytrain_216000 = all_pos[216000:219000], x[216000:219000], y[216000:219000]
    postrain_219000, xtrain_219000, ytrain_219000 = all_pos[219000:222000], x[219000:222000], y[219000:222000]
    postrain_222000, xtrain_222000, ytrain_222000 = all_pos[222000:225000], x[222000:225000], y[222000:225000]
    postrain_225000, xtrain_225000, ytrain_225000 = all_pos[225000:228000], x[225000:228000], y[225000:228000]
    postrain_228000, xtrain_228000, ytrain_228000 = all_pos[228000:231000], x[228000:231000], y[228000:231000]
    postrain_231000, xtrain_231000, ytrain_231000 = all_pos[231000:], x[231000:], y[231000:]
    
    np.savez_compressed('final.split.up.seln.big.npz', xfinal=xfinal, yfinal=yfinal, posfinal=posfinal, xtest=xtest, ytest=ytest, postest=postest,
                                                      postrain_0=postrain_0, xtrain_0= xtrain_0, ytrain_0=ytrain_0,
                                                      postrain_3000=postrain_3000, xtrain_3000= xtrain_3000, ytrain_3000=ytrain_3000,
                                                      postrain_6000=postrain_6000, xtrain_6000= xtrain_6000, ytrain_6000=ytrain_6000,
                                                      postrain_9000=postrain_9000, xtrain_9000= xtrain_9000, ytrain_9000=ytrain_9000,
                                                      postrain_12000=postrain_12000, xtrain_12000= xtrain_12000, ytrain_12000=ytrain_12000,
                                                      postrain_15000=postrain_15000, xtrain_15000= xtrain_15000, ytrain_15000=ytrain_15000,
                                                      postrain_18000=postrain_18000, xtrain_18000= xtrain_18000, ytrain_18000=ytrain_18000,
                                                      postrain_21000=postrain_21000, xtrain_21000= xtrain_21000, ytrain_21000=ytrain_21000,
                                                      postrain_24000=postrain_24000, xtrain_24000= xtrain_24000, ytrain_24000=ytrain_24000,
                                                      postrain_27000=postrain_27000, xtrain_27000= xtrain_27000, ytrain_27000=ytrain_27000,
                                                      postrain_30000=postrain_30000, xtrain_30000= xtrain_30000, ytrain_30000=ytrain_30000,
                                                      postrain_33000=postrain_33000, xtrain_33000= xtrain_33000, ytrain_33000=ytrain_33000,
                                                      postrain_36000=postrain_36000, xtrain_36000= xtrain_36000, ytrain_36000=ytrain_36000,
                                                      postrain_39000=postrain_39000, xtrain_39000= xtrain_39000, ytrain_39000=ytrain_39000,
                                                      postrain_42000=postrain_42000, xtrain_42000= xtrain_42000, ytrain_42000=ytrain_42000,
                                                      postrain_45000=postrain_45000, xtrain_45000= xtrain_45000, ytrain_45000=ytrain_45000,
                                                      postrain_48000=postrain_48000, xtrain_48000= xtrain_48000, ytrain_48000=ytrain_48000,
                                                      postrain_51000=postrain_51000, xtrain_51000= xtrain_51000, ytrain_51000=ytrain_51000,
                                                      postrain_54000=postrain_54000, xtrain_54000= xtrain_54000, ytrain_54000=ytrain_54000,
                                                      postrain_57000=postrain_57000, xtrain_57000= xtrain_57000, ytrain_57000=ytrain_57000,
                                                      postrain_60000=postrain_60000, xtrain_60000= xtrain_60000, ytrain_60000=ytrain_60000,
                                                      postrain_63000=postrain_63000, xtrain_63000= xtrain_63000, ytrain_63000=ytrain_63000,
                                                      postrain_66000=postrain_66000, xtrain_66000= xtrain_66000, ytrain_66000=ytrain_66000,
                                                      postrain_69000=postrain_69000, xtrain_69000= xtrain_69000, ytrain_69000=ytrain_69000,
                                                      postrain_72000=postrain_72000, xtrain_72000= xtrain_72000, ytrain_72000=ytrain_72000,
                                                      postrain_75000=postrain_75000, xtrain_75000= xtrain_75000, ytrain_75000=ytrain_75000,
                                                      postrain_78000=postrain_78000, xtrain_78000= xtrain_78000, ytrain_78000=ytrain_78000,
                                                      postrain_81000=postrain_81000, xtrain_81000= xtrain_81000, ytrain_81000=ytrain_81000,
                                                      postrain_84000=postrain_84000, xtrain_84000= xtrain_84000, ytrain_84000=ytrain_84000,
                                                      postrain_87000=postrain_87000, xtrain_87000= xtrain_87000, ytrain_87000=ytrain_87000,
                                                      postrain_90000=postrain_90000, xtrain_90000= xtrain_90000, ytrain_90000=ytrain_90000,
                                                      postrain_93000=postrain_93000, xtrain_93000= xtrain_93000, ytrain_93000=ytrain_93000,
                                                      postrain_96000=postrain_96000, xtrain_96000= xtrain_96000, ytrain_96000=ytrain_96000,
                                                      postrain_99000=postrain_99000, xtrain_99000= xtrain_99000, ytrain_99000=ytrain_99000,
                                                      postrain_102000=postrain_102000, xtrain_102000= xtrain_102000, ytrain_102000=ytrain_102000,
                                                      postrain_105000=postrain_105000, xtrain_105000= xtrain_105000, ytrain_105000=ytrain_105000,
                                                      postrain_108000=postrain_108000, xtrain_108000= xtrain_108000, ytrain_108000=ytrain_108000,
                                                      postrain_111000=postrain_111000, xtrain_111000= xtrain_111000, ytrain_111000=ytrain_111000,
                                                      postrain_114000=postrain_114000, xtrain_114000= xtrain_114000, ytrain_114000=ytrain_114000,
                                                      postrain_117000=postrain_117000, xtrain_117000= xtrain_117000, ytrain_117000=ytrain_117000,
                                                      postrain_120000=postrain_120000, xtrain_120000= xtrain_120000, ytrain_120000=ytrain_120000,
                                                      postrain_123000=postrain_123000, xtrain_123000= xtrain_123000, ytrain_123000=ytrain_123000,
                                                      postrain_126000=postrain_126000, xtrain_126000= xtrain_126000, ytrain_126000=ytrain_126000,
                                                      postrain_129000=postrain_129000, xtrain_129000= xtrain_129000, ytrain_129000=ytrain_129000,
                                                      postrain_132000=postrain_132000, xtrain_132000= xtrain_132000, ytrain_132000=ytrain_132000,
                                                      postrain_135000=postrain_135000, xtrain_135000= xtrain_135000, ytrain_135000=ytrain_135000,
                                                      postrain_138000=postrain_138000, xtrain_138000= xtrain_138000, ytrain_138000=ytrain_138000,
                                                      postrain_141000=postrain_141000, xtrain_141000= xtrain_141000, ytrain_141000=ytrain_141000,
                                                      postrain_144000=postrain_144000, xtrain_144000= xtrain_144000, ytrain_144000=ytrain_144000,
                                                      postrain_147000=postrain_147000, xtrain_147000= xtrain_147000, ytrain_147000=ytrain_147000,
                                                      postrain_150000=postrain_150000, xtrain_150000= xtrain_150000, ytrain_150000=ytrain_150000,
                                                      postrain_153000=postrain_153000, xtrain_153000= xtrain_153000, ytrain_153000=ytrain_153000,
                                                      postrain_156000=postrain_156000, xtrain_156000= xtrain_156000, ytrain_156000=ytrain_156000,
                                                      postrain_159000=postrain_159000, xtrain_159000= xtrain_159000, ytrain_159000=ytrain_159000,
                                                      postrain_162000=postrain_162000, xtrain_162000= xtrain_162000, ytrain_162000=ytrain_162000,
                                                      postrain_165000=postrain_165000, xtrain_165000= xtrain_165000, ytrain_165000=ytrain_165000,
                                                      postrain_168000=postrain_168000, xtrain_168000= xtrain_168000, ytrain_168000=ytrain_168000,
                                                      postrain_171000=postrain_171000, xtrain_171000= xtrain_171000, ytrain_171000=ytrain_171000,
                                                      postrain_174000=postrain_174000, xtrain_174000= xtrain_174000, ytrain_174000=ytrain_174000,
                                                      postrain_177000=postrain_177000, xtrain_177000= xtrain_177000, ytrain_177000=ytrain_177000,
                                                      postrain_180000=postrain_180000, xtrain_180000= xtrain_180000, ytrain_180000=ytrain_180000,
                                                      postrain_183000=postrain_183000, xtrain_183000= xtrain_183000, ytrain_183000=ytrain_183000,
                                                      postrain_186000=postrain_186000, xtrain_186000= xtrain_186000, ytrain_186000=ytrain_186000,
                                                      postrain_189000=postrain_189000, xtrain_189000= xtrain_189000, ytrain_189000=ytrain_189000,
                                                      postrain_192000=postrain_192000, xtrain_192000= xtrain_192000, ytrain_192000=ytrain_192000,
                                                      postrain_195000=postrain_195000, xtrain_195000= xtrain_195000, ytrain_195000=ytrain_195000,
                                                      postrain_198000=postrain_198000, xtrain_198000= xtrain_198000, ytrain_198000=ytrain_198000,
                                                      postrain_201000=postrain_201000, xtrain_201000= xtrain_201000, ytrain_201000=ytrain_201000,
                                                      postrain_204000=postrain_204000, xtrain_204000= xtrain_204000, ytrain_204000=ytrain_204000,
                                                      postrain_207000=postrain_207000, xtrain_207000= xtrain_207000, ytrain_207000=ytrain_207000,
                                                      postrain_210000=postrain_210000, xtrain_210000= xtrain_210000, ytrain_210000=ytrain_210000,
                                                      postrain_213000=postrain_213000, xtrain_213000= xtrain_213000, ytrain_213000=ytrain_213000,
                                                      postrain_216000=postrain_216000, xtrain_216000= xtrain_216000, ytrain_216000=ytrain_216000,
                                                      postrain_219000=postrain_219000, xtrain_219000= xtrain_219000, ytrain_219000=ytrain_219000,
                                                      postrain_222000=postrain_222000, xtrain_222000= xtrain_222000, ytrain_222000=ytrain_222000,
                                                      postrain_225000=postrain_225000, xtrain_225000= xtrain_225000, ytrain_225000=ytrain_225000,
                                                      postrain_228000=postrain_228000, xtrain_228000= xtrain_228000, ytrain_228000=ytrain_228000,
                                                      postrain_231000=postrain_231000, xtrain_231000= xtrain_231000, ytrain_231000=ytrain_231000)
