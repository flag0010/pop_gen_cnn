from common import *
from os import popen
import numpy as np
import cPickle

rmse = lambda x,y: np.mean([(float(iii)-float(jjj))**2. for iii,jjj in zip(x,y)])**0.5
r2   = lambda x,y: np.corrcoef(map(float, x),map(float,y))[0][1]**2

def find_key(k):
    for i in nn_out:
        if abs(i-k) < .00000001: return i
    

res = get_file('test.data.results.csv', ',')
nn_out = {}
for index, theta, num_sites, predict_rho, real_rho in res[1:]:
    if real_rho in nn_out: print 'collision', real_rho #doesn't print
    nn_out[float(real_rho)] = map(float, [real_rho, predict_rho, theta])

theta2lkfile = {12.0006: "ldhat.data/new_lk_t12.txt", 6.0003: "ldhat.data/new_lk_t6.txt", 
              18.0009: "ldhat.data/new_lk_t18.txt",	60.003: "ldhat.data/new_lk_t60.txt",
              24.0012: "ldhat.data/new_lk_t24.txt"}


print len(nn_out)

p = get_file('ldhat.data/pairs.txt', ',')
myidx = len(p)
X = {}
for i in p:
    #print X
    z = i[0].split('_')
    theta = float(z[0].replace('sites/', ''))
    rho = float(z[1])
    #print z, theta, rho
    lkfile = theta2lkfile[theta]
    if theta not in X: X[theta] = []
    b = open('infile', 'w')
    b.write('ldhat.data/'+i[0]+'\n')
    b.write('ldhat.data/'+i[1]+'\n')
    b.write('1\n'+lkfile+'\n')
    b.write('0\n0\n0\n0\n0\n0\n')
    b.close()
    k = popen('LDhat-master/pairwise < infile')
    for j in k:
        if "Maximum at 4Ner(region) =" in j:
            pred_rho = j.strip().split()[4]
            mykey = find_key(rho)
            print myidx, i, rho, pred_rho, nn_out[mykey][1]
            X[theta].append([rho, float(pred_rho), nn_out[mykey][1]])
    myidx-=1
    if not myidx % 10:
        real = [u[0] for u in X[theta]]
        ldhat = [u[1] for u in X[theta]]
        nn = [u[2] for u in X[theta]]
        print rmse(real, ldhat), rmse(real, nn), rmse(ldhat, nn)
        print r2(real, ldhat), r2(real, nn), r2(ldhat, nn)
cPickle.dump(X, open('realrho.ldhatrho.nnrho.pickle', 'w'))

# for i in X:
#     z = X[i]
#     real_rho = [u[0] for u in z]
#     ldhat_rho = [u[1] for u in z]
#     keyd_rho = [find_key(i) for i in real_rho]
#     nn_rho = [nn_out[r][1] for r in keyd_rho]
#     print 'theta=', i
#     print rmse(real_rho, ldhat_rho)
#     print rmse(real_rho, nn_rho)
#     print rmse(nn_rho, ldhat_rho)
#     print
#     print r2(real_rho, ldhat_rho)
#     print r2(real_rho, nn_rho)
#     print r2(nn_rho, ldhat_rho)
#     print 
        #print i 
# ldhat.data/sites/24.0012_48.8999304542_fas
# ldhat.data/locs/24.0012_48.8999304542_loc
# 1
# LDhat-master/lk_files/lk_n50_t0.001
# 0
# 0
# 0
# 0
# 0
# 0


