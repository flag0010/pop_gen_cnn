from common import *
from itertools import izip
def tajD(N,S,k):
    '''N is the number of indv
       S is the number of seg. sites
       k is the (ave number of pairwise nuc. diffs), e.g. total number of pairwise diffs divided by choose(N,2)
       so:
       i1   AT|T|GGCG|A|CAG|T
       i2   AT|G|GGCG|C|CAG|A
       i3   AT|G|GGCG|A|CAG|A
       i4   AT|G|GGCG|C|CAG|T
       N = 4, S = 3, k = 11/6 = 1.83333
       D = 1.08976
       '''
    a1, a2 = 0.0, 0.0
    for b in range(1,N):
        a1 += 1. /  b
        a2 += 1. / (b * b)
    b1 = (N + 1.) / (3. * (N - 1.))
    b2 = (2. * N * N + 2. * N + 6.) / (9. * N * (N - 1.))
    c1 = b1 - 1. / a1
    c2 = b2 - (N + 2.) / (a1 * N) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    num = k - S / a1
    den = (e1 * S + e2 * S * (S - 1))**0.5
    D = num / den
    return D

def calc_S_and_k_from_seqs(list_of_seqs):
    ###first a function of taking the "sequence configuration" and converting it to a tally of pairwise diffs
    def pairwise_diffs2(x):
        if len(x) == 1:
            return 0
        else:
            q = 0.0
            for i,j in pairwise(x):
                q += i*j
            return q 
    #now loop through sites and tally S and k
    S, p = 0.0, 0.0
    for n in izip(*list_of_seqs):
        config = count_all(n).values()
        diffs = pairwise_diffs2(config)
        if diffs:  ## all the non-seg sites return zero, which is False
            p += diffs  
            S += 1
    k =  p / choose(len(list_of_seqs), 2)
    return S, k

def seq_boot(seqs):
    n = zip(*seqs)
    idxs = range(len(n))
    while 1:
        k = defaultdict(list)
        idx_list = sampler(idxs, len(n), replacement=True)
        for idx in idx_list:
            for pos, val in enumerate(n[idx]):
                k[pos].append(val)
        yield map(lambda s: ''.join(s), k.values())

def permute(seqs, reps = 1000):
    n = []
    repnum=0
    boot = seq_boot(seqs)
    while repnum < reps:
        rep = boot.next()
        N = len(rep)
        S, k = calc_S_and_k_from_seqs(rep)
        try:
            val = tajD(N,S,k)
        except:
            val = 'NA'
        if val != 'NA':
            n.append(val)
            repnum+=1
        #print repnum
    n.sort()
    #print len(n)
    lower_idx, upper_idx = int(len(n) * 0.025)-1, len(n) - int(len(n) * 0.025)
    mean = lambda s: sum(s) / float(len(s))
    return mean(n), n[lower_idx], n[upper_idx], len(n)
    
            
if __name__ == '__main__':
    seqs = ['ATTGGCGACAGT', 'ATGGGCGCCAGA', 'ATGGGCGACAGA', 'ATGGGCGCCAGT']
    ##works with non DNA seq data too
    #seqs = ['00001001000101'*100,
    #        '01101001000101'*100,
    #        '00001000010111'*100,
    #        '01101001010101'*100]
    N = len(seqs)
    S, k = calc_S_and_k_from_seqs(seqs)
    print N, S, k
    real = tajD(N,S,k)
    mean_perm, lower_ci, upper_ci, perms = permute(seqs)
    print "real Tajima's D", real
    print 'mean of bootstraps', mean_perm
    print 'lower 0.025 CI', lower_ci
    print 'upper 0.975 CI', upper_ci
    print 'num. replicates', perms 
    if lower_ci <= 0 <= upper_ci: print 'not sig. different from zero'
    else: print 'sig. different from zero'
    #from matplotlib import pyplot as plt
    #plt.hist(n)
    #plt.show()
    #idx = 0
    #for i in seqs:
    #    print '>'+str(idx)
    #    print i
    #    idx+=1
