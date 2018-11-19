def make_minor_allele_0(m):
    '''converts in place. assumes transposed snp matrix (indv on columns)
    In case of ties, defaults to doing nothing, which may not be ideal...
    possibly better to randomly decide to conv. or not for ties to destroy info.
    but only matters with even # of indv.'''
    q = m.shape
    cutoff = q[1] * 0.5
    i = np.where(m.sum(axis=1) < cutoff)
    m[i]*=-1
    m[i]+=1
