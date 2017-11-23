from sklearn.neighbors import NearestNeighbors
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
    return amat*2-1 
