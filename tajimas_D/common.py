from string import maketrans
import operator
import random
from collections import defaultdict

###Stats and math functions
def weighted_sampler(pop_dict):
    """randomly sample a dictionary's keys based on weights stored as values
       example:
           m = {'a':3, 'b':2, 'c':5}
           samps = [weighted_sampler(m) for _ in xrange(1000)]
           #samps should be a ~ 300, b ~ 200, and c ~ 500
           >>> samps.count('a')
           304
           >>> samps.count('b')
           211
           >>> samps.count('c')
           485
       of course, being a random sampler your results will vary"""
    ch = random.random() * sum(pop_dict.values())
    f = sorted(pop_dict.keys())
    for i, w in enumerate([pop_dict[x] for x in f]):
        ch -= w
        if ch < 0: return f[i]

def choose(n,k):
    '''implements binomial coefficient function
       see: https://en.wikipedia.org/wiki/Binomial_coefficient 
       performance not tested on really large values'''
    return reduce(lambda a,b: a*(n-b)/(b+1),xrange(k),1)

def sampler(pop, size, replacement=False):
    '''a quick re-implementation of the python random sampler that
       allows for sampling with or without replacement (pythons builtin only
       allows without replacement)'''
    if replacement:
        return [random.choice(pop) for i in xrange(size)]
    else:
        return random.sample(pop, size)          

def rank(x):
    '''returns the sample rank of the elements in a list'''
    out={}
    idx=0
    for i in x:
        out[idx] = i
        idx+=1
    p1 =  (j[0] for j in sorted(sort_dict_by_val(out), key=lambda s: s[1]))
    p2 = range(len(x))
    idx=0
    for i in p1:
        p2[i] = idx
        idx+=1
    return p2

def order(x):
    '''returns the sample indeces that would return the list in sorted order
       ie: 
       x = (4,3,406,5)
       sorted(x) == [x[i] for i in order(x)]'''
    out={}
    idx=0
    for i in x:
        out[idx] = i
        idx+=1
    p1 =  [j[0] for j in sorted(sort_dict_by_val(out), key=lambda s: s[1])]
    return p1

###Useful functions for bioinformatics
###NOTE: biopython offers more robust versions, but sometimes you just need something quick and dirty
def revcom (s):
    '''returns the reverse complement of a DNA sequence string
       only accepts ACGT, upper or lowercase'''
    rv_s = s[::-1] #strange python string reversal, it works!
    trans_table = maketrans("atcgATCG", "tagcTAGC")
    rv_comp_s = rv_s.translate(trans_table)
    return rv_comp_s

def get_fasta(file_name):
    '''read a properly formated fasta and return a dict
       with key=readname and value=sequence
       reads the whole file in'''
    d = [i.strip() for i in open(file_name,'r')]
    out={}
    for i in d:
        if i.startswith('>'):
            curr_seq = i[1:] 
            out[curr_seq] = []
        else:
            out[curr_seq].append(i)
    for i in out:
        out[i] = ''.join(out[i])
    return out

def get_fasta_buffer(file_name):
    '''An efficient fasta reader that is buffered and therefore
       useful for big fasta files.  It returns each fasta one by 
       as a tuple -> (name, sequence). '''
    file_iter = open(file_name)
    current_seq = [] # a dummy, needed to get through the 1st read only
    for line in file_iter:
        if not line.startswith('>'):
            current_seq.append(line.strip())
        else:
            if len(current_seq) != 0:
                yield (current_name, ''.join(current_seq))
            current_name = line[1:].strip()
            current_seq = []
    yield (current_name, ''.join(current_seq))

print_fasta = lambda s: ('>'+i+'\n' + s[i] for i in s)

###Set functions
def intersection(sets):
    """Get the intersection of all input sets"""
    if all((type(i)==type(set()) for i in sets)):
        return reduce(set.intersection, sets)
    else:
        sets = map(set, sets)
        return reduce(set.intersection, sets)

def union(sets):
    """Get the union of all input sets"""
    if all((type(i)==type(set()) for i in sets)):
        return reduce(set.union, sets)
    else:
        sets = map(set, sets)
        return reduce(set.union, sets)

def join(seqs):
    """Join any input sequences that support concatenation"""
    return reduce(operator.concat, seqs)

#Misc
def get_file(filename, splitchar = 'NA', buffered = False):
    if not buffered:
        if splitchar == 'NA':
            return [i.strip().split() for i in open(filename)]
        else: return [i.strip().split(splitchar) for i in open(filename)]
    else:
        if splitchar == 'NA':
            return (i.strip().split() for i in open(filename))
        else: return (i.strip().split(splitchar) for i in open(filename))

def sort_dict_by_val(aDict):
    '''returns a list of tuples sorted by the dict values'''
    return sorted(aDict.iteritems(), key=lambda (k,v): (v,k))

def pairwise(li):  
    '''a convienience function that produces all pairwise comparisons from a list'''
    for i in range(len(li)):
        j = i+1
        while j < len(li):
            yield (li[i], li[j])
            j += 1

class Hash(defaultdict):
    '''works like a perl hash, auto-initializing
       Note: be careful with "if x in Hash"'''
    def __init__(self):
        defaultdict.__init__(self, Hash)
    def __reduce__(self):
        r = defaultdict.__reduce__(self)
        # override __init__ args
        return (r[0], (), r[2], r[3], r[4]) 

def count_all(xlist, proportions=False):
    '''Count all the items in a list, return a dict
       with the item as key and counts as value.
       If proportions are set to True, the values 
       are the proportions not counts'''
    out =  defaultdict(int)
    for i in xlist: out[i]+=1
    if proportions:
        out2 = {}
        tot_sz = float(sum(out.values()))
        for i in out: out2[i] = out[i] / tot_sz
        return out2
    else: return out

###Combinatorial functions
def product(*args, **kwds):
    ''' product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
        product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111'''
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def permutations(iterable, r=None):
    ''' permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
        permutations(range(3)) --> 012 021 102 120 201 210 '''
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    for indices in product(range(n), repeat=r):
        if len(set(indices)) == r:
            yield tuple(pool[i] for i in indices)

def combinations(iterable, r):
    ''' combinations('ABCD', 2) --> AB AC AD BC BD CD
        combinations(range(4), 3) --> 012 013 023 123 '''
    pool = tuple(iterable)
    n = len(pool)
    for indices in permutations(range(n), r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

def combinations_with_replacement(iterable, r):
    '''combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC'''
    pool = tuple(iterable)
    n = len(pool)
    for indices in product(range(n), repeat=r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)

