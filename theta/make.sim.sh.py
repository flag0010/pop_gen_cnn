from random import *
import json
import numpy as np

v = np.random.uniform(10, 50, 5000)
#x = range(10,101)
#v = [choice(x) for i in xrange(40000)]
json.dump([i for i in v], open('true.thetas.json', 'w'))

for i in v:
    print "./ms 40 1 -t " + str(i) + " >> theta.sims.txt"
#print v[:10]
