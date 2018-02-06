# generates pareto front reference sets for problems. 
# william la cava 2018
import numpy as np
from tqdm import tqdm
from math import isclose
import sys

problems = ["dtlz1_m3", "dtlz1_m5","dtlz1_m25","dtlz1_m50","dtlz1_m75","dtlz1_m100" ,
            "dtlz2_m3", "dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100" ,
            "dtlz3_m3", "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100",
            "dtlz4_m3", "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]
if len(sys.argv)>1:
    problems = [sys.argv[1]]
# set up weight vectors w evenly distributed in objective space

# then f_i(x) = 0.5 * w_i / sum(w) for dtlz1
# and f_i(x) = w_i / sqrt(sum(w)) for dtlz2-4

n_samples = 10000
np.random.seed(42)
for p in problems:
    d = int(p.split('m')[1])
    print(p)
    out_file = 'pareto_fronts/' + p + '.pf'

    for n in tqdm(np.arange(n_samples)):
        # randomly sample weights
        w = np.random.rand(d)

        if 'dtlz1' in p:
            f = np.array([0.5 * w_i / np.sum(w) for w_i in w])
            
            if not isclose(np.sum(f),0.5,abs_tol=10**-6):
                print('f:',f,'sum(f):',np.sum(f))
   
            assert (isclose(np.sum(f),0.5,abs_tol=10**-6))
        else:
            f = np.array([ w_i / np.sqrt(np.sum([w_j**2 for w_j in w])) for w_i in w])
            if not isclose(np.sum([f_i**2 for f_i in f]),1.0,abs_tol=10**-6):
                print('f:',f,'sum f^2 :', np.sum([f_i**2 for f_i in f]))
            assert (isclose(np.sum([f_i**2 for f_i in f]),1.0,abs_tol=10**-6))
            

        o = 'w' if n == 0 else 'a'
        with open(out_file,o) as out:
            out.write(','.join([str(f_i) for f_i in f]) + '\n')


