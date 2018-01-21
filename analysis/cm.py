# script to calculate the convergence measure for a set of dtlz experiments.
# william la cava 2018

import numpy as np
import sys

# inputs: 1) file, 2) output file (default: appends _cm to file)

#in_file=sys.argv[1]
#filename = in_file.split('/')[-1]
#out_file=''.join(['/'.join(in_file.split('/')[:-1]),'/'] + filename.split('.')[:-1]+
#                    ['_ave.',filename.split('.')[-1]])

def cm(X,one=False):
    if one:
        return np.mean([np.sum([i for i in x])-0.5 for x in X])
    else:
        return np.mean([np.abs(np.sum([i**2 for i in x])-1.0) for x in X])
#X = []
#c=0
#cms = []
#with open(in_file,'r') as f:
#    for lines in f:
#        if lines!='\n':        
#            x = [float(i) for i in lines.strip().split(' ')]            
#            X.append(x)
#        else:
#            c+=1
#            #print('set ',c)
#            one = 'dtlz1' in in_file
#            result = cm(X,one) 
#            
#            with open(out_file,'a') as out:
#                out.write(str(result)+'\n')
#            cms.append(result)
#
#print(c,' trials processed. output written to ', out_file)
#print('cms:',cms)
