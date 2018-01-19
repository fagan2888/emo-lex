# william la cava 2018
# bound, normalize and filter results
# compute hypervolume indicators
# store as csv file for plotting

import numpy as np
import os 

problems = ["dtlz1_m5","dtlz1_m25","dtlz1_m50","dtlz1_m75","dtlz1_m100" 
            "dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100" 
            "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100"
            "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]

selectors = ["lex","nsga2","hype"] 

g='200'

X = []
c=0

for p in problems:
    res_files = [p + '_' + s + '.' + g for s in selectors]
    all_res = p + '.' +  g
    os.system("cat {P}_*.{G} > {RES}".format(P=p,G=g,RES=all_res))

    with open(all_res,'r') as f:
        for lines in f:
            if lines!='\n':        
                x = [float(i) for i in lines.strip().split(' ')]            
                X.append(x)   
    max_obj = np.max(np.max([i for i in X]))
    min_obj = np.min(np.min([i for i in X]))

    #normalize results
    for rf in res_files:
        with open(rf,'r') as f:
            for lines in f:
                if lines!='\n':
                    x = [(float(i)-min_obj)/(max_obj-min_obj)+1.0 for i in lines.strip().split(' ')]
                    out = [str(i) + ' ' for i in x] 
                else: 
                    out = '\n'
                with open(rf+'.norm','w') as f:
                    f.write(out + '\n')


