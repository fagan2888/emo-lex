# william la cava 2018
# collect results
# calculate convergence measure and inverted generation distance
# store as csv file for plotting

import numpy as np
import os 
from cm import cm
from igd import igd

data_path = '../code/runs/'
problems = ["dtlz1_m5","dtlz1_m25","dtlz1_m50","dtlz1_m75","dtlz1_m100", 
            "dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100", 
           "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100",
            "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]
#problems = ["dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100", 
#            "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100",
#            "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]


selectors = ["lex","nsga2","hype"] 

g='1000'
out_file = 'summary_' +'_'.join(selectors) + '_' + g + '.csv'
X = []
c=0

# write header
with open(out_file,'w') as out:
    out.write('dataset,m,method,trial,cm,igd\n')

for p in problems:
    for s in selectors:
        try:
            #print(p,s)
            rf = data_path + p + '_' + s + '.' + g
            trial = 0
            try: 
                f = open(rf,'r')
            except IOError:
                rf = data_path + p + '_' + s + '_.' + g
                f = open(rf,'r')
            
            print('===\n',rf,'===\n')
            with f:
                for lines in f:
                    if lines!='\n':        
                        x = [float(i) for i in lines.strip().split(' ')]            
                        X.append(x)
                    else:
                        c+=1
                        #print('set ',c)
                        one = 'dtlz1' in p
                        CM = cm(X,one)
                        IGD = igd(X,p)
                        dataset = p.split('_')[0]
                        m = 'm=' + p.split('m')[-1]

                        output = ','.join([dataset, m, s, str(trial), str(CM), str(IGD)]) + '\n'
                        
                        with open(out_file,'a') as out:
                            out.write(output)
                        print(output,end='')                     
                        trial = trial + 1
                        X = []
        except:
            print(p,s,'failed')
