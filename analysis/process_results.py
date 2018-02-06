# william la cava 2018
# collect results
# calculate convergence measure and inverted generation distance
# store as csv file for plotting

import numpy as np
import os 
from cm import cm
from igd import igd

data_paths = (['../code/runs/']+ ['../code/trial_'+str(i)+'/runs/' for i in np.arange(30)+1] +
             ['../code/runs/old_runs_g1000/'])
root = "dtlz"
probs=[1,2,3,4]
ms=[3,5,25,50,75,100]
problems = [root+str(p)+'_m'+str(m) for p in probs for m in ms]
print('problems:',problems)
#problems = ["dtzl1_m3","dtlz1_m5","dtlz1_m25","dtlz1_m50","dtlz1_m75","dtlz1_m100", 
#           "dtzl2_m3", "dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100", 
#           "dtlz3_m3", "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100",
#            "dtlz4_m3", "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]
#problems = ["dtlz2_m5","dtlz2_m25","dtlz2_m50","dtlz2_m75","dtlz2_m100", 
#            "dtlz3_m5","dtlz3_m25","dtlz3_m50","dtlz3_m75","dtlz3_m100",
#            "dtlz4_m5","dtlz4_m25","dtlz4_m50","dtlz4_m75","dtlz4_m100"]


selectors = ["lex","nsga2","hype"] 

g='1000'

out_file = 'summary_' +'_'.join(selectors) + '_new.csv'
X = []
c=0
trials_per_exp = {}

# write header
with open(out_file,'w') as out:
    out.write('dataset,m,method,trial,cm,igd\n')

for p in problems:
    m = float(p.split('m')[-1])
    if 'dtlz2' in problems:
        if m == 3 or m == 5:
            g = 500
        elif m == 25:
            g = 1000
        elif m == 50:
            #g = 2000
            g = 800
        elif m == 75:
            #g = 3000
            g = 1000
        elif m == 100:
            #g = 4000
            g = 1000
    else:
        if m == 3: 
            g = 500
        elif m == 5:
            g = 1000
        elif m == 25:
            #g = 2000
            g = 1000
        elif m == 50:
            #g = 3000
            g = 1000
        elif m == 75:
            #g = 4000
            g = 1000
        elif m == 100:
            #g = 5000
            g = 1000
    g = str(g)
    for s in selectors: 
        trial = 1
        for d in data_paths:
            rf = d + p + '_' + s + '.' + g
            print('===',rf,'===')
#            if trial>35: continue
            try:
                try:
                    f = open(rf,'r')
                except IOError:
                    f = open(d + p + '_' + s + '_.' + g)

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
                            print('',output,end='')                     
                            trial = trial + 1
                            X = []
                            print('trial',trial)
                         
            except IOError:    
                print('','!! cant open',rf)
        trials_per_exp.update({p+'_'+s:trial})


for k,v in sorted(trials_per_exp.items()):
    
    if v < 30:
        print(k,v, '(need more trials)')
    else:
        print(k,v)

