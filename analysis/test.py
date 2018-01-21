import numpy as np
X = []
with open('../code/runs/dtlz2_m25_lex_.1000','r') as f:
    for lines in f:
        if lines != '\n':
            print(lines)
            X.append([float(x) for x in lines.strip().split(' ')])
X = np.array(X)        
print('X shape:' , X.shape)

