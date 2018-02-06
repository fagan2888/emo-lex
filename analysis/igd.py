# calculates the inverted generational distance of a set of results compared to a reference pareto
# front. 
# william la cava 2018

import numpy as np
from sklearn.neighbors import NearestNeighbors
import pandas as pd

def igd(X, problem):
    # X: data
    # problem: name of problem
    X = np.asarray(X)
    PF = pd.read_csv('pareto_fronts/'+problem +'.pf',header=None).values
    m = int(problem.split('m')[-1])
    #print('PF shape:',PF.shape, 'X shape:',X.shape)
    nbrs= NearestNeighbors(n_neighbors=1,algorithm='ball_tree').fit(X)

    d,i = nbrs.kneighbors(PF)
    assert(len(d)==PF.shape[0])    
    return np.mean(d)
