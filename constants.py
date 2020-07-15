import numpy as np
import cmath

global sigma1
global sigma2
global sigma3
global I
 
sigma1 = np.matrix([[0,1],[1,0]])
sigma2 = np.matrix([[0,-1j],[1j,0]])
sigma3 = np.matrix([[1,0],[0,-1]])
I = np.matrix([[1,0],[0,1]])
