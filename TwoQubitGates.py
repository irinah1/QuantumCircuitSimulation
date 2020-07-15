import numpy as np
import cmath
from constants import *

def operator_sum(rho, list):
res = 0
for mat in list:
    if type(mat)!=np.matrix:
        print("Cannot dagger non numpy.matrix type")
    else:
        res = res + mat*rho*(mat.H)
return res
