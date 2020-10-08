import numpy as np
import cmath
from constants import *

class Qudit:
    """
    input:
    d = dimension
    r = [r1, r2, r3, ...]: Qudit
    r = [r1, r2, r3, ...]: Qudit
    """
    def __init__(self, d, vec1, vec2):
        
        self.dimension = d
        scale1 = 1.
        scale2 = 1.
        sum = 0
        for i in range(d):
            sum += vec1[i]**2
        if (sum > 1):
            #print("Blochvector r1 has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
            scale1 = 1/np.sqrt(sum)
        for i in range(d):
            sum += vec2[i]**2
        if (sum > 1):
            #print("Blochvector r1 has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
            scale2 = 1/np.sqrt(sum)
        if type(vec1)==list:
            vec = np.array(vec1)
        if type(vec2)==list:
            vec = np.array(vec2)
            
        vec = np.asmatrix(vec)
        if vec.shape != (1, self.dimension):
            raise ValueError(" Dimension must match array shape or list length!")
        
        matrix1 = scale1**2*vec1.H*vec1
        matrix2 = scale2**2*vec2.H*vec2
        self.matrix = np.kron(matrix1, matrix2)
    
    def get_matrix(self):
        return self.matrix
    
    def _set_matrix(self, mat):
        if type(mat)!=np.matrix:
            mat = np.asmatrix(mat)
        self.matrix = mat
        
    def get_dimension(self):
        return self.dimension
        
        

