import numpy as np
import cmath
from constants import *

class Qudit:
    """
    input:
    d = dimension
    r = [r1, r2, r3, ...]: Qudit
    """
    def __init__(self, d, vec):
        
        self.dimension = d
        scale = 1.
        
        if type(vec)==list:
            vec = np.array(vec)
            
        if (np.sum(vec**2) > 1):
            print("Length of vector r ist scaled to 1.")
            scale = 1/np.sqrt(np.sum(vec**2))
            
        vec = np.asmatrix(vec)
        if vec.shape != (1, self.dimension):
            raise ValueError(" Dimension must match array shape or list length!")
        self.matrix = scale**2*vec.H*vec
    
    def get_matrix(self):
        return self.matrix
    
    def _set_matrix(self, mat):
        if type(mat)!=np.matrix:
            mat = np.asmatrix(mat)
        self.matrix = mat
        
    def get_dimension(self):
        return self.dimension
        
        
