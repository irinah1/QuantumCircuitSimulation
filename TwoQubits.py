import numpy as np
import cmath
from constants import *

class TwoQubits:
    """
    input:
    
    r = [rx, ry, rz]: Blochvector entries: rho = 0.5*(I + r*sigma_vector)
    """
    def __init__(self, vec1, vec2):
        scale1 = 1.
        scale2 = 1.
        if (vec1[0]**2+vec1[1]**2+vec1[2]**2 > 1):
            print("Blochvector r1 has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
            scale1 = 1/np.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
        if (vec2[0]**2+vec2[1]**2+vec2[2]**2 > 1):
            print("Blochvector r2 has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
            scale2 = 1/np.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
        if type(vec1)==list:
            vec = np.array(vec1)
        if type(vec2)==list:
            vec = np.array(vec2)
        
        matrix1 = 0.5*(I + scale1*vec1[0]*sigma1 + scale1*vec1[1]*sigma2 + scale1*vec1[2]*sigma3)
        matrix2 = 0.5*(I + scale2*vec2[0]*sigma1 + scale2*vec2[1]*sigma2 + scale2*vec2[2]*sigma3)
    
        self.matrix = np.matrix([[matrix1[0,0]*matrix2[0,0], matrix1[0,0]*matrix2[0,1],                                                matrix1[0,1]*matrix2[0,0], matrix1[0,1]*matrix2[0,1]],
    
                                 [matrix1[0,0]*matrix2[1,0], matrix1[0,0]*matrix2[1,1],
                                  matrix1[0,1]*matrix2[1,0], matrix1[0,1]*matrix2[1,1]],
                              
                                 [matrix1[1,0]*matrix2[0,0], matrix1[1,0]*matrix2[0,1],
                                 matrix1[1,1]*matrix2[0,0], matrix1[1,1]*matrix2[0,1]],
                              
                                 [matrix1[1,0]*matrix2[1,0], matrix1[1,0]*matrix2[1,1],
                                  matrix1[1,1]*matrix2[1,0], matrix1[1,1]*matrix2[1,1]]])
    def get_matrix(self):
        return self.matrix
