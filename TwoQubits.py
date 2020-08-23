import numpy as np
import cmath
from constants import *
from helper_functions_2Qubits import *

class TwoQubits:
    """
    input:
    
    r = [rx, ry, rz]: Blochvector entries: rho = 0.5*(I + r*sigma_vector)
    
    uses:
                  a00 ... a03
    coef_matrix =  .  .    .
                   .     . .
                  a30 ... a33
              
    where rho = a00*(I⊗I) + a01*(I⊗sigmaX) + ... + a33*(sigmaZ⊗sigmaZ)
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
        
        #matrix1 = 0.5*(I + scale1*vec1[0]*sigma1 + scale1*vec1[1]*sigma2 + scale1*vec1[2]*sigma3)
        #matrix2 = 0.5*(I + scale2*vec2[0]*sigma1 + scale2*vec2[1]*sigma2 + scale2*vec2[2]*sigma3)
    
        self.coef_matrix = 0.25*np.matrix([[   1,        vec2[0],         vec2[1],          vec2[2]],
                                        [vec1[0],vec1[0]*vec2[0], vec1[0]*vec2[1],  vec1[0]*vec2[2]],
                                        [vec1[1],vec1[1]*vec2[0], vec1[1]*vec2[1],  vec1[1]*vec2[2]],
                                        [vec1[2],vec1[2]*vec2[0], vec1[2]*vec2[1],  vec1[2]*vec2[2]]])
        
    
    def get_matrix(self):
        return get_rho_from_Pauli_basis(self.coef_matrix)
    
    def set_Pauli_basis_matrix(self, matrix):
        self.coef_matrix = matrix
