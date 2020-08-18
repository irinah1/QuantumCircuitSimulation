import numpy as np
import cmath
from constants import *

def get_rho_from_Pauli_basis(coef_matrix):
    """
    input:
    Pauli basis coefficient matrix
    
    output:
    rho
    """
    res_mat = np.matrix([[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j], [0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j]])
    i = 0
    for mat1 in [I, sigma1, sigma2, sigma3]:
        j = 0
        for mat2 in [I, sigma1, sigma2, sigma3]:
            res_mat += coef_matrix[i,j]*np.kron(mat1, mat2)
            j += 1
        i += 1
    return res_mat

def get_Pauli_basis_from_rho(res):
    """
    input:
    rho = res
    
    output:
    Pauli basis coefficient matrix
    """
    # bring into pauli basis (16 real coefficients)
    pauli_basis_mat = np.matrix([[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j]])

    Sigma = [I, sigma1, sigma2, sigma3]

    for i in range(len(Sigma)):
        for j in range(len(Sigma)):
            pauli_basis_mat[i,j] = 0.25*np.trace(np.kron(Sigma[i], Sigma[j])*res)
    
    return pauli_basis_mat
