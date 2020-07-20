import numpy as np
import cmath
from constants import *

def tensorProdtoMatrix2D(mat1, mat2):
    """
    calculates tensor product in matrix form
    """
    matrix = np.matrix([[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j],[0.+0.*1j,0.+0.*1j,0.+0.*1j,0.+0.*1j]])
    for i in range(mat1.shape[0]):
        for j in range(mat1.shape[1]):
            for k in range(mat2.shape[0]):
                for l in range(mat2.shape[1]):
                    matrix[2*i+k, 2*j+l] = mat1[i,j]*mat2[k,l]
    return matrix

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
            res_mat += coef_matrix[i,j]*tensorProdtoMatrix2D(mat1, mat2)
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
    pauli_basis_mat = np.matrix([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
    
    pauli_basis_mat[0,0] = np.real(0.25*(res[0,0] + res[1,1] + res[2,2] + res[3,3]))    #a00
    pauli_basis_mat[0,1] = np.real(0.25*(res[0,1] + res[1,0] + res[2,3] + res[3,2]))    #a01
    pauli_basis_mat[0,2] = np.real(0.25*1j*(res[0,1] - res[1,0] + res[2,3] - res[3,2])) #a02
    pauli_basis_mat[0,3] = np.real(0.25*(res[0,0] - res[1,1] + res[2,2] - res[3,3]))    #a03
    
    pauli_basis_mat[1,0] = np.real(0.25*(res[0,2] + res[1,3] + res[2,0] + res[3,1]))    #a10
    pauli_basis_mat[1,1] = np.real(0.25*(res[0,3] + res[1,2] + res[2,1] + res[3,0]))    #a11
    pauli_basis_mat[1,2] = np.real(0.25*1j*(res[0,3] - res[1,2] + res[2,1] - res[3,0])) #a12
    pauli_basis_mat[1,3] = np.real(0.25*(res[0,2] - res[1,3] + res[2,0] - res[3,1]))    #a13
    
    pauli_basis_mat[2,0] = np.real(0.25*1j*(res[0,2] + res[1,3] - res[2,0] - res[3,1])) #a20
    pauli_basis_mat[2,1] = np.real(0.25*1j*(res[0,3] + res[1,2] - res[2,1] - res[3,0])) #a21
    pauli_basis_mat[2,2] = np.real(0.25*(-res[0,3] + res[1,2] + res[2,1] - res[3,0]))   #a22
    pauli_basis_mat[2,3] = np.real(0.25*1j*(res[0,2] - res[1,3] - res[2,0] + res[3,1])) #a23
    
    pauli_basis_mat[3,0] = np.real(0.25*(res[0,0] + res[1,1] - res[2,2] - res[3,3]))    #a30
    pauli_basis_mat[3,1] = np.real(0.25*(res[0,1] + res[1,0] - res[2,3] - res[3,2]))    #a31
    pauli_basis_mat[3,2] = np.real(0.25*1j*(res[0,1] - res[1,0] - res[2,3] + res[3,2])) #a32
    pauli_basis_mat[3,3] = np.real(0.25*(res[0,0] - res[1,1] - res[2,2] + res[3,3]))    #a33
    
    return pauli_basis_mat