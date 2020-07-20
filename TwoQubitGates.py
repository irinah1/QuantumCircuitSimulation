import numpy as np
import cmath
from constants import *
from helper_functions_2Qubits import *

def I2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    return coef_matrix

def X2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    factor_matrix = np.matrix([[1,1,1,1],[1,1,1,1,],[-1,-1,-1,-1],[-1,-1,-1,-1]])
    if qubitnumber == 2:
        factor_matrix = factor_matrix.T
    
    for i in range(coef_matrix.shape[0]):
        for j in range(coef_matrix.shape[1]):
            coef_matrix[i,j] = coef_matrix[i,j]*factor_matrix[i,j]
    return coef_matrix

def Y2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    factor_matrix = np.matrix([[1,1,1,1],[-1,-1,-1,-1,],[1,1,1,1],[-1,-1,-1,-1]])
    if qubitnumber == 2:
        factor_matrix = factor_matrix.T
    
    for i in range(coef_matrix.shape[0]):
        for j in range(coef_matrix.shape[1]):
            coef_matrix[i,j] = coef_matrix[i,j]*factor_matrix[i,j]
    return coef_matrix
        
def Z2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    factor_matrix = np.matrix([[1,1,1,1],[-1,-1,-1,-1,],[-1,-1,-1,-1],[1,1,1,1]])
    if qubitnumber == 2:
        factor_matrix = factor_matrix.T
    
    for i in range(coef_matrix.shape[0]):
        for j in range(coef_matrix.shape[1]):
            coef_matrix[i,j] = coef_matrix[i,j]*factor_matrix[i,j]
    return coef_matrix

def H2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    trans_matrix = np.matrix([[1,0,0,0],[0,0,0,1],[0,0,-1,0],[0,1,0,0]])
    if qubitnumber == 2:
        coef_matrix = coef_matrix*trans_matrix
    else:
        coef_matrix = trans_matrix*coef_matrix
    
    return coef_matrix

def S2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    trans_matrix = np.matrix([[1,0,0,0],[0,0,-1,0],[0,1,0,0],[0,0,0,1]])
    if qubitnumber == 2:
        coef_matrix = coef_matrix*trans_matrix
    else:
        coef_matrix = trans_matrix*coef_matrix
    
    return coef_matrix

def T2_gate(coef_matrix, qubitnumber):
    """
    qubitnumber element {1,2}
    """
    if qubitnumber<1 or qubitnumber>2:
        print("Qubit number must be either 1 or 2!")
    if qubitnumber == 2:
        trans_matrix = tensorProdtoMatrix2D(I,np.matrix([[1,0],[0,np.exp(1j*np.pi*0.25)]]))
    else:
        trans_matrix = tensorProdtoMatrix2D(np.matrix([[1,0],[0,np.exp(1j*np.pi*0.25)]]),I)
    rho = get_rho_from_Pauli_basis(coef_matrix)
    rho = operator_sum(rho, [trans_matrix])
    return get_Pauli_basis_from_rho(rho)

def CNOTgate(coef_matrix):
    rho = get_rho_from_Pauli_basis(coef_matrix)
    rho = operator_sum(rho, [np.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])])
    return get_Pauli_basis_from_rho(rho)

def CZgate(coef_matrix):
    rho = get_rho_from_Pauli_basis(coef_matrix)
    rho = operator_sum(rho, [np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])])
    return get_Pauli_basis_from_rho(rho)

def SWAPgate(coef_matrix):
    rho = get_rho_from_Pauli_basis(coef_matrix)
    rho = operator_sum(rho, [np.matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])])
    return get_Pauli_basis_from_rho(rho)

def operator_sum(rho, list):
    res = 0
    for mat in list:
        if type(mat)!=np.matrix:
            print("Cannot dagger non numpy.matrix type")
        else:
            res = res + mat*rho*(mat.H)
    
    return res
