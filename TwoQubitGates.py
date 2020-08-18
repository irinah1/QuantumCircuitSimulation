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
        trans_matrix = np.kron(I,np.matrix([[1,0],[0,np.exp(1j*np.pi*0.25)]]))
    else:
        trans_matrix = np.kron(np.matrix([[1,0],[0,np.exp(1j*np.pi*0.25)]]),I)
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

def individual_amplitude_damping(coef_matrix, Gamma, t):
    # performs amplitude damping for each qubit individually
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E10 = (1/np.sqrt(2))*np.kron(np.matrix([[0,np.sqrt(1-np.exp(-Gamma*t))],[0,0]]), I)
    E01 = (1/np.sqrt(2))*np.kron(I, np.matrix([[0,np.sqrt(1-np.exp(-Gamma*t))],[0,0]]))
    E00 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E01.H)*E01-(E10.H)*E10)
    #print("E01:\n",E01*rho*E01.H, "\n")
    #print("E10:\n",E10*rho*E10.H, "\n")
    #print("E00:\n",E00*rho*E00.H, "\n")
    rho = operator_sum(rho, [E10, E01, E00])
    return get_Pauli_basis_from_rho(rho)

def fully_correlated_amplitude_damping(coef_matrix, Gamma, t):
    # performs correlated amplitude damping
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E1 = np.kron(np.matrix([[0,np.sqrt(1-np.exp(-Gamma*t))],[0,0]]), np.matrix([[0,np.sqrt(1-np.exp(-Gamma*t))],[0,0]]))
    E0 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E1.H)*E1)

    rho = operator_sum(rho, [E1, E0])
    return get_Pauli_basis_from_rho(rho)

def individual_phase_damping(coef_matrix, Gamma, t):
    # performs phase damping for each qubit individually
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E10 = (1/np.sqrt(2))*np.kron(np.matrix([[0,0],[0,np.sqrt(1-np.exp(-2*Gamma*t))]]), I)
    E01 = (1/np.sqrt(2))*np.kron(I, np.matrix([ [0,0],[0,np.sqrt(1-np.exp(-2*Gamma*t))]]))
    E00 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E01.H)*E01-(E10.H)*E10)
    
    rho = operator_sum(rho, [E10, E01, E00])
    return get_Pauli_basis_from_rho(rho)

def fully_correlated_phase_damping(coef_matrix, Gamma, t):
    # performs correlated phase damping
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E1 = np.kron(np.matrix([[0,0],[0,np.sqrt(1-np.exp(-2*Gamma*t))]]), np.matrix([[0,0],[0,np.sqrt(1-np.exp(-2*Gamma*t))]]))
    E0 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E1.H)*E1)
    
    rho = operator_sum(rho, [E1, E0])
    return get_Pauli_basis_from_rho(rho)

def individual_phase_flip(coef_matrix, p):
    # performs phase flip for each qubit individually
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E10 = (1/np.sqrt(2))*np.kron(np.sqrt(p)*np.matrix([[1,0],[0,-1]]), I)
    E01 = (1/np.sqrt(2))*np.kron(I, np.sqrt(p)*np.matrix([[1,0],[0,-1]]))
    E00 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E01.H)*E01-(E10.H)*E10)
    
    rho = operator_sum(rho, [E10, E01, E00])
    return get_Pauli_basis_from_rho(rho)

def fully_correlated_phase_flip(coef_matrix, p):
    # performs correlated phase flip
    rho = get_rho_from_Pauli_basis(coef_matrix)
    
    E1 = (1/np.sqrt(2))*np.kron(np.sqrt(p)*np.matrix([[1,0],[0,-1]]), np.sqrt(p)*np.matrix([[1,0],[0,-1]]))
    E0 = np.sqrt(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])-(E1.H)*E1)
    
    rho = operator_sum(rho, [E1, E0])
    return get_Pauli_basis_from_rho(rho)

def operator_sum(rho, list):
    res = 0
    for mat in list:
        if type(mat)!=np.matrix:
            print("Cannot dagger non numpy.matrix type")
        else:
            res = res + mat*rho*(mat.H)
    
    return res
