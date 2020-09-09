import numpy as np
import cmath
import math
from constants import *
from SingleQudit import *
import scipy.special
    
# each operation takes a Bloch vector as input and outputs a new Bloch vector
    
# quantum gates

def I_gate(mat):
    #return operator_sum(blochvector, [I])
    return mat
    
def X_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 0
    gate_matrix[0,1] = 1
    gate_matrix[1,0] = 1
    gate_matrix[1,1] = 0
    return qudit_operator_sum(mat, [gate_matrix])
    
def Y_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 0
    gate_matrix[0,1] = -1j
    gate_matrix[1,0] = 1j
    gate_matrix[1,1] = 0
    return qudit_operator_sum(mat, [gate_matrix])
    
def Z_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 1
    gate_matrix[0,1] = 0
    gate_matrix[1,0] = 0
    gate_matrix[1,1] = -1
    return qudit_operator_sum(mat, [gate_matrix])
    
def H_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 1/np.sqrt(2)
    gate_matrix[0,1] = 1/np.sqrt(2)
    gate_matrix[1,0] = 1/np.sqrt(2)
    gate_matrix[1,1] = -1/np.sqrt(2)
    return qudit_operator_sum(mat, [gate_matrix])
    
def S_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 1
    gate_matrix[0,1] = 0
    gate_matrix[1,0] = 0
    gate_matrix[1,1] = np.exp(1j*np.pi/2)
    return qudit_operator_sum(mat, [gate_matrix])
    
def T_gate(mat):
    dim = mat.shape[0]
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = 1
    gate_matrix[0,1] = 0
    gate_matrix[1,0] = 0
    gate_matrix[1,1] = np.exp(1j*np.pi/4)
    return qudit_operator_sum(mat, [gate_matrix])
    
# noise operations
    
def general_bit_flip(mat, p):
    dim = mat.shape[0]
    operatorlist=[np.sqrt(p)*X(dim), np.sqrt(1-p)*np.asmatrix(np.identity(dim, dtype=complex))]
    return qudit_operator_sum(mat, operatorlist)

def general_phase_flip(mat, p):
    dim = mat.shape[0]
    operatorlist=[np.sqrt(p)*Z(dim), np.sqrt(1-p)*np.asmatrix(np.identity(dim, dtype=complex))]
    return qudit_operator_sum(mat, operatorlist)

def depolarization(mat, p):
    # applies Weyl channels with probability p
    dim = mat.shape[0]
    operatorlist = []
    for n in range(dim):
        for m in range(dim):
            operatorlist.append(np.sqrt(p)/dim*(X(dim)**n)*(Z(dim)**m))
    operatorlist.append(np.sqrt(1-p)*np.asmatrix(np.identity(dim, dtype=complex)))
    return qudit_operator_sum(mat, operatorlist)
    
def single_Weyl_channel(mat, n, m, p):
    """
    n, m : int : X^n Z^m
    p : float : probability of W_nm
    """
    dim = mat.shape[0]
    operatorlist = [np.sqrt(p)*(X(dim)**n)*(Z(dim)**m), np.sqrt(1-p)*np.asmatrix(np.identity(dim, dtype=complex))]
    return qudit_operator_sum(mat, operatorlist)

def amplitude_damping(mat, gammalist):
    """
    input: gammalist = [gamma_01, gamma_02, ..., gamma_0d-1,
                        gamma_12,..., gamma_1d-1,
                        ...,
                        gamma_(d-2)(d-1)]
    """
    dim = mat.shape[0]
    try:
        gammalen = len(gammalist)
    except:
        gammalen = gammalist.shape[0]
    if gammalen != round(0.5*dim*(dim-1)):
        raise ValueError("Dimension of gamma list {} must match         0.5*system_dimension*(system_dimension) = {}.".format(gammalen,round(0.5*dim*(dim-1))))
    
    operatorlist = []
    j = 0
    k = 1
    minus = (dim-2-j)
    for i in range(gammalen):
        if (i - minus > 0):
            minus += (dim-2-j)
            j += 1
            k = j+1
    
        # |j><i|
        basisvec1 = np.zeros(dim)
        basisvec1[k] = 1
        basisvec1 = np.asmatrix(basisvec1)
        basisvec2 = np.zeros(dim)
        basisvec2[j] = 1
        basisvec2 = np.asmatrix(basisvec2)
        operatorlist.append(np.sqrt(gammalist[i])*basisvec1.T*basisvec2)
        #print(i, ":", j, k)
        #print(np.sqrt(gammalist[i])*basisvec1.T*basisvec2)
        k += 1

    sum = np.asmatrix(np.zeros((dim,dim), dtype=complex))
    for op in operatorlist:
        sum += op.H*op
    operatorlist.append(np.sqrt(np.identity(dim, dtype=complex)-sum))

    return qudit_operator_sum(mat, operatorlist)

def general_phase_damping(mat, p):
    dim = mat.shape[0]
    operatorlist = []
    for n in range(dim):
        operatorlist.append(np.sqrt(scipy.special.binom(dim-1, n)*((0.5*(1-p))**n)*((0.5*(1+p))**(dim-1-n)))*(Z(dim)**n))
    return qudit_operator_sum(mat, operatorlist)
    
def physical_phase_damping(mat, gamma):
    #daming of off diagonal elements
    """
    gamma : damping strength exp(-gamma)
    """
    dim = mat.shape[0]
    for i in range(dim):
        for j in range(dim):
            mat[i,j] = mat[i,j]*((np.exp(-gamma))**((i-j)**2))
    return mat


def qudit_operator_sum(density, list):
    res = 0
    for mat in list:
        if type(mat)!=np.matrix:
            print("Cannot dagger non numpy.matrix type")
        else:
            res = res + mat*density*(mat.H)
    return res
