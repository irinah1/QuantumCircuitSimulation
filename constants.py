import numpy as np
import cmath

global sigma1
global sigma2
global sigma3
global I
 
sigma1 = np.matrix([[0,1],[1,0]])
sigma2 = np.matrix([[0,-1j],[1j,0]])
sigma3 = np.matrix([[1,0],[0,-1]])
I = np.asmatrix(np.identity(2))

def X(dim):
    return np.asmatrix(np.concatenate((np.concatenate((np.zeros((1,dim-1), dtype=complex), np.identity(1, dtype=complex)),axis=1), np.concatenate((np.identity(dim-1,dtype=complex), np.zeros((dim-1,1),dtype=complex)),axis=1))))
    
def Z(dim):
    mat = np.zeros((dim,dim), dtype=complex)
    np.fill_diagonal(mat, np.power(np.exp(2*np.pi*1j/dim), np.arange(dim)))
    return np.asmatrix(mat)
