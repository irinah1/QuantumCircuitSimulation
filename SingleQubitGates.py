import numpy as np
import cmath
from constants import *
    
# each operation takes a Bloch vector as input and outputs a new Bloch vector
    
# quantum gates

def I_gate(blochvector):
    #return operator_sum(blochvector, [I])
    return blochvector
    
def X_gate(blochvector):
    #return operator_sum(blochvector, [np.matrix([[0,1],[1,0]])])
    return np.array([blochvector[0], -blochvector[1], -blochvector[2]])
    
def Y_gate(blochvector):
    #return operator_sum(blochvector, [np.matrix([[0,-1j],[1j,0]])])
    return np.array([-blochvector[0], blochvector[1], -blochvector[2]])
    
def Z_gate(blochvector):
    #return operator_sum(blochvector, [np.matrix([[1,0],[0,-1]])])
    return np.array([-blochvector[0], -blochvector[1], blochvector[2]])
    
def H_gate(blochvector):
    #return operator_sum(blochvector, [(1/np.sqrt(2)*np.matrix([[1,1],[1,-1]]))])
    return np.array([blochvector[2], -blochvector[1], blochvector[0]])

def S_gate(blochvector):
    return operator_sum(blochvector, [np.matrix([[1,0],[0,np.exp(1j*np.pi/2)]])])
    
def T_gate(blochvector):
    return operator_sum(blochvector, [np.matrix([[1,0],[0,np.exp(1j*np.pi/4)]])])
    
# noise operations
    
def amplitude_damping(blochvector, Gamma, t):
    p=1
    E0 = np.sqrt(p)*np.matrix([ [1,     0          ],
                                [0,np.exp(-0.5*Gamma*t)]])
    E1 = np.sqrt(p)*np.matrix([ [0,np.sqrt(1-np.exp(-Gamma*t))],
                                [0,     0        ]])
    # for general amplitude damping
    #E2 = np.sqrt(1-p)*np.matrix([[np.sqrt(1-gamma),0],
    #                             [      0         ,1]])
    #E3 = np.sqrt(1-p)*np.matrix([[      0       ,0],
    #                             [np.sqrt(gamma),0]])
    return operator_sum(blochvector, [E0, E1])
    
def phase_damping(blochvector, Gamma, dw, t):
    E0 = np.matrix([ [1,                 0               ],
                     [0,np.exp(-Gamma*t)*np.exp(-1j*dw*t)]])
    E1 = np.matrix([ [0,               0              ],
                     [0, np.sqrt(1-np.exp(-2*Gamma*t))]])
    return operator_sum(blochvector, [E0, E1])
    

def operator_sum(blochvector, list):
    x = 0
    y = 0
    z = 0
    res = 0
    for mat in list:
        if type(mat)!=np.matrix:
            print("Cannot dagger non numpy.matrix type")
        else:
            res = res + mat*0.5*(I + blochvector[0]*sigma1 + blochvector[1]*sigma2 + blochvector[2]*sigma3)*(mat.H)
    x = np.real(2*res[0,1])
    y = -np.imag(2*res[0,1])
    z = 2*res[0,0] - 1
    return x, y, z
