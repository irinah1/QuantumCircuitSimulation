import numpy as np
import cmath
from constants import *
from helper_functions_2Qubits import *
from TwoQubits import *
from TwoQubitGates import *
from SingleQudit import *
from SingleQuditGates import *
import random

def ideal_QFT_circuit(input_pauli_basis):
    """
    input:
    input_pauli_basis : coef. matrix in Pauli basis
    
    performs quantum Fuorier transform and outputs coef. matrix in Pauli basis
    """
    input_pauli_basis = H2_gate(input_pauli_basis, 1)
    input_pauli_basis = R2_gate(input_pauli_basis)
    input_pauli_basis = H2_gate(input_pauli_basis, 2)
    input_pauli_basis = SWAPgate(input_pauli_basis)

    return input_pauli_basis
    
def noisy_QFT_circuit(input_pauli_basis, set_seed1=1, set_seed2=234):
    """
    input:
    input_pauli_basis : coef. matrix in Pauli basis
    
    performs quantum Fuorier transform and outputs coef. matrix in Pauli basis
    """
    prob_list = get_probability_list(5, 0, 0.2/np.exp(1), set_seed=set_seed1)
    random.seed(set_seed2)
    
    # circuit
    number = random.randint(0,15)
    #print("Noise gate number ", number)
    input_pauli_basis = get_gate_from_list(input_pauli_basis, number, prob_list[0])
    input_pauli_basis = H2_gate(input_pauli_basis, 1)
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_pauli_basis = get_gate_from_list(input_pauli_basis, number, prob_list[1])
    input_pauli_basis = R2_gate(input_pauli_basis)
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_pauli_basis = get_gate_from_list(input_pauli_basis, number, prob_list[2])
    input_pauli_basis = H2_gate(input_pauli_basis, 2)
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_pauli_basis = get_gate_from_list(input_pauli_basis, number, prob_list[3])
    input_pauli_basis = SWAPgate(input_pauli_basis)
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_pauli_basis = get_gate_from_list(input_pauli_basis, number, prob_list[4])
    
    return input_pauli_basis
    
def leaked_noisy_QFT_circuit(input_density_matrix, set_seed1=1, set_seed2=234):
    """
    input:
    input_pauli_basis : coef. matrix in Pauli basis
    
    performs quantum Fuorier transform and outputs coef. matrix in Pauli basis
    """
    prob_list = get_probability_list(5, 0, 0.2/np.exp(1), set_seed=set_seed1)
    random.seed(set_seed2)
    dim = int(np.round(np.sqrt(input_density_matrix.shape[0])))
    # circuit
    # --- noise 1 ---
    number = random.randint(0,15)
    #print("Noise gate number ", number)
    input_density_matrix = get_qudit_gate_from_list(input_density_matrix, number, prob_list[0])
    # --- Hadamard Gate 1 ---
    input_density_matrix = qudit_operator_sum(input_density_matrix, [np.kron(qudit_operation_matrix(1/np.sqrt(2)*np.matrix([[1,1],[1,-1]]), dim), np.identity(dim,dtype=complex))])
    # --- noise 2 ---
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_density_matrix = get_qudit_gate_from_list(input_density_matrix, number, prob_list[1])
    # --- R2 Gate ---
    input_density_matrix = qudit_operator_sum(input_density_matrix, [two_qubit_to_qudit_gate(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1j]]), dim)])
    # --- noise 3 ---
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_density_matrix = get_qudit_gate_from_list(input_density_matrix, number, prob_list[2])
    # --- Hadamard Gate 2 ---
    input_density_matrix = qudit_operator_sum(input_density_matrix, [np.kron(np.identity(dim,dtype=complex), qudit_operation_matrix(1/np.sqrt(2)*np.matrix([[1,1],[1,-1]]), dim))])
    # --- noise 4 ---
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_density_matrix = get_qudit_gate_from_list(input_density_matrix, number, prob_list[3])
    # --- SWAP Gate ---
    input_density_matrix = qudit_operator_sum(input_density_matrix, [two_qubit_to_qudit_gate(np.matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]), dim)])
    # --- noise 5 ---
    number = random.randint(0, 15)
    #print("Noise gate number ", number)
    input_density_matrix = get_qudit_gate_from_list(input_density_matrix, number, prob_list[4])
    return input_density_matrix
    
def get_probability_list(length, min_val, max_val, set_seed=1):
    intervall = max_val - min_val
    # seed random number generator
    random.seed(set_seed)
    # generate random numbers between 0-1 * intervall
    checkpoint_list = []
    for i in range(length-1):
        checkpoint_list.append(intervall*random.random())
    checkpoint_list.sort()
    probability_list = []
    for x in range(len(checkpoint_list)):
        if x == 0:
            probability_list.append(min_val+checkpoint_list[x])
        else:
            probability_list.append(min_val+checkpoint_list[x]-checkpoint_list[x-1])
    probability_list.append(min_val+intervall-checkpoint_list[-1])
    #print(np.sum(np.array(probability_list)))
    return probability_list

def get_gate_from_list(coef_matrix, number, p):
    """
    coef_matrix : Pauli basis
    number      : number in list
    p           : probability
    """
    number = number % 16
    
    if number == 0:
        coef_matrix = (1-p)*coef_matrix + p*X2_gate(coef_matrix, 1)
    elif number == 1:
        coef_matrix = (1-p)*coef_matrix + p*X2_gate(coef_matrix, 2)
    elif number == 2:
        coef_matrix = (1-p)*coef_matrix + p*Y2_gate(coef_matrix, 1)
    elif number == 3:
        coef_matrix = (1-p)*coef_matrix + p*Y2_gate(coef_matrix, 2)
    elif number == 4:
        coef_matrix = (1-p)*coef_matrix + p*Z2_gate(coef_matrix, 1)
    elif number == 5:
        coef_matrix = (1-p)*coef_matrix + p*Z2_gate(coef_matrix, 2)
    elif number == 6:
        Gamma = np.sqrt(-np.log(1-p))
        coef_matrix = individual_amplitude_damping(coef_matrix, Gamma, Gamma)
    elif number == 7:
        Gamma = np.sqrt(-np.log(1-p))
        coef_matrix = fully_correlated_amplitude_damping(coef_matrix, Gamma, Gamma)
    elif number == 8:
        Gamma = np.sqrt(-np.log(1-p))
        coef_matrix = individual_phase_damping(coef_matrix, Gamma, Gamma)
    elif number == 9:
        Gamma = np.sqrt(-np.log(1-p))
        coef_matrix = fully_correlated_phase_damping(coef_matrix, Gamma, Gamma)
    elif number == 10:
        coef_matrix = individual_phase_flip(coef_matrix, p)
    elif number == 11:
        coef_matrix = fully_correlated_phase_flip(coef_matrix, p)
    elif number == 12:
        coef_matrix = individual_bit_flip(coef_matrix, p)
    elif number == 13:
        coef_matrix = fully_correlated_bit_flip(coef_matrix, p)
    elif number == 14:
        coef_matrix = individual_bit_phase_flip(coef_matrix, p)
    elif number == 15:
        coef_matrix = fully_correlated_bit_phase_flip(coef_matrix, p)
    
    return coef_matrix

def get_qudit_gate_from_list(density_matrix, number, p):
    """
    density_matrix : Pauli basis
    number      : number in list
    p           : probability
    """
    number = number % 30
    dim = int(np.round(np.sqrt(density_matrix.shape[0])))
    
    if number == 0:
        # X @ I
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma1, dim), np.identity(dim, dtype=complex)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2, dtype=complex))])
    elif number == 1:
        # I @ X
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), qudit_operation_matrix(sigma1, dim)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2, dtype=complex))])
    elif number == 2:
        # X @ X
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma1, dim), qudit_operation_matrix(sigma1, dim)), np.sqrt(1-p)*np.sqrt(np.identity(dim**2, dtype=complex)-np.kron(qudit_operation_matrix(sigma1, dim), qudit_operation_matrix(sigma1, dim)))])
    if number == 3:
        # Y @ I
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma2, dim), np.identity(dim, dtype=complex)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2, dtype=complex))])
    elif number == 4:
        # I @ Y
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), qudit_operation_matrix(sigma2, dim)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2, dtype=complex))])
    elif number == 5:
        # Y @ Y
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma2, dim), qudit_operation_matrix(sigma2, dim)), np.sqrt(1-p)*np.sqrt(np.identity(dim**2, dtype=complex)-np.kron(qudit_operation_matrix(sigma1, dim), qudit_operation_matrix(sigma1, dim)))])
    if number == 6:
        # Z @ I
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma3, dim), np.identity(dim, dtype=complex)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2))])
    elif number == 7:
        # I @ Z
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), qudit_operation_matrix(sigma3, dim)), np.sqrt(1-p)*np.asmatrix(np.identity(dim**2))])
    elif number == 8:
        # Z @ Z
        density_matrix = qudit_operator_sum(density_matrix, [np.sqrt(p)*np.kron(qudit_operation_matrix(sigma3, dim), qudit_operation_matrix(sigma3, dim)), np.sqrt(1-p)*np.sqrt(np.identity(dim**2, dtype=complex) - np.kron(qudit_operation_matrix(sigma3, dim), qudit_operation_matrix(sigma3, dim)))])
    elif number == 9:
        # general_bit_flip @ I
        E1 = np.sqrt(p)*np.kron(X(dim), np.identity(dim, dtype=complex))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 10:
        # I @ general_bit_flip
        E1 = np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), X(dim))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 11:
        # general_bit_flip @ general_bit_flip
        E1 = np.sqrt(p)*np.kron(X(dim), X(dim))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 12:
        # general_phase_flip @ I
        E1 = np.sqrt(p)*np.kron(Z(dim), np.identity(dim, dtype=complex))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 13:
        # I @ general_phase_flip
        E1 = np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), Z(dim))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 14:
        # general_phase_flip @ general_phase_flip
        E1 = np.sqrt(p)*np.kron(Z(dim), Z(dim))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 15:
        # depolarization @ I
        operatorlist = []
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for n in range(dim):
            for m in range(dim):
                operatorlist.append(np.sqrt(p)/dim*(np.kron(X(dim)**n, np.identity(dim)))*(np.kron(Z(dim)**m, np.identity(dim, dtype=complex))))
                sum += operatorlist[-1].H*operatorlist[-1]
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 16:
        # I @ depolarization
        operatorlist = []
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for n in range(dim):
            for m in range(dim):
                operatorlist.append(np.sqrt(p)/dim*(np.kron(np.identity(dim), X(dim)**n))*(np.kron(np.identity(dim, dtype=complex), Z(dim)**m)))
                sum += operatorlist[-1].H*operatorlist[-1]
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 17:
        # depolarization @ depolarization
        operatorlist = []
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for n in range(dim):
            for m in range(dim):
                operatorlist.append(np.sqrt(p)/dim*(np.kron(X(dim)**n, X(dim)**n))*(np.kron(Z(dim)**m, Z(dim)**m)))
                sum += operatorlist[-1].H*operatorlist[-1]
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 18:
        # single_Weyl_channel @ I
        random.seed(p)
        n = random.randint(0, dim)
        m = random.randint(0, dim)
        E1 = np.sqrt(p)*np.kron((X(dim)**n)*(Z(dim)**m), np.identity(dim, dtype=complex))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 19:
        # I @ single_Weyl_channel
        random.seed(p)
        n = random.randint(0, dim)
        m = random.randint(0, dim)
        E1 = np.sqrt(p)*np.kron(np.identity(dim, dtype=complex), (X(dim)**n)*(Z(dim)**m))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 20:
        # single_Weyl_channel @ single_Weyl_channel
        random.seed(p)
        n = random.randint(0, dim)
        m = random.randint(0, dim)
        E1 = np.sqrt(p)*np.kron((X(dim)**n)*(Z(dim)**m), (X(dim)**n)*(Z(dim)**m))
        E0 = np.sqrt(np.identity(dim**2, dtype=complex)-(E1.H)*E1)
        density_matrix = qudit_operator_sum(density_matrix, [E1, E0])
    elif number == 21:
        # amplitude_damping @ I
        random.seed(p)
        gammalist = []
        for i in range(round(0.5*dim*(dim-1))):
            gammalist.append(p*random.random())
        operatorlist = []
        j = 0
        k = 1
        minus = (dim-2-j)
        for i in range(round(0.5*dim*(dim-1))):
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
            operatorlist.append(np.kron(np.sqrt(gammalist[i])*basisvec1.T*basisvec2, np.identity(dim, dtype=complex)))
            k += 1
            sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 22:
        # amplitude_damping @ I
        random.seed(p)
        gammalist = []
        for i in range(round(0.5*dim*(dim-1))):
            gammalist.append(p*random.random())
        operatorlist = []
        j = 0
        k = 1
        minus = (dim-2-j)
        for i in range(round(0.5*dim*(dim-1))):
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
            operatorlist.append(np.kron(np.identity(dim, dtype=complex), np.sqrt(gammalist[i])*basisvec1.T*basisvec2))
            k += 1
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 23:
        # amplitude_damping @ I
        random.seed(p)
        gammalist = []
        for i in range(round(0.5*dim*(dim-1))):
            gammalist.append(p*random.random())
        operatorlist = []
        j = 0
        k = 1
        minus = (dim-2-j)
        for i in range(round(0.5*dim*(dim-1))):
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
            operatorlist.append(np.kron(np.sqrt(gammalist[i])*basisvec1.T*basisvec2, np.sqrt(gammalist[i])*basisvec1.T*basisvec2))
            k += 1
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 24:
        # general_phase_damping @ I
        operatorlist = []
        for n in range(dim):
            operatorlist.append(np.kron(np.sqrt(scipy.special.binom(dim-1, n)*((0.5*(1-p))**n)*((0.5*(1+p))**(dim-1-n)))*(Z(dim)**n), np.identity(dim, dtype=complex)))
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 25:
        # I @ general_phase_damping
        operatorlist = []
        for n in range(dim):
            operatorlist.append(np.kron(np.identity(dim, dtype=complex), np.sqrt(scipy.special.binom(dim-1, n)*((0.5*(1-p))**n)*((0.5*(1+p))**(dim-1-n)))*(Z(dim)**n)))
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 26:
        # general_phase_damping @ general_phase_damping
        operatorlist = []
        for n in range(dim):
            operatorlist.append(np.kron(
            np.sqrt(scipy.special.binom(dim-1, n)*((0.5*(1-p))**n)*((0.5*(1+p))**(dim-1-n)))*(Z(dim)**n), np.sqrt(scipy.special.binom(dim-1, n)*((0.5*(1-p))**n)*((0.5*(1+p))**(dim-1-n)))*(Z(dim)**n)))
        sum = np.asmatrix(np.zeros((dim**2,dim**2), dtype=complex))
        for op in operatorlist:
            sum += op.H*op
        operatorlist.append(np.sqrt(np.identity(dim**2, dtype=complex)-sum))
        density_matrix = qudit_operator_sum(density_matrix, operatorlist)
    elif number == 27:
        # physical_phase_damping @ I
        for i in range(dim**2):
            for j in range(dim**2):
                i1 = i%dim
                i2 = int(i/dim)
                j1 = j%dim
                j2 = int(j/dim)
                density_matrix[i,j] = density_matrix[i,j]*((np.exp(-p))**((i1-j1)**2))
    elif number == 28:
        # physical_phase_damping @ I
        for i in range(dim**2):
            for j in range(dim**2):
                i1 = i%dim
                i2 = int(i/dim)
                j1 = j%dim
                j2 = int(j/dim)
                density_matrix[i,j] = density_matrix[i,j]*((np.exp(-p))**((i2-j2)**2))
    elif number == 29:
        # physical_phase_damping @ I
        for i in range(dim**2):
            for j in range(dim**2):
                i1 = i%dim
                i2 = int(i/dim)
                j1 = j%dim
                j2 = int(j/dim)
                density_matrix[i,j] = density_matrix[i,j]*((np.exp(-p))**((i1-j1)**2)+(i2-j2)**2)
    return density_matrix

def qudit_operation_matrix(mat, dim):
    gate_matrix = np.asmatrix(np.identity(dim, dtype=complex))
    gate_matrix[0,0] = mat[0,0]
    gate_matrix[0,1] = mat[0,1]
    gate_matrix[1,0] = mat[1,0]
    gate_matrix[1,1] = mat[1,1]
    return gate_matrix

def two_qubit_to_qudit_gate(mat, dim):
    gate_matrix = np.asmatrix(np.identity(dim**2, dtype=complex))
    gate_matrix[0,0] = mat[0,0]
    gate_matrix[0,1] = mat[0,1]
    gate_matrix[1,0] = mat[1,0]
    gate_matrix[1,1] = mat[1,1]
    gate_matrix[dim,dim] = mat[2,2]
    gate_matrix[dim,dim+1] = mat[2,3]
    gate_matrix[dim+1,dim] = mat[3,2]
    gate_matrix[dim+1,dim+1] = mat[3,3]
    gate_matrix[0,dim] = mat[0,2]
    gate_matrix[0,dim+1] = mat[0,3]
    gate_matrix[dim,0] = mat[2,0]
    gate_matrix[dim+1,0] = mat[3,0]
    gate_matrix[1,dim] = mat[1,2]
    gate_matrix[1,dim+1] = mat[1,3]
    gate_matrix[dim,1] = mat[2,1]
    gate_matrix[dim+1,1] = mat[3,1]
    return gate_matrix
