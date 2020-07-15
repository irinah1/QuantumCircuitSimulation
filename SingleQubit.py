import numpy as np
import cmath
from constants import *

class Qubit:
    """
    input:
    
    r = [rx, ry, rz]: Blochvector entries: rho = 0.5*(I + r*sigma_vector)
    """
    def __init__(self, vec):
        scale = 1.
        if (vec[0]**2+vec[1]**2+vec[2]**2 > 1):
            print("Blochvector r has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
            scale = 1/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
        if type(vec)==list:
            vec = np.array(vec)
        self.bloch_vector = scale*vec
        self.matrix = 0.5*(I + self.bloch_vector[0]*sigma1 + self.bloch_vector[1]*sigma2 + self.bloch_vector[2]*sigma3)
    
    #def propagate(self, t):
    #    # t: time (in seconds)
    #    self.matrix = np.matrix([[1+(self.matrix[0,0]-1)*np.exp(-self.G1*t), self.matrix[0,1]*np.exp(1j*self.dw*t)*np.exp(-self.G2*t)],
    #              [self.matrix[1,0]*np.exp(-1j*self.dw*t)*np.exp(-self.G2*t), self.matrix[1,1]*np.exp(-self.G1*t)]])
    
    def get_matrix(self):
        return self.matrix
    
    def _set_Bloch_vector(self, vec):
        if (vec[0]**2+vec[1]**2+vec[2]**2 > 1):
            print("Blochvector r has to fullfill rx^2 + ry^2 + rz^2 <= 1.")
        if type(vec)==list:
            vec = np.array(vec)
        self.bloch_vector = vec
        self.matrix = 0.5*(I + vec[0]*sigma1 + vec[1]*sigma2 + vec[2]*sigma3)
        
    def get_Bloch_vector(self):
        return self.bloch_vector[0], self.bloch_vector[1], self.bloch_vector[2]

def plot_qubit(qubit):
    if type(qubit) != Qubit:
        print("Input of one_qubit_gate_operation must be of type Qubit.")
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    import numpy as np

    a,b,c = qubit.get_Bloch_vector()

    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111, projection='3d')

    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = 1 * np.outer(np.cos(u), np.sin(v))
    y = 1 * np.outer(np.sin(u), np.sin(v))
    z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z, color='gray', alpha=0.1)
    ax.quiver(0, 0, 0, a, b, c, length=1, normalize=True)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    plt.show()
