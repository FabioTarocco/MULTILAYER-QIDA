
import numpy as np
from qiskit.quantum_info import Statevector
from qiskit.circuit.quantumcircuit import Parameter, QuantumCircuit, Gate

class particle_conserving_singleU(Gate):
    def __init__(self, parameters):
        super().__init__('U', 2, parameters)
        
    def _u_gate(self):
        qc = QuantumCircuit(2)
        qc.unitary(self.to_matrix(), [0, 1])
        self.definition = qc
        
    def to_matrix(self):
    
        theta_1 = float(self.params[0])
        theta_2 = float(self.params[1])
        
        a = np.cos(theta_1)
        b = np.exp(1j*theta_2)*np.sin(theta_1)
        c = np.exp(1j*(-theta_2))*np.sin(theta_1)
        d = -np.cos(theta_1)
        
        return np.array([[1, 0, 0, 0], 
                         [0, a, b, 0], 
                         [0, c, d, 0],
                         [0, 0, 0, 1]])


def parity_conserving_layer(counter, qc, ent_map, M):
    params = []
    for c,t in ent_map:
        counter = counter + 1
        theta_1 = Parameter(parameterBitString(M,counter))
        counter = counter + 1
        theta_2 = Parameter(parameterBitString(M,counter))
        params =  [theta_1, theta_2]
        qc.append(particle_conserving_singleU(params), [c, t])
        print(counter)
        qc.barrier()
    return counter
    
    
    
class parity_conserving_singleU(Gate):
    def __init__(self, parameters):
        super().__init__('U', 2, parameters)
        
        
    def _u_gate(self):
        qc = QuantumCircuit(2)
        qc.unitary(self.to_matrix(), [0, 1])
        self.definition = qc
        
    def to_matrix(self):
    
        theta_1 = float(self.params[0])
        theta_2 = float(self.params[1])
        
        a = np.cos(theta_1)
        b = np.exp(1.j*theta_2)*np.sin(theta_1)
        c = np.exp(1.j*(-theta_2))*np.sin(theta_1)
        d = -np.cos(theta_1)
        
        return np.array([[a, 0, 0, b], 
                         [0, 1, 0, 0], 
                         [0, 0, 1, 0],
                         [c, 0, 0, d]])


def particle_conserving_layer(counter, qc, ent_map, M):
    params = []
    for c,t in ent_map:
        counter = counter + 1
        theta_1 = Parameter(parameterBitString(M,counter))
        counter = counter + 1
        theta_2 = Parameter(parameterBitString(M,counter))
        params =  [theta_1, theta_2]
        qc.append(particle_conserving_singleU(params), [c, t])
        print(counter)
        qc.barrier()
    return counter
     
        
def parameterBitString(M,i):
    return format(i,'b').zfill(M)
    
    
counter = -1


depth = 1
ENTANGLER_MAP = [[0,1],[1,2],[2,5],[5,8],[8,7],[7,6],[6,3]]


Nx = 3
Ny = 3
qc = QuantumCircuit(Nx*Ny)
M = 10


qc.barrier()
print("BEFORE LAYER", counter)
counter = parity_conserving_layer(counter, qc, ENTANGLER_MAP, M)
print("AFTER LAYER", counter)
    

print(qc.draw())
print(qc.num_parameters)
print(qc.num_nonlocal_gates())

params = np.random.uniform(size = qc.num_parameters, high = np.pi, low = -np.pi)
qc.assign_parameters(params, inplace=True)
state = Statevector(qc)
psi = np.array(state)
print(qc.draw())
print(params)
print(np.real(psi.T@psi))






