
import numpy as np
from qiskit.quantum_info import Statevector
import qiskit.quantum_info as qi
from qiskit.circuit.quantumcircuit import Parameter, QuantumCircuit, Gate


class Orth4(Gate):

    def __init__(self, parameters):
        super().__init__('U', 2, parameters)
        
    def _u_gate(self):
        qc = QuantumCircuit(2)
        qc.unitary(self.to_matrix(), [0, 1])
        self.definition = qc

    
    def compose_orth4(self,params):

        CNOT2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        S = np.array([[1,0],[0,1.j]])
        Sdag = np.array([[1,0],[0,-1j]])
        
        I = np.eye(2)

        def ry(theta):
            a = np.cos(theta/2)
            b = np.sin(theta/2)
            return np.array([[a,-b],[b,a]])

        def rz(theta):
            a = np.exp(-1.j*theta/2)
            b = np.exp(1.j*theta/2)
            return np.array([[a,0],[0,b]])
        
        def central_su2s(x):
            return np.kron(rz(x[5])@ry(x[4])@rz(x[3]),I) @ np.kron(I, rz(x[2])@ry(x[1])@rz(x[0]))
        
        return np.kron(I, Sdag) @ np.kron(Sdag, I) @ np.kron(ry(-np.pi/2), I) @ CNOT2 @ central_su2s(params) @ CNOT2 @ np.kron(ry(np.pi/2), I) @ np.kron(S, I) @ np.kron(I, S)
    
    def to_matrix(self):
        thetas = []
        for i in range(0,6):
            thetas.append(float(self.params[i]))
        return np.array(self.compose_orth4(thetas))

def orth4_layer(counter, qc, ent_map, M):
    params = []
    for c,t in ent_map:
        for _ in range(0,6):
            counter = counter + 1
            params.append(Parameter(parameterBitString(M,counter)))
        qc.append(Orth4(params), [c, t])
        qc.barrier()
        params = []
    return counter

class SU4(Gate):

    def __init__(self, parameters):
        super().__init__('U', 2, parameters)
        
    def _u_gate(self):
        qc = QuantumCircuit(2)
        qc.unitary(self.to_matrix(), [0, 1])
        self.definition = qc

    
    def compose_su4(self,params):

        CNOT2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        CNOT1 = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

        I = np.eye(2)
        def ry(theta):
            a = np.cos(theta/2)
            b = np.sin(theta/2)
            return np.array([[a,-b],[b,a]])

        def rz(theta):
            a = np.exp(-1j*theta/2)
            b = np.exp(1J*theta/2)
            return np.array([[a,0],[0,b]])

        def swap_block(x):
            return CNOT2 @ np.kron(ry(x[2]),I) @ CNOT1 @ np.kron(ry(x[1]), I) @ np.kron(I, rz(x[0])) @ CNOT2
        
        def general_su2(x):
            return np.kron(rz(x[5])@ry(x[4])@rz(x[3]),I) @ np.kron(I, rz(x[2])@ry(x[1])@rz(x[0]))
        
        return general_su2(params[9:15:]) @ swap_block(params[6:9:]) @ general_su2(params[:6:])
    def to_matrix(self):
        thetas = []
        for i in range(0,15):
            thetas.append(float(self.params[i]))
        return np.array(self.compose_su4(thetas))

def su4_layer(counter, qc, ent_map, M):
    params = []
    for c,t in ent_map:
        for _ in range(0,15):
            counter = counter + 1
            params.append(Parameter(parameterBitString(M,counter)))
        qc.append(SU4(params), [c, t])
        qc.barrier()
        params = []
    return counter


class real_SU4(Gate):

    def __init__(self, parameters):
        super().__init__('U', 2, parameters)
        
    def _u_gate(self):
        qc = QuantumCircuit(2)
        qc.unitary(self.to_matrix(), [0, 1])
        self.definition = qc

    
    def compose_real_su4(self,params):

        CNOT2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        CNOT1 = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

        I = np.eye(2)
        def Ry(theta):
            a = np.cos(theta/2)
            b = np.sin(theta/2)
            return np.array([[a,-b],[b,a]])
        
        def su_2(x):
            return np.kron(Ry(x[1]),I) @ np.kron(I, Ry(x[0]))
        
        def swap_block(x):
            return CNOT2 @ np.kron(Ry(x[1]),I) @ CNOT1 @ np.kron(Ry(x[0]), I) @ CNOT2
        
        return su_2([params[4],params[5]]) @ swap_block([params[2], params[3]]) @ su_2([params[0], params[1]])

        #return np.kron(Ry(params[5]),I) @ np.kron(I, Ry(params[4]))@ CNOT2 @ np.kron(Ry(params[3]),I)@CNOT1@np.kron(Ry(params[2]),I) @ CNOT2 @ np.kron(Ry(params[1]),I) @ np.kron(I, Ry(params[0]))

    def to_matrix(self):
        thetas = []
        for i in range(0,6):
            thetas.append(float(self.params[i]))
        return np.array(self.compose_real_su4(thetas))

def real_su4_layer(counter, qc, ent_map, M):
    params = []
    for c,t in ent_map:
        for i in range(0,6):
            counter = counter + 1
            params.append(Parameter(parameterBitString(M, counter)))
        qc.append(real_SU4(params), [c, t])
        qc.barrier()
        params = []
    return counter



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
        qc.barrier()
    return counter
     
        
def parameterBitString(M,i):
    return format(i,'b').zfill(M)


def general_SU4(counter, qc, q0,q1,M):
    def add_single_SU2 (counter, qc,q, M):
        counter = counter + 1
        qc.rz(Parameter(parameterBitString(M,counter)),q)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q)
        counter = counter + 1
        qc.rz(Parameter(parameterBitString(M,counter)),q)
        return counter


    def add_N_block(counter, qc, q0,q1, M):
        qc.cx(q1, q0)
        counter = counter + 1
        qc.rz(Parameter(parameterBitString(M,counter)),q0)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q1)
        qc.cx(q0, q1)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q1)
        qc.cx(q1, q0)

        return counter

    counter = add_single_SU2(counter, qc, q0, M)
    counter = add_single_SU2(counter, qc, q1, M)
    qc.barrier()
    counter = add_N_block(counter, qc, q0,q1, M)
    qc.barrier()
    counter = add_single_SU2(counter, qc, q0, M)
    counter = add_single_SU2(counter, qc, q1, M)
    qc.barrier()
    return counter

def general_real_SU4(counter, qc, q0,q1,M):
    def add_single_SU2 (counter, qc,q, M):
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q)
        return counter
    
    def add_N_block(counter, qc, q0,q1, M):
        qc.cx(q1, q0)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q1)
        qc.cx(q0, q1)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q1)
        qc.cx(q1, q0)
        return counter

    counter = add_single_SU2(counter, qc, q0, M)
    counter = add_single_SU2(counter, qc, q1, M)
    qc.barrier()
    counter = add_N_block(counter, qc, q0,q1, M)
    qc.barrier()
    counter = add_single_SU2(counter, qc, q0, M)
    counter = add_single_SU2(counter, qc, q1, M)
    qc.barrier()
    return counter

def general_orth4(counter, qc, q0,q1,M):
    def add_single_SU2 (counter, qc,q, M):
        counter = counter + 1
        qc.rz(Parameter(parameterBitString(M,counter)),q)
        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q)
        counter = counter + 1
        qc.rz(Parameter(parameterBitString(M,counter)),q)
        return counter
    
    qc.s(q0)
    qc.s(q1)
    qc.ry(np.pi/2,q1)
    qc.barrier()
    qc.cx(q1,q0)
    counter = add_single_SU2(counter, qc, q0, M)
    counter = add_single_SU2(counter, qc, q1, M )
    qc.cx(q1,q0)
    qc.barrier()
    qc.id(q0)
    qc.ry(-np.pi/2,q1)
    qc.sdg(q0)
    qc.sdg(q1)
    return counter





params = np.random.randint(size=6, low= -np.pi, high=np.pi) 


qc = QuantumCircuit(2)
counter = -1
ent = [[0,1]]
M = 10

general_orth4(counter=counter, qc=qc, q0=0,q1=1,M=10)
print(qc.draw())
print(params)
qc.assign_parameters(params, inplace=True)

mat_qiskit = np.array(qi.Operator(qc))
state = Statevector(qc)
psi = np.array(state)



qc2 = QuantumCircuit(2)
counter = -1
ent = [[0,1]]
M = 10
mat_fabio = Orth4(params).to_matrix()

print(np.round(mat_qiskit,4))
print("\n\n")
print(np.round(mat_fabio,4))
print("\n\n")
"""

params = np.random.randint(size=15, low= -np.pi, high=np.pi)
qc = QuantumCircuit(2)
counter = -1
ent = [[0,1]]
M = 10

general_SU4(counter=counter, qc=qc, q0=0,q1=1,M=10)
print(params)
qc.assign_parameters(params, inplace=True)

mat_qiskit = np.array(qi.Operator(qc))
state = Statevector(qc)
psi = np.array(state)

#print(state)
#print(np.real(psi.T@psi))

print("\n\n")

qc2 = QuantumCircuit(2)
counter = -1
ent = [[0,1]]
M = 10
counter = su4_layer(counter=counter,qc = qc2, ent_map=ent, M=M )
mat_fabio = SU4(params).to_matrix()
qc2.assign_parameters(params, inplace=True)
state = Statevector(qc2)
psi = np.array(state)
print(mat_qiskit)
print("\n\n")
print(mat_fabio)
print("\n\n")

print(np.linalg.det(mat_fabio))
print(np.linalg.det(mat_qiskit))


counter = -1


depth = 1
ENTANGLER_MAP = [[1,0]]#,[1,2]],[2,5],[5,8],[8,7],[7,6],[6,3]]


Nx = 1
Ny = 2
qc = QuantumCircuit(Nx*Ny)
M = 10


qc.barrier()
qc.x(1)
print("BEFORE LAYER", counter)
counter = real_su4_layer(counter, qc, ENTANGLER_MAP, M)
print("AFTER LAYER", counter)
print(qc.draw())
print(qc.num_parameters)
print(qc.num_nonlocal_gates())
print("AFTER PARAMETERS ASSIGNEMT")

params = [np.pi]*qc.num_parameters
qc.assign_parameters(params, inplace=True)
state = Statevector(qc)
psi = np.array(state)
print(state)
print(qc.draw())
print(params)
print(np.real(psi.T@psi))
"""
    






