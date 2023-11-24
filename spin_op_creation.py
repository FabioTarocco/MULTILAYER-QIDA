import numpy as np
import networkx as nx
from qiskit.quantum_info.operators import SparsePauliOp
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit

def create_single_string(op,pair,N):
    temp = list("I"*N)
    temp[pair[1]] = op
    temp[pair[0]] = op
    return "".join(temp)

def create_ham_dict(p, ops, J, h, N):
    ham_dict = dict()
    for a,b in p:
        for op_i in range(len(ops)):
            ham_dict[create_single_string(ops[op_i],(a,b),N)] = J[op_i]
    if h != 0.0:
        for n in range(0,N):
            ham_dict[create_single_string("Z",(n,n),N)] = h
    return ham_dict

def ham_as_matrix(H):
    H_tot = 0
    for SpinString in H.keys():
        H_tot = H_tot + SparsePauliOp.from_list([(SpinString,H[SpinString])]).to_matrix()
    return H_tot

def create_pairs(nrow, ncol):
    pairs = []
    grid = nx.grid_2d_graph(nrow, ncol)
    adj_matrix = nx.to_numpy_array(grid)
    for i in range (0,ncol*nrow):
        for j in range(i,ncol*nrow):
            adj_matrix[j,i]=0

    for i in range(0,nrow*ncol):
        for j in range(0,nrow*ncol):
            if adj_matrix[i, j] == 1:
                pairs.append((i, j))
    return pairs


def neel_ref_en(N,O):
    circuit = QuantumCircuit(N)
    for i in range(0,N,2):
        circuit.x(i)
    state = Statevector(circuit)
    s = np.array(Statevector(circuit))
    return np.real(s.conjugate()@O@s)
"""

#-----------------
#LATTICE SETTINGS
#-----------------
Nrow = 3
Ncol = 4

#--------------------------------------
#-------Hamiltonian Settings-----------
#--------------------------------------
interaction = ["X","Y","Z"]
interaction_term = [0.25, 0.25, 0.25]
mag_field = 0.0

#--------------------------------------
#----------------MAIN------------------
#--------------------------------------

couple = create_pairs(Nrow, Ncol)
h_dict = create_ham_dict(couple, interaction, interaction_term, mag_field, Nrow*Ncol)
h_tot = ham_as_matrix(h_dict)
print(couple)
print(h_dict)
print(h_tot)
#-------NEEL STATE FOR REFERENCE------
print(neel_ref_en(Nrow*Ncol,h_tot))
"""
qc = QuantumCircuit(9)
qc.cx(0,1)
print(qc.num_nonlocal_gates())
