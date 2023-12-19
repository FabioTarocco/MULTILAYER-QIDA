import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from qiskit.quantum_info.operators import SparsePauliOp
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit

def create_single_string(op,pair,N):
    temp = list("I"*N)
    temp[pair[1]] = op
    temp[pair[0]] = op
    return "".join(temp)

def heisenberg_model(p, ops, J, h, N):
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

def create_pairs(dim, model_lattice="grid"):
    if len(dim)==1:
        n_sites = dim[0]
    else:
        nrow =  dim[0]
        ncol = dim[1]
    pairs = []
    if model_lattice =="grid":
        lattice = nx.grid_2d_graph(nrow, ncol)
        adj_matrix = nx.to_numpy_array(lattice)
        for i in range (0,ncol*nrow):
            for j in range(i,ncol*nrow):
                adj_matrix[j,i]=0

        for i in range(0,nrow*ncol):
            for j in range(0,nrow*ncol):
                if adj_matrix[i, j] == 1:
                    pairs.append((i, j))
        
    if model_lattice == "triangle":
        sites = np.arange(.5*n_sites*(n_sites+1), dtype=int)[::-1]
        lists = []
        for i in range(1,n_sites+1):
            lists.append(list(sites[:i:])[::-1])
            sites = sites[i:]
        lists = lists[::-1]
        for i in range(1,n_sites):
            for j in range(len(lists[i])):
                pairs.append((lists[i-1][j],lists[i][j]))
                pairs.append((lists[i-1][j+1],lists[i][j]))
                pairs.append((lists[i-1][j], lists[i-1][j+1]))

    if model_lattice =="triangular_lattice":
        sites = np.arange(nrow*ncol, dtype=int)[::-1]
        lists = []
        for i in range(0,ncol):
            lists.append(list(sites[:nrow:])[::-1])
            sites = sites[nrow:]
        lists = lists[::-1]
        print(lists)
        for i in range(1,nrow):
            for j in range(len(lists[i])):
                pairs.append((lists[i-1][j],lists[i][j]))
                if j != nrow-1:
                    pairs.append((lists[i-1][j+1],lists[i][j]))
                    pairs.append((lists[i-1][j], lists[i-1][j+1]))
                if i == nrow-1 and j != len(lists[i])-1:
                    pairs.append((lists[i][j], lists[i][j+1]))
        
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

couple = create_pairs([Nrow, Ncol], model_lattice = "grid")
h_dict = heisenberg_model(couple, interaction, interaction_term, mag_field, Nrow*Ncol)
h_tot = ham_as_matrix(h_dict)
print(couple)
print(h_dict)
print(h_tot)
#-------NEEL STATE FOR REFERENCE------
print(neel_ref_en(Nrow*Ncol,h_tot))
"""
couple = create_pairs([3,3], model_lattice = "triangular_lattice")
print(couple)
exit(0)
g = nx.from_edgelist(couple)
plt.figure()
random_pos = nx.random_layout(g, seed=19)
pos = nx.spring_layout(g, pos=random_pos)
nx.draw(g, with_labels = True, pos=pos)
plt.show()
print(couple)
#qc = QuantumCircuit(9)
#qc.cx(0,1)
#print(qc.num_nonlocal_gates())
