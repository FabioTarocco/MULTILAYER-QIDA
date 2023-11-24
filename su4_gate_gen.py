
import numpy as np
from qiskit.quantum_info import Statevector
from qiskit.circuit.quantumcircuit import Parameter, QuantumCircuit



def general_SU4(counter, qc, q0,q1,M):
    def add_single_SU2 (counter, qc,q, M):

        #counter = counter + 1
        #qc.rz(Parameter(parameterBitString(M,counter)),q)

        counter = counter + 1
        qc.ry(Parameter(parameterBitString(M,counter)),q)

        #counter = counter + 1
        #qc.rz(Parameter(parameterBitString(M,counter)),q)
        return counter


    def add_N_block(counter, qc, q0,q1, M):
        qc.cx(q1, q0)

        #counter = counter + 1
        #qc.rz(Parameter(parameterBitString(M,counter)),q0)
        
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


def parameterBitString(M,i):
    return format(i,'b').zfill(M)
    
def extract_N_block(paramslist,ent,depth):
    ct = len(ent)
    n = len(set(np.array(ent).flatten()))
    paramslist = paramslist[n::]
    sliced = []
    for _ in range(depth):
        for _ in range(ct):
            sliced.append(paramslist[:6:])
            paramslist = paramslist[6::]
        paramslist = paramslist[n::]
    
    relevant = []
    for s in sliced:
        relevant.append(s[2:4:])
    relevant = list(np.concatenate(relevant).flat)
    x_ = relevant[::2]
    y_ = relevant[1::2]
    return x_,y_
    

count = -1

depth = 1
ENTANGLER_MAP = [[0,1],[1,2]]#[2,5],[5,8],[8,7],[7,6],[6,3]]


Nx = 4
Ny = 1
qc = QuantumCircuit(Nx*Ny)
M = 8

for d in range(0,depth):
    for i in range(0,Nx*Ny):
        if i in(set(np.array(ENTANGLER_MAP).flatten())):
            count = count + 1
            qc.rz(Parameter(parameterBitString(M,count)),i)
    
    qc.barrier()

    for a,b in ENTANGLER_MAP:

        count = count + 1
        count = general_SU4(count, qc, a, b,M)
        
for i in range(0,Nx*Ny):
    count = count + 1
    qc.rz(Parameter(parameterBitString(M,count)),i)

print(qc.num_parameters)
print(qc.num_nonlocal_gates())
p = np.random.rand(qc.num_parameters)
print(qc.assign_parameters(p, inplace=False).draw())
print(extract_N_block(p,ENTANGLER_MAP, depth))








