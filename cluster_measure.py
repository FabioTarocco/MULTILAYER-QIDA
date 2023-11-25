import pandas as pd
import numpy as np
from scipy import linalg as LA
from scipy.optimize import minimize
from scipy.optimize import BFGS


# 1-QUBIT GATES
I = np.eye(2)
Z = np.array([[1., 0.],[0. ,-1.]])
X = np.array([[0., 1.],[1. ,0.]])
Y = np.array([[0., -1.j],[1.j ,0.]])

# 2-QUBITs GATES
I4 = np.eye(4)
CNOT1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
CNOT2 = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]])

CZ = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])

CY1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,-1j],[0,0,1j,0]])
CY2 = np.array([[1,0,0,0],[0,0,0,-1j],[0,0,1,0],[0,1j,0,0]])

CH1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1/np.sqrt(2),1/np.sqrt(2)],[0,0,1/np.sqrt(2),-1/np.sqrt(2)]])
CH2 = np.array([[1,0,0,0],[0,1/np.sqrt(2),0,1/np.sqrt(2)],[0,0,1,0],[0,1/np.sqrt(2),0,-1/np.sqrt(2)]])

SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
iSWAP = np.array([[1,0,0,0],[0,0,1j,0],[0,1j,0,0],[0,0,0,1]])

M_LIST = [I4, 
          CNOT1, CNOT2,
          CZ,
          CY1,
          CY2,
          CH1,
          CH2,
          SWAP, iSWAP]

def ry(theta):
    a = np.cos(theta/2)
    b = np.sin(theta/2)
    return np.array([[a,-b],[b,a]])

def rz(theta):
    a = np.exp(1j*theta/2)
    b = np.exp(-1J*theta/2)
    return np.array([[a,0],[0,b]])
   
def compose_N_block(t1,t2):
    return CNOT2 @ np.kron(I, ry(t1)) @ CNOT1 @ np.kron(I, ry(t2)) @ CNOT2

def compose_su4(params):
    return np.kron(ry(params[0]),I) @ np.kron(I, ry(params[1])) @ compose_N_block(params[2],params[3]) @ np.kron(ry(params[4]),I) @ np.kron(I, ry(params[5]))

def measure(a):
    all_diff = []
    for b in M_LIST:
        dist = LA.norm(a - b, 'fro')
        all_diff.append(dist)
    return all_diff

def distance (x):
    r = compose_N_block(x[0],x[1])
    return np.real(LA.norm(I4 - r))

def distance_su4(x):
    r = compose_su4(x)
    return np.real(LA.norm(CH1 - r))


convergence = []
parameters = []

def optimize():
    print("Starting Optimization!")
    def callback(x_params):
        convergence.append(distance_su4(x_params))
        print(convergence[-1])
        parameters.append(x_params)

    res = minimize(fun = distance_su4, 
                   x0 = np.random.uniform(size = 6, low=-np.pi, high=np.pi),
                   callback=callback,
                   tol=1e-6,
                   method='BFGS'
                   )
    parameters.append(res.x)
    return res.x, distance_su4(res.x)

par, dist = optimize()
print(par, dist)

print(np.round(compose_su4(par),3))

#print("\n")
#print(np.round(compose_N_block(np.mod(np.pi,2*np.pi),np.mod(0, 2*np.pi))))
#print("\n")
#halfs = [np.pi*k/2 for k in range(0,5)]

#print(halfs)
#or i in halfs:
#    print(np.round(ry(np.mod(i,2*np.pi)),4))