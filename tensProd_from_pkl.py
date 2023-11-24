import numpy as np
import pandas as pd
import pickle
from pprint import pprint

ket0 = np.array([[1],[0]])
ket1 = np.array([[0],[1]])
ket0_t = (1,0)
ket1_t = (0,1)

print(type(ket0_t))

test = dict()
test = {  '|000011100>': 0.76, 
        '|000001101>': 0.12,
        '|000101011>': 0.12}


def gs_to_vector(gs,N):
    wf = 0.0
    for v,coeff in gs.items():
        res = coeff
        for i in v[1:-1]:
            res = np.kron(res,ket0 if i =='0' else ket1)
        wf += res
    assert len(wf) == 2**N, "\nERROR\n\tLength of the vector is {}. For {}-qubit the length is {}.".format(len(wf), N, 2**N)
    return wf
    
def gs_to_tensProd(gs):
    coeffs = list()
    wf = list()
    for v,coeff in gs.items():
        coeffs.append(coeff)
        res = list()
        for i in v[1:-1]:
            res.append(ket0_t if i =='0' else ket1_t)
        wf.append(res)
    return wf, coeffs
    
N = 6
filet = pd.read_pickle("./pickle/gs2x3_dict.pkl")[0]
print(filet)

gs, coeffs = gs_to_tensProd(filet)


try:
    with open("./gs2x3.pkl", 'rb') as file:
        data = pickle.load(file)
except FileNotFoundError:
    data = []

with open("./gs2x3.pkl", 'wb') as file:
    pickle.dump(data, file)

new_entry = {'gs':gs, 'coeffs':coeffs, 'gs_dict':filet}
data.append(new_entry)

with open("./gs2x3.pkl", 'wb') as file:
        pickle.dump(data, file)



