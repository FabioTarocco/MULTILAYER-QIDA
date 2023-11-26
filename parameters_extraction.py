
import numpy as np
from qiskit.quantum_info import Statevector
from qiskit.circuit.quantumcircuit import Parameter, QuantumCircuit
import pandas as pd
import sys
from pprint import pprint
if len(sys.argv) ==1:
    flag_plot = False
else:
    flag_plot = True
    t = sys.argv[1]
    bins = sys.argv[2]
    prec = sys.argv[3]
    if len(sys.argv) == 5:
        s = float(sys.argv[4])
    else:
        s = 100

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
    
def extract_N_block(paramslist,ent,depth,offset):
    ct = len(ent)
    #print(ct)
    n = len(set(np.array(ent).flatten()))
    #print(n)
    paramslist = paramslist[n::]
    sliced = []
    for i in range(depth):
        for _ in range(ct):
            if i > (offset-1):
                #print("Depth {}".format(i+1))
                sliced.append(paramslist[:6:])
            paramslist = paramslist[6::]
        paramslist = paramslist[n::]
    print(sliced)
    relevant = []
    for s in sliced:
        relevant.append(s[2:4:])
    if relevant == []:
        return [],  []
    else:
        relevant = list(np.concatenate(relevant).flat)
        x_ = relevant[::2]
        y_ = relevant[1::2]
        return x_,y_
    

count = -1

depth = 1
ladderl = 6
qidal = 1
fullyl=1

ENTANGLER_MAP = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]#[2,5],[5,8],[8,7],[7,6],[6,3]]
ENTANGLER_MAP_F = [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7],[0,8],
                   [1,2],[1,3],[1,4],[1,5],[1,6],[1,7],[1,8],
                   [2,3],[2,4],[2,5],[2,6],[2,7],[2,8],
                   [3,4],[3,5],[3,6],[3,7],[3,8],
                   [4,5],[4,6],[4,7],[4,8],
                   [5,6],[5,7],[5,8],
                   [6,7],[6,8],
                   [7,8]]
#ENTANGLER_MAP_A= [[0,1],[1,2],[2,5],[5,8],[8,7],[7,6],[6,3]]
ENTANGLER_MAP_A= [[0,1],[1,2],[2,5],[5,8],[7,8],[6,7],[3,6]]
ENTANGLER_MAP_3x4= [[0,1],[2,3],[8,9],[10,11]]
Nx  = 3
Ny = 3
qc = QuantumCircuit(Nx*Ny)
M = 10

for d in range(0,1):
    for i in range(0,Nx*Ny):
        if i in(set(np.array(ENTANGLER_MAP_A).flatten())):
            count = count + 1
            qc.rz(Parameter(parameterBitString(M,count)),i)
    
    qc.barrier()

    for a,b in ENTANGLER_MAP_A:

        count = count + 1
        count = general_SU4(count, qc, a, b,M)
        
for i in range(0,Nx*Ny):
    count = count + 1
    qc.rz(Parameter(parameterBitString(M,count)),i)

print(qc.draw())
print(qc.num_parameters)
p = pd.read_pickle("./pickle/vqe1_3x3_1A+_SU4.pkl")[0]
print(qc.assign_parameters(p["Optimal_params"]))
x_,y_ = extract_N_block(p["Optimal_params"],ENTANGLER_MAP_A, 1,0)
#print(x_,y_)
"""
import random
def randomColor():
    color = random.randrange(0, 2**24)
    hex_color = hex(color)
    std_color = "#" + hex_color[2:]
    return std_color

colors = []
color_label = []
total_x = []
total_y = []
en = []
min_params = []
print("Ladder SU4")
for j in range(1,ladderl+1):
    print("Depth {}".format(j))
    color = randomColor()
    color_label.append([color,"Ladder with {} layers".format(j)])
    en.append([])
    for i in range(1,51):
        p = pd.read_pickle("./pickle/vqe{}_3x3_ladder_{}SU4.pkl".format(i,j))[0]
        x_,y_ = extract_N_block(p["Optimal_params"],ENTANGLER_MAP, j,0)
        en[j-1].append(p["Energy"])
        total_x.append(x_)
        total_y.append(y_)
        for c in range(len(x_)):
            colors.append(color)
    print("Min Energy: {}, Vqe #: {}, Mean: {}, Std: {}".format(np.min(en[j-1]), 
                                                                np.where(en[j-1] == np.min(en[j-1])), 
                                                                np.mean(en[j-1]), np.std(en[j-1])))

print("\nQIDA SU(4)")
for j in range(1,qidal+1):

    print("Depth {}".format(j))
    color = randomColor()
    color_label.append([color,"{} QIDA".format(j)])
    en.append([])
    for i in range(1,51):
        p = pd.read_pickle("./pickle/vqe{}_3x3_{}A_SU4.pkl".format(i,j))[0]
        x_,y_ = extract_N_block(p["Optimal_params"],ENTANGLER_MAP_A, j,0)
        en[j +ladderl-1].append(p["Energy"])
        total_x.append(x_)
        total_y.append(y_)
        for c in range(len(x_)):
            colors.append(color)
    print("Min Energy: {}, Vqe #: {}, Mean: {}, Std: {}".format(np.min(en[j +ladderl-1]),
                                                                np.where(en[j +ladderl-1] == np.min(en[j +ladderl-1])), 
                                                                np.mean(en[j +ladderl-1]), np.std(en[j +ladderl-1])))

print("\nALL-to-ALL Su(4)")
for j in range(1,fullyl+1):
    print("Depth {}".format(j))
    color = randomColor()
    color_label.append([color,"{} Full entanglement".format(j)])
    en.append([])
    for i in range(1,51):
        p = pd.read_pickle("./pickle/vqe{}_3x3_ladder_1F.pkl".format(i))[0]
        x_,y_ = extract_N_block(p["Optimal_params"],ENTANGLER_MAP_F, j,0)
        en[j + ladderl + qidal-1].append(p["Energy"])
        total_x.append(x_)
        total_y.append(y_)
        for c in range(len(x_)):
            colors.append(color)
    print("Min Energy: {}, Vqe #: {}, Mean: {}, Std: {}".format(np.min(en[j + ladderl + qidal-1]),
                                                                np.where(en[j + ladderl + qidal-1] == np.min(en[j + ladderl + qidal-1])),
                                                                np.mean(en[j + ladderl + qidal-1]), np.std(en[j + ladderl + qidal-1])))
print(len(en))

print(len(p["Optimal_params"]))
print(qc.assign_parameters(p["Optimal_params"], inplace=False).draw())
total_x = list(np.concatenate(total_x).flat)
total_y = list(np.concatenate(total_y).flat)

if flag_plot:
    import matplotlib.pyplot as plt
    #plt.scatter()



    plt.figure(figsize=(8, 8))

    if t == "scatter":
        plt.scatter(np.mod(total_x,2*np.pi), np.mod(total_y, 2*np.pi), c=colors, alpha=0.5, s=s)  # Scatter plot of (x, y) pairs
        plt.title('Distribution of (x, y) Pairs: total {} pairs'.format(len(total_y)))

    if t == "hist":
        plt.hist2d(np.mod(total_x,2*np.pi),np.mod(total_y, 2*np.pi), bins=int(bins), cmap='PuRd')
        plt.colorbar()
        plt.title('2D Histogram')

    # Setting x and y ticks to be periodic between -π and π

    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    if int(prec) == 2:
        plt.xticks(np.linspace(0, 2*np.pi, 5), [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
        plt.yticks(np.linspace(0, 2*np.pi, 5), [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])

    if int(prec) == 4:

        plt.xticks(np.linspace(0, 2*np.pi, 9), [r'$0$',   r'$\frac{\pi}{4}$',  r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',r'$\pi}$', 
                                                    r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$',r'$\frac{7\pi}{4}$',r'$2\pi$'])
        plt.yticks(np.linspace(0, 2*np.pi, 9), [r'$0$',   r'$\frac{\pi}{4}$',  r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$',r'$\pi}$', 
                                                    r'$\frac{5\pi}{4}$', r'$\frac{3\pi}{2}$',r'$\frac{7\pi}{4}$',r'$2\pi$'])

    if int(prec) == 8:
        plt.xticks(np.linspace(0, 2*np.pi, 17), [r'$0$',  r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$', r'$\frac{3\pi}{8}$', r'$\frac{\pi}{2}$', r'$\frac{5\pi}{8}$',
                                                r'$\frac{3\pi}{4}$', r'$\frac{7\pi}{8}$', r'$\pi}$', r'$\frac{9\pi}{8}$', r'$\frac{5\pi}{4}$', r'$\frac{11\pi}{8}$',
                                                r'$\frac{3\pi}{2}$', r'$\frac{13\pi}{8}$', r'$\frac{7\pi}{4}$', r'$\frac{15\pi}{8}$', r'$2\pi$'])
        plt.yticks(np.linspace(0, 2*np.pi, 17), [r'$0$',  r'$\frac{\pi}{8}$', r'$\frac{\pi}{4}$', r'$\frac{3\pi}{8}$', r'$\frac{\pi}{2}$', r'$\frac{5\pi}{8}$',
                                                r'$\frac{3\pi}{4}$', r'$\frac{7\pi}{8}$', r'$\pi}$', r'$\frac{9\pi}{8}$', r'$\frac{5\pi}{4}$', r'$\frac{11\pi}{8}$',
                                                r'$\frac{3\pi}{2}$', r'$\frac{13\pi}{8}$', r'$\frac{7\pi}{4}$', r'$\frac{15\pi}{8}$', r'$2\pi$'])



    # Display the plot

    from matplotlib.lines import Line2D
    def createLegend(cl):
        legend = []
        for item,i in enumerate(cl):
            print(item,i)
            legend.append(Line2D([0], [0], marker='o', color=cl[item][0], label='{}'.format(i[1]), markersize=15))
        return legend

    plt.legend(handles=createLegend(color_label))
    plt.grid(True)
    plt.show()
"""