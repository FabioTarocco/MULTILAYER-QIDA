import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import random

depth = 1
#file_qida = "./MULTIQIDA/VQE_res/3x3_{}/vqe{}_3x3_{}A.pkl"
#file_ladder = "./ladder/3x3/ladder_{}/vqe{}_3x3_ladder_{}SU4.pkl"

file_qida = "./pickle/vqe{}_3x3_{}A_SU4_ord.pkl"
file_ladder = "./pickle/vqe{}_3x3_ladder_{}SU4.pkl"
file_7 = "./pickle/vqe{}_3x3_{}A_SU4.pkl"

mode = 'QIDA'#'QIDA'
gates = 'SU4'#'SU(4)'
size = 500

f_save  = "./plots/{}-1-Layer-{}-{}VQEs.png".format(mode, depth,gates, size)

def createLegend(cl):
    legend = []

    for item,i in enumerate(cl):
        legend.append(Line2D([0], [0], marker='o', color=cl[item], label="{} layers".format(item+1), markersize=15))
    return legend


def randomColor():
    color = random.randrange(0, 2**24)
    hex_color = hex(color)
    std_color = "#" + hex_color[2:]
    return std_color

def df_creation(f,size,depth):
    df = []
    for i in range(1,size+1):
        #p = pd.read_pickle(f.format(depth,i,depth))[0]
        p = pd.read_pickle(f.format(i,depth))[0]
        data = dict()
        data.update({"Params": p["Optimal_params"], "Energy":p['Energy'], "Index": i })
        df.append(data)
    return df

def stats(e):

    return np.min(e), np.where(e == np.min(e)), np.std(e), np.var(e), np.mean(e)

def plt_energies(data:dict):
    """
        data -> indices of vqes, vqes, min_vqe, index_best_vqe, std, var, mean
    """
    i = []
    en = []
    for item in data:
        i.append(item['Index'])
        en.append(item['Energy'])

    min_vqe, min_vqe_idx, vqe_std, vqe_var, vqe_mean = stats(en)
    return i, en, min_vqe, min_vqe_idx, vqe_std, vqe_var, vqe_mean

def size_plot (idx, s_min, s_all, l):
    s = []
    for i in range(l):
        if i in idx:
            s.append(s_min)
        else:
            
            s.append(s_all)
    return s


def create_plot(size,title,n_elem):
    plt.figure(figsize=size)
    plt.title(title)
    plt.xlim(0, n_elem+1)
    plt.xlabel("Index of VQE")
    plt.ylabel("Energy")
    plt.axhline(-4.7493272585529, c='b')
    plt.axhline(-3, c='r')

def add_things(x_,y_,s,c):
    plt.scatter(x=x_, y=y_, s=s, alpha=1,c=c)

mins = ""
title = "{}-(1-3)Layer-{}-{}VQEs".format(mode,gates,size)

create_plot(size=(8,8), title=title, n_elem = size)
colors = ['black', 'brown', 'red', 'orange', 'yellow', 'lime', 'green'][::-1]
for d in range(1,depth+1):
    df = df_creation(f=file_qida, size=size, depth = d)
    i,en, min_vqe, min_vqe_idx, vqe_std, vqe_var, vqe_mean = plt_energies(data = df)
    s = size_plot(min_vqe_idx, s_min = 40, s_all=10, l = len(en))
    mins+="Depth:{}\nmin_vqe: {} (i={})\nstd: {} / var: {}\n(mean: {})\n%QE: {}\n-------\n".format(d,min_vqe, i[min_vqe_idx[0][0]], vqe_std, vqe_var, vqe_mean,100*((min_vqe - (-3.0))/(-4.7493272585529 - (-3.0))))

    add_things(x_=i, y_=en, s=s, c=colors[d-1])

plt.legend(handles=createLegend(colors[:1:]))
plt.savefig(f_save)

f_save_txt  = "./plots/{}-1-Layer-{}-{}VQEs.txt".format(mode,gates,size)
f = open(f_save_txt, "w")
f.write(str(mins))
f.close()
