import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



file_qida = "./pickle/vqe{}_3x3_{}A_SU4_ord.pkl"
file_ladder = "./pickle/vqe{}_3x3_ladder_{}SU4.pkl"

def df_creation(f,size):
    df = []
    for i in range(1,size+1):
        p = pd.read_pickle(f.format(i,1))[0]
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

def create_plot(size,title,x_,y_,s, file_png):
    plt.figure(figsize=size)
    plt.title(title)
    plt.xlim(0, len(en)+1)
    plt.xlabel("Index of VQE")
    plt.ylabel("Energy")
    plt.axhline(-4.7493272585529, c='b')
    plt.axhline(-3, c='r')
    #plt.scatter(x=x_, y = [-4.7493272585529]*len(x_), c='r')
    #plt.scatter(x=x_, y = [-3]*len(x_), c='b')
    plt.scatter(x=x_, y=y_, s=s, alpha=1)
    plt.savefig(file_png)

mode = 'LADDER'
depth = 1
size = 50

f_save  = "./plots/{}-{}Layer-SU(4)-{}VQEs.png".format(mode, depth, size)

df = df_creation(f=file_ladder, size=size)
i,en, min_vqe, min_vqe_idx, vqe_std, vqe_var, vqe_mean= plt_energies(data = df)
s = size_plot(min_vqe_idx, s_min = 40, s_all=10, l = len(en))

print(min_vqe, min_vqe_idx)

title = "{}-{}Layer-SU(4)-{}VQEs\nmin_vqe: {} (i={})\nstd: {} / var: {}\n(mean = {})".format(mode,depth,size,min_vqe, i[min_vqe_idx[0][0]], vqe_std, vqe_var, vqe_mean)

create_plot(size=(8,8), title=title,x_=i, y_=en, s=s, file_png=f_save)
