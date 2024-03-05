import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt 
from matplotlib.patches import Patch
import seaborn as sns

base_path = "./MULTIQIDA/VQE_res/"
add_path = ["3x3/{}_3x3_4LCX.pkl","3x3/{}_3x3_5LCX.pkl","3x3/{}_3x3_1ACX_1VBCX_1VLCX.pkl","3x3/{}_3x3_1ASO4_1BSO4_1LSO4.pkl",
            "2x6/{}_2x6_4LCX.pkl","2x6/{}_2x6_5LCX.pkl","2x6/{}_2x6_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","2x6/{}_2x6_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4/{}_3x4_4LCX.pkl","3x4/{}_3x4_5LCX.pkl","3x4/{}_3x4_6LCX.pkl","3x4/{}_3x4_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4/{}_3x4_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_mag/{}_3x4_mag_4LCX.pkl","3x4_mag/{}_3x4_mag_5LCX.pkl","3x4_mag/{}_3x4_mag_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_mag/{}_3x4_mag_1ASO4_1BSO4_1CSO4_1LSO4.pkl",
            "3x4_delta/{}_3x4_delta_4LCX.pkl","3x4_delta/{}_3x4_delta_5LCX.pkl","3x4_delta/{}_3x4_delta_6LCX.pkl","3x4_delta/{}_3x4_delta_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4_delta/{}_3x4_delta_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_delta_0_1/{}_3x4_delta_0_1_4LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_5LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ASO4_1BSO4_1CSO4_1LSO4.pkl"]
conf_layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
conf_ed = [1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6]

ansatz = ["$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$",
          "$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$",
          "$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$(L)_{6}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$",
          "$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$",
          "$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$(L)_{6}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$",
          "$(L)_{4}^{CX}$","$(L)_{5}^{CX}$","$QIDA^{CX}$","$QIDA^{SO4}$"]
#TOLTO 3x4_6LCX un 3 in conf_ed ed un 1 in layers

NVQE = 50
systems = ["3x3","2x6","3x4","3x4_mag","3x4_delta","3x4_delta_0_1"]
label_system = ["3x3", "2x6", "3x4", "3x4 $h=2.0$", "3x4 $\Delta=2/3$","3x4 $\Delta=1/10$"]
e_ed = [-4.7493272585529,-6.60347247538703,-6.69168019351493,-9.50847255610575,-5.33875124065748,-4.272670379839620]
neel = [-3.0,-4.0,-4.25,-4.25,-4.25,-4.25]

cm = [sns.color_palette("Spectral")[0],sns.color_palette("mako")[-2], sns.color_palette("magma")[1], sns.color_palette("magma")[3], sns.color_palette("viridis")[-2], sns.color_palette("Set2")[5]]

cm = [sns.color_palette("Set2")[0],sns.color_palette("Set2")[1], sns.color_palette("Set2")[2], sns.color_palette("Set2")[3], sns.color_palette("Set2")[4], sns.color_palette("Set2")[5]]
systems_configs = []
for item in zip(add_path, conf_layers, conf_ed,ansatz):
    systems_configs.append((base_path+item[0],item[1],e_ed[item[2]-1],neel[item[2]-1],item[3],cm[item[2]-1]))

#pprint(systems_configs)

out = ""

fig, axs = plt.subplots(figsize=(6,4))
plt.title(" MRED vs. RQE", fontsize = 20)
axs.set_xlabel("RQE(%)", fontsize =18)
axs.set_ylabel("MRED(%)", fontsize =18)
axs.set_yscale("log")
axs.yaxis.tick_right()
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
axs.yaxis.set_label_position("right")
l_elems = []

min_ladder = 1000
max_qida = -1
result = pd.DataFrame(columns=['Config','E_ed','E_neel','E_mean','Std','Var','E_min','AQE_mean','RQE_mean','AQE_min','RQE_min','MED','MAED','MRED'])
for c,l,e,n,label,cmap in systems_configs:
    energies = []
    rqe = []
    aqe = []
    med = 0.0
    mred = 0.0
    maed = 0.0

    print("OPENING FILE {}".format(c))
    print("SYSTEM ED:{}".format(e))
    for i in range(1,NVQE+1):
        vqe = pd.read_pickle(c.format(i))[0]
        if l == 1:
            en_i = vqe["Energy"]
        if l == 2:
            en_i = vqe["Energy_12"]
        if l == 3:
            en_i = vqe["Energy_123"]
        if l == 4:
            en_i = vqe["Energy_1234"]
        if l == 5:
            en_i = vqe["Energy_12345"]
        if en_i > n:
            energies.append(n)
            rqe.append(0.0)
            aqe.append(100*(n/e))
        else:
            energies.append(en_i)
            rqe.append(100*((en_i - n)/(e - n)))
            aqe.append(100*(en_i/e))
        

    e_min = np.min(energies)
    rqe_min = np.max(rqe)
    aqe_min = np.max(aqe)

    for j in energies:
        med += (np.abs((j - e_min)))/NVQE
    for j in rqe:
        mred += (np.abs((j - rqe_min)))/NVQE
    for j in aqe:
        maed += (np.abs((j - aqe_min)))/NVQE
    
    o = "RESULTS for {}\n\tED(E): {}\tNeel(E):{}\n\tMEAN VQE(E): {}\n\tStd: {}\n\tVar: {}\n\tMIN VQE(E): {}\n\tMEAN AQE(%): {}\n\tAQE BEST(%): {}\n\tMEAN RQE(%): {}\n\tRQE BEST(%): {}\n\tMED(E): {}\n\tMAED(%): {}\n\tMRED(%): {}\n\n".format(c.format(i),e,n, np.mean(energies), np.std(energies), np.var(energies), e_min,np.mean(aqe),aqe_min, np.mean(rqe), rqe_min, med,maed, mred)
    result.loc[len(result.index)] = {'Config':c,'E_ed':e,'E_neel':n,'E_mean':np.mean(energies),'Std':np.std(energies),'Var':np.var(energies),'E_min':np.min(energies),'AQE_mean':np.mean(aqe),'RQE_mean':np.mean(rqe),'AQE_min':aqe_min,'RQE_min':rqe_min,'MED':med,'MAED':maed,'MRED':mred}
    print(o)
    out += o
    if "QIDA" in label:
        max_qida = mred if max_qida < mred else max_qida
            
        if "SO4" in label:
            axs.scatter(x=np.mean(rqe),y=mred, s=35, marker='x',color=cmap)
            print("x: {}//{}\n--------------------\n".format(np.mean(rqe), mred))
        else:
            axs.scatter(x=np.mean(rqe),y=mred, s=35, marker='+',color=cmap)
            print("x: {}//{}\n--------------------\n".format(np.mean(rqe), mred))
    else:

        min_ladder = mred if min_ladder > mred else min_ladder
        axs.scatter(x=np.mean(rqe),y=mred, s=30, marker='o',color=cmap)
        print("o: {}//{}\n---------------------\n".format(np.mean(rqe), mred))
    #axs.scatter(x=rqe_min,y=mred_error, s=20, marker='|',color=cmap)
    #axs.annotate(label,xy=(np.mean(rqe),mred_error/2))
    
print(min_ladder, max_qida)
for label_name, patch_col in zip(label_system, cm):
    l_elems.append(Patch(facecolor=patch_col, label = label_name))
#axs.hlines(xmin=40, xmax=100.1,y=1, linestyles={'dotted'}, colors='black')
axs.hlines(xmin=40, xmax=100.1,y=min_ladder, linestyles={'dotted'}, colors='gray',alpha = 0.2)
axs.hlines(xmin=40, xmax=100.1,y=max_qida, linestyles={'dotted'}, colors='gray', alpha = 0.2)
plt.legend(handles = l_elems,loc='lower left',fontsize = 12)
plt.tight_layout()
plt.savefig("comparison_mixed.pdf")
plt.show()
exit(0)
with open("aned_recap.txt", 'w') as file:
    file.write(out)

with open("aned_recap_df.txt", 'w') as file:
    result.to_csv('aned_recap_df.txt', sep='\t',index =False)







