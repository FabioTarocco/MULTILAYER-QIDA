import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pprint import pprint
from matplotlib.patches import Patch




base_path = "./MULTIQIDA/VQE_res/"
add_path = ["3x3/{}_3x3_4LCX.pkl","3x3/{}_3x3_5LCX.pkl","3x3/{}_3x3_conv_1ACX_1VBCX_1VLCX.pkl","3x3/{}_3x3_conv_1ASO4_1BSO4_1LSO4.pkl",
            "2x6/{}_2x6_4LCX.pkl","2x6/{}_2x6_5LCX.pkl","2x6/{}_2x6_conv_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","2x6/{}_2x6_conv_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4/{}_3x4_4LCX.pkl","3x4/{}_3x4_5LCX.pkl","3x4/{}_3x4_6LCX.pkl","3x4/{}_3x4_conv_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4/{}_3x4_conv_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_mag/{}_3x4_mag_4LCX.pkl","3x4_mag/{}_3x4_mag_5LCX.pkl","3x4_mag/{}_3x4_mag_conv_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_mag/{}_3x4_mag_conv_1ASO4_1BSO4_1CSO4_1LSO4.pkl",
            "3x4_delta/{}_3x4_delta_4LCX.pkl","3x4_delta/{}_3x4_delta_5LCX.pkl","3x4_delta/{}_3x4_delta_6LCX.pkl","3x4_delta/{}_3x4_delta_conv_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4_delta/{}_3x4_delta_conv_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_delta_0_1/{}_3x4_delta_0_1_4LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_5LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_conv_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_conv_1ASO4_1BSO4_1CSO4_1LSO4.pkl"]
path_complete = []
for p in add_path:
    path_complete.append(base_path+p)
#print(path_complete)
NVQE  = 50
n_configs = [4,4,5,4,5,4]
n_systems = ["3x3","2x6","3x4_i","3x4_mag","3x4_delta","3x4_delta_0_1"]
conf_layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
conf_ed = [1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6]

titles = {  "3x3":"System size: 3x3\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "2x6":"System size: 2x6\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "3x4_i":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "3x4_mag":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=2.0, Anisotropic term: $\Delta$=0.0",
            "3x4_delta":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=2/3",
            "3x4_delta_0_1":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$, Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=1/10"}
c1 = ['4LCX','5LCX','QIDACX','QIDASO4']
c1_map = ['inferno_r','rocket_r','crest','viridis_r']
c2 = ['4LCX','5LCX','6LCX','QIDACX','QIDASO4']

system_complete = []
for ns in range(len(n_systems)):
    for _ in range(NVQE*n_configs[ns]):
        system_complete.append(n_systems[ns])



# Data creation
        
e_ed = {"3x3":-4.7493272585529,"2x6":-6.60347247538703,"3x4_i":-6.69168019351493,"3x4_mag":-9.50847255610575, "3x4_delta":-5.33875124065748,"3x4_delta_0_1":-4.272670379839620}
layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
configs = ['4LCX','5LCX','QIDACX','QIDASO4','4LCX','5LCX','QIDACX','QIDASO4','4LCX','5LCX','6LCX','QIDACX','QIDASO4','4LCX','5LCX','QIDACX','QIDASO4','4LCX','5LCX','6LCX','QIDACX','QIDASO4','4LCX','5LCX','QIDACX','QIDASO4']

conv =dict()
for cp in range(len(path_complete)):
    experiments = []
    if "3x3" not in path_complete[cp]:
        continue
    #print("OPENING FILE {}".format(path_complete[cp]))
    energies = []
    sum = []
    for i in range(1,NVQE+1):
        runs = []
        vqe = pd.read_pickle(path_complete[cp].format(i))[0]
        if layers[cp] == 1:
            #sum.append(len(vqe["Convergence"]))
            runs.append(vqe["Convergence"])
            energies.append(vqe["Energy"])
        if layers[cp]==2:
            temp = [item for inner_list in vqe["Convergence_12"] for item in inner_list]
            runs.append(temp)
            #sum.append(len(temp))
            energies.append(vqe["Energy_12"])
        if layers[cp]==3:
            temp = [item for inner_list in vqe["Convergence_123"] for item in inner_list]
            runs.append(temp)
            #sum.append(len(temp))
            energies.append(vqe["Energy_123"])
        if layers[cp]==4:
            temp = [item for inner_list in vqe["Convergence_1234"] for item in inner_list]
            runs.append(temp)
            #sum.append(len(temp))
            energies.append(vqe["Energy_1234"])
        if layers[cp]==5:
            temp = [item for inner_list in vqe["Convergence_12345"] for item in inner_list]
            runs.append(temp)
            #sum.append(len(temp))
            energies.append(vqe["Energy_12345"])
        experiments.append(runs[0])
    conv.update({n_systems[conf_ed[cp]-1]+"_"+configs[cp]:[experiments, np.min(energies), np.argwhere(np.array(energies)==np.min(energies))[0][0],energies,e_ed[n_systems[conf_ed[cp]-1]], titles[n_systems[conf_ed[cp]-1]]]})


#Data Pre-Process

selection = ['3x3_4LCX','3x3_5LCX','3x3_QIDACX', '3x3_QIDASO4']

fig, ax = plt.subplots(len(selection),1,figsize=(8,6), sharex=True)
plt.subplots_adjust(hspace=0)


for s_i,sel in enumerate(selection):

    ax[s_i].set_ylabel("AQE(%)",fontsize=10)
    ax[s_i].set_xlabel("Iterations", fontsize=12)
    #ax[s_i].set_yticklabels([60,80,100], fontsize = 10)
    data,min,min_index,en,ed_min,title = conv[sel]

    max_iters = np.max(np.array([len(i) for i in data]))
    min_iters = np.min(np.array([len(i) for i in data]))
    max_iter_list = np.argmax(np.array([len(i) for i in data]))
    min_iter_list = np.argmin(np.array([len(i) for i in data]))
    #print(max_iters, max_iter_list, min_iters, min_iter_list, len(data[max_iter_list]), len(data[min_iter_list]))
    ticksx = np.arange(0,max_iters,1)
    all_ticksx = []

    for i in data:
        all_ticksx.append(np.arange(0,len(i),1))

    data_p = []
    for i in data:
        temp = []
        for r in i:
            temp.append((r/ed_min)*100)
        data_p.append(temp)
    #for i in range(len(data)):
    #    print(len(data_p[i]), len(data[i]), len(all_ticksx[i]))
    #print(len(data_p), len(data), len(all_ticksx))

    filled_nan_data = np.full((len(data_p),max_iters),np.nan)
    for i, exp in enumerate(data_p):
        filled_nan_data[i,:len(exp)] = exp

    for i in filled_nan_data:
        print(len(i))

    #print("MINIMUM: {}({}) -{}({})".format(min, min_index,data[min_iter_list][-1],min_iter_list))

    #Data visualization


    #print(max_iters, len(ticksx))
    #ax[s_i].set_xlim(0,max_iters) DA CONTROLLARE

    #print(len(data[min_index]))
    #print("Min Iters: {} - Max Itets: {}".format(min_iters,max_iters))
    #print("Min En: {} - Min (NP): {} - Index: {}\n\tMean: {}".format(min, np.min(en),min_index, np.mean(en)))
    #print("ED {}".format(ed_min))
    #ax[s_i].axhline(y=(ed_min/ed_min)*100, linestyle='--', c='black')
    #ax[s_i].axhline(y=(min/ed_min)*100, linestyle='--')


    l = max_iters#len(all_ticksx[min_index])
    env_min = np.nanmin(filled_nan_data,axis=0)
    env_max = np.nanmax(filled_nan_data,axis=0)
    ax[s_i].fill_between(ticksx[:l],env_min,env_max,alpha=0.2,edgecolor =sns.color_palette("mako")[-3],color=sns.color_palette("mako")[-2])
    # for y,x in zip(data_p,all_ticksx):
    #     print(len(x),len(y))
    #     ax[s_i].plot(x[:l],y[:l], c=sns.color_palette("mako")[-1],alpha = 0.2)
    #ax[s_i].plot(all_ticksx[min_iter_list],data_p[min_iter_list],c=c_min_iter)
    ax[s_i].plot(ticksx[:l],np.nanmean(filled_nan_data,axis=0)[:l],c=sns.color_palette("viridis")[2])
    ax[s_i].plot(all_ticksx[min_index][:l],data_p[min_index][:l], c=sns.color_palette("viridis")[0])
    ax[s_i].label_outer()
    ax[s_i].set_ylim(60,102)
    ax[s_i].invert_yaxis()


    l_elems = [Patch(facecolor=sns.color_palette("mako")[-1],alpha = 0.2, label = 'All runs'),
               Patch(facecolor=sns.color_palette("viridis")[2], label = '$VQE_{Mean}$'),
               Patch(facecolor=sns.color_palette("viridis")[0], label = '$VQE_{Min}$')]
    ax[s_i].legend(handles = l_elems,loc='upper right',fontsize =7)
plt.suptitle("Convergence\n"+title, fontsize=12,y=0.95)
plt.tight_layout()
plt.savefig("3x3.pdf", format='pdf')
plt.show()





        
