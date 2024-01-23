import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Patch
from pprint import pprint


base_path = "./MULTIQIDA/VQE_res/"
add_path = ["3x3/{}_3x3_4LCX.pkl","3x3/{}_3x3_5LCX.pkl","3x3/{}_3x3_1ACX_1VBCX_1VLCX.pkl","3x3/{}_3x3_1ASO4_1BSO4_1LSO4.pkl",
            "2x6/{}_2x6_4LCX.pkl","2x6/{}_2x6_5LCX.pkl","2x6/{}_2x6_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","2x6/{}_2x6_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4/{}_3x4_4LCX.pkl","3x4/{}_3x4_5LCX.pkl","3x4/{}_3x4_6LCX.pkl","3x4/{}_3x4_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4/{}_3x4_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_mag/{}_3x4_mag_4LCX.pkl","3x4_mag/{}_3x4_mag_5LCX.pkl","3x4_mag/{}_3x4_mag_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_mag/{}_3x4_mag_1ASO4_1BSO4_1CSO4_1LSO4.pkl",
            "3x4_delta/{}_3x4_delta_4LCX.pkl","3x4_delta/{}_3x4_delta_5LCX.pkl","3x4_delta/{}_3x4_delta_6LCX.pkl","3x4_delta/{}_3x4_delta_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4_delta/{}_3x4_delta_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_delta_0_1/{}_3x4_delta_0_1_4LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_5LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ASO4_1BSO4_1CSO4_1LSO4.pkl"]
path_complete = []
for p in add_path:
    path_complete.append(base_path+p)
print(path_complete)
NVQE  = 50
n_configs = [4,4,5,4,5,4]
n_systems = ["3x3","2x6","3x4","3x4_mag","3x4_delta","3x4_delta_0_1"]
c1 = ['4LCX','5LCX','QIDACX','QIDASO4']
c2 = ['4LCX','5LCX','6LCX','QIDACX','QIDASO4']

system_complete = []
for ns in range(len(n_systems)):
    for _ in range(NVQE*n_configs[ns]):
        system_complete.append(n_systems[ns])
    
print(system_complete.count("3x3"))

    
configs_complete = []
for ns in range(len(n_systems)):
    if n_configs[ns]==4:
        for cf in c1:
            for _ in range(NVQE):
                configs_complete.append(cf)
    else:
        for cf in c2:
            for _ in range(NVQE):
                configs_complete.append(cf)

print(configs_complete)

value_complete = np.repeat(np.random.exponential(size=NVQE),np.sum(n_configs))
layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
value_complete = []
for cp in range(len(path_complete)):
    print("OPENING FILE {}".format(path_complete[cp]))
    for i in range(1,NVQE+1):
        vqe = pd.read_pickle(path_complete[cp].format(i))[0]
        if layers[cp] == 1:
            en_i = vqe["Energy"]
            value_complete.append(en_i)
        if layers[cp]==2:
            en_i = vqe["Energy_12"]
            value_complete.append(en_i)
        if layers[cp]==3:
            en_i = vqe["Energy_123"]
            value_complete.append(en_i)
        if layers[cp]==4:
            en_i = vqe["Energy_1234"]
            value_complete.append(en_i)
        if layers[cp]==5:
            en_i = vqe["Energy_12345"]
            value_complete.append(en_i)


df = pd.DataFrame(data={'En':value_complete, 'Config':configs_complete,'System':system_complete})



e_ed = {"3x3":-4.7493272585529,"2x6":-6.60347247538703,"3x4":-6.69168019351493,"3x4_mag":-9.50847255610575, "3x4_delta":-5.33875124065748,"3x4_delta_0_1":-4.272670379839620}

e_neel = {"3x3":-3.0,"2x6":-4-0,"3x4":-4.25,"3x4_mag":-4.25, "3x4_delta":-4.25,"3x4_delta_0_1":-4.25}
titles = {  "3x3":"System size: 3x3\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "2x6":"System size: 2x6\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "3x4":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=0.0",
            "3x4_mag":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=2.0, Anisotropic term: $\Delta$=0.0",
            "3x4_delta":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=2/3",
            "3x4_delta_0_1":"System size: 3x4\nInteraction term: $J_{xx}=J_{yy}=J_{zz}=1.0$\n Mag. Field: $h$=0.0, Anisotropic term: $\Delta$=1/10"}
print(df)
for sel in n_systems:
    df_sel = df[df['System']==sel]
    sel_min = e_ed[sel]
    print(sel_min)
    df_sel_normalized = df_sel.copy()
    
    
    cm1=sns.color_palette("viridis")
    cm2=sns.color_palette("inferno")
    if len(df_sel_normalized)==200:
        df_x1 = pd.DataFrame({'E_p':(df_sel_normalized['En']/(sel_min))*100})#,'t':np.repeat('AQE',200),'Config':df_sel_normalized['Config'],'Pos':np.repeat([0.0,3.0,6.0,9.0],50)})
        df_x2 = pd.DataFrame({'E_p':((df_sel_normalized['En'] - e_neel[sel])/(e_ed[sel]-e_neel[sel]))*100})#,'t':np.repeat('RQE',200),'Config':df_sel_normalized['Config'],'Pos':np.repeat([1.0,4.0,7.0,10.0],50)})

    else:
        df_x1 = pd.DataFrame({'E_p':(df_sel_normalized['En']/(sel_min))*100})#,'t':np.repeat('AQE',250),'Config':df_sel_normalized['Config'],'Pos':np.repeat([0.0,3.0,6.0,9.0,12.0],50)})
        df_x2 = pd.DataFrame({'E_p':((df_sel_normalized['En'] - e_neel[sel])/(e_ed[sel]-e_neel[sel]))*100})#,'t':np.repeat('RQE',250),'Config':df_sel_normalized['Config'],'Pos':np.repeat([0.0,3.0,6.0,9.0,12.0],50)})

    
    fig, ax1 = plt.subplots(figsize=(8,6))
    ax1.set_ylabel("AQE(%)",fontsize=12)
    ax1.set_ylim(55,100.5)
    ax1.set_yticks([70,80,85,90,92.5,95,97.5,100])
    ax1.set_xlabel("Configurations",fontsize = 12)
    if len(df_sel)==200:
        data = [df_x1['E_p'].iloc[0:49],
                df_x1['E_p'].iloc[50:99],
                df_x1['E_p'].iloc[100:149],
                df_x1['E_p'].iloc[150:199],
                ]
        
        pos = np.linspace(0,20,num=6)
        ax1.set_xlim(1.75,18.25)
        ticks = []
        offset = 0.8
        for item in pos[1:-1]:
            ticks.append(item-offset)
            ticks.append(item+offset)
        ax1.set_xticks(pos[1:-1])

        ax1.set_xticklabels(c1)
        ax1.xaxis.grid(True)
        ax1=plt.violinplot(dataset=data, positions=pos[1:-1],widths=3, showextrema=True, showmeans=True, showmedians=False)
        colors = []
        for nv in range(len(c1)):
            colors.append(cm1[nv])
    else:
        data = [df_x1['E_p'].iloc[0:49],
                df_x1['E_p'].iloc[50:99],
                df_x1['E_p'].iloc[100:149],
                df_x1['E_p'].iloc[150:199],
                df_x1['E_p'].iloc[200:249]
                ]
        
        pos = np.linspace(0,20,num=7)
        ax1.set_xlim(1.75,18.25)
        ticks = []
        offset = 0.7
        for item in pos[1:-1]:
            ticks.append(item-offset)
            ticks.append(item+offset)
        ax1.set_xticks(pos[1:-1])

        ax1.set_xticklabels(c2)
        ax1.xaxis.grid(True)
        ax1=plt.violinplot(dataset=data, positions=pos[1:-1],widths=1.3, showextrema=True, showmeans=True, showmedians=False)
        
        colors = []
        for nv in range(len(c2)):
            colors.append(cm1[nv])

    for pc, color in zip(ax1['bodies'],colors):
        pc.set_facecolor(color)
    ax1['cmeans'].set_colors(colors)
    ax1['cmaxes'].set_colors(colors)
    ax1['cmins'].set_colors(colors)
    ax1['cbars'].set_colors(colors)
    line = ax1['cmaxes'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([1.4,1])+center
        new_line.append(l)
    ax1['cmaxes'].set_segments(new_line)
    line = ax1['cmeans'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([0.8,1])+center
        new_line.append(l)
    ax1['cmeans'].set_segments(new_line)
    if len(df_sel_normalized)==200:
        plt.legend(c1,loc='lower right')
    else:
        plt.legend(c2, loc='lower right')
    title = "System:{}".format(sel)+"\n$E_{ED}$:"+" {}".format(sel_min)
    plt.title(titles[sel],fontsize=12)
    plt.tight_layout()
    plt.show()
    plt.savefig(sel+"_aqe.pdf")
    plt.close('all')

