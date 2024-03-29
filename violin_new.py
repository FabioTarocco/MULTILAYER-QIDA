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
#c1 = ['4LCX','5LCX','QIDACX','QIDASO4']
#c2 = ['4LCX','5LCX','6LCX','QIDACX','QIDASO4']

c1 = ['$(L)_{4}^{CX}$','$(L)_{5}^{CX}$','$QIDA^{CX}$','$QIDA^{SO4}$']
c2 = ['$(L)_{4}^{CX}$','$(L)_{5}^{CX}$','$(L)_{6}^{CX}$','$QIDA^{CX}$','$QIDA^{SO4}$']

system_complete = []
for ns in range(len(n_systems)):
    for _ in range(NVQE*n_configs[ns]):
        system_complete.append(n_systems[ns])
    

    
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

layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
value_complete = []
for cp in range(len(path_complete)):

    print("OPENING FILE {}".format(path_complete[cp]))

    for i in range(1,NVQE+1):
        vqe = pd.read_pickle(path_complete[cp].format(i))[0]
        if layers[cp] == 1:
            en_i = vqe["Energy"]
        if layers[cp]==2:
            en_i = vqe["Energy_12"]
        if layers[cp]==3:
            en_i = vqe["Energy_123"]
        if layers[cp]==4:
            en_i = vqe["Energy_1234"]
        if layers[cp]==5:
            en_i = vqe["Energy_12345"]
        value_complete.append(en_i)

print(len(value_complete), len(configs_complete), len(system_complete))
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
    print(len(df_x2[df_x2['E_p']<0]))
    df_x2['E_p'] = df_x2['E_p'].clip(lower=0)


    fig, ax1 = plt.subplots(figsize=(8,6))
    plt.yticks(fontsize=15)
    ax2 = ax1.twinx()
    ax1.set_ylabel("AQE(%)",fontsize=18,color=cm1[0])
    ax1.set_yticks([50,75,80,85,90,92.5,95,97.5,100])
    ax1.set_ylim(np.min(df_x1)-1,np.max(df_x1)+1)
    ax1.set_xlabel("Configurations",fontsize = 18)
    ax2.set_ylabel("RQE(%)",fontsize=18,labelpad=17,color=cm1[3],rotation=-90)
    ax2.set_yticks([0,25,50,75,80,85,90,95,100])
    ax2.set_ylim(np.min(df_x2)-1,np.max(df_x2)+1)
    plt.yticks(fontsize=15)
    if len(df_sel)==200:
        data1 = [df_x1['E_p'].iloc[0:49],
                df_x1['E_p'].iloc[50:99],
                df_x1['E_p'].iloc[100:149],
                df_x1['E_p'].iloc[150:199]]
        
        data2 = [df_x2['E_p'].iloc[0:49],
                df_x2['E_p'].iloc[50:99],
                df_x2['E_p'].iloc[100:149],
                df_x2['E_p'].iloc[150:199]]
        
        pos = np.linspace(0,20,num=6)
        ax1.set_xlim(1.75,18.25)
        ticks1 = []
        ticks2 = []
        offset = 0.9
        for item in pos[1:-1]:
            ticks1.append(item-offset)
            ticks2.append(item+offset)
        ax1.set_xticks(pos[1:-1])

        ax1.set_xticklabels(c1, fontsize = 17)
        ax1.xaxis.grid(True)
        v1=ax1.violinplot(dataset=data1, positions=ticks1,widths=1.7, showextrema=True, showmeans=True, showmedians=False)
        ax1.tick_params(axis='y', labelcolor=cm1[0])
        v2=ax2.violinplot(dataset=data2, positions=ticks2,widths=1.7, showextrema=True, showmeans=True, showmedians=False)
        ax2.tick_params(axis='y', labelcolor=cm1[3])
        colors1 = []
        colors2 = []
        for nv in range(len(c1)):
            colors1.append(cm1[0])
            colors2.append(cm1[3])
        #colors = ['Purple','Orange','Purple','Orange','Purple','Orange','Purple','Orange']
    else:
        data1 = [df_x1['E_p'].iloc[0:49],
                df_x1['E_p'].iloc[50:99],
                df_x1['E_p'].iloc[100:149],
                df_x1['E_p'].iloc[150:199],
                df_x1['E_p'].iloc[200:249]]
        data2 = [
                df_x2['E_p'].iloc[0:49],
                df_x2['E_p'].iloc[50:99],
                df_x2['E_p'].iloc[100:149],
                df_x2['E_p'].iloc[150:199],
                df_x2['E_p'].iloc[200:249]]
        pos = np.linspace(0,20,num=7)
        ax1.set_xlim(1.75,18.25)
        ticks1 = []
        ticks2 = []
        offset = 0.8
        for item in pos[1:-1]:
            ticks1.append(item-offset)
            ticks2.append(item+offset)
        ax1.set_xticks(pos[1:-1])
        ax1.set_xticklabels(c2, fontsize = 17)
        ax1.xaxis.grid(True)
        v1=ax1.violinplot(dataset=data1, positions=ticks1,widths=1.5, showextrema=True, showmeans=True, showmedians=False)
        ax1.tick_params(axis='y', labelcolor=cm1[0])
        v2=ax2.violinplot(dataset=data2, positions=ticks2,widths=1.5, showextrema=True, showmeans=True, showmedians=False)
        ax2.tick_params(axis='y', labelcolor=cm1[3])
        
        colors1 = []
        colors2 = []
        for nv in range(len(c2)):
            colors1.append(cm1[0])

            colors2.append(cm1[3])
        #colors = ['Purple','Orange','Purple','Orange','Purple','Orange','Purple','Orange','Purple','Orange']
            
    for pc, color in zip(v1['bodies'],colors1):
        pc.set_facecolor(color)

    for pc, color in zip(v2['bodies'],colors2):
        pc.set_facecolor(color)
    v1['cmeans'].set_colors(colors1)
    v1['cmaxes'].set_colors(colors1)
    v1['cmins'].set_colors(colors1)
    v1['cbars'].set_colors(colors1)
    v2['cmeans'].set_colors(colors2)
    v2['cmaxes'].set_colors(colors2)
    v2['cmins'].set_colors(colors2)
    v2['cbars'].set_colors(colors2)


    line = v1['cmaxes'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([1.7,1])+center
        new_line.append(l)
    v1['cmaxes'].set_segments(new_line)
    line = v1['cmeans'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([0.8,1])+center
        new_line.append(l)
    v1['cmeans'].set_segments(new_line)
    v1['cmeans'].set_visible(False)
    ax1.scatter(y=np.mean(np.array(data1), axis=1),x=ticks1, s=10, marker='o',c=colors1)


    
    line = v2['cmaxes'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([1.7,1])+center
        new_line.append(l)
    v2['cmaxes'].set_segments(new_line)
    line = v2['cmeans'].get_segments()
    new_line = []
    for l in line:
        center = l.mean(axis=0)
        l = (l-center)* np.array([0.8,1])+center
        new_line.append(l)
    v2['cmeans'].set_segments(new_line)
    v2['cmeans'].set_visible(False)
    ax2.scatter(y=np.mean(np.array(data2), axis=1),x=ticks2, s=10, marker='o',c=colors2)

    l_elems = [Patch(facecolor=cm1[0], label = 'AQE'),Patch(facecolor=cm1[3], label = 'RQE')]
    plt.legend(handles = l_elems,loc='upper right',fontsize = 14)
    title = "System:{}".format(sel)+"\n$E_{ED}$:"+" {}".format(sel_min)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    plt.title(titles[sel],fontsize=17)
    plt.tight_layout()
    plt.show()
    plt.savefig(sel+"_bigger.pdf")
    plt.close('all')



    """
    print(df_sel_normalized)
    print(len(df_sel))
    print(len(df_sel_normalized))

    fig, ax1 = plt.subplots(figsize=(10,10))
    pos1 = ax1.get_position()
    ax1=sns.violinplot(y="AQE", x="Config",width=0.2, hue = "Config", data=df_sel_normalized,inner=None, palette="Set2",cut=0)
    ax1.set_ylabel("AQE%")
    ax1.set_ylim(60,100.5)
    ax1.set_yticks([60,75,80,85,90,92.5,95,97.5,100])
    ax2 = ax1.twinx()
    ax2=sns.violinplot(y="RQE", x="Config",width=0.2, hue="Config", data=df_sel_normalized,inner=None, palette="Pastel2",cut=0)
    pos2=ax2.get_position()
    ax2.set_p
    print(pos1, pos2)
    ax2.set_ylabel("RQE%")
    ax2.set_ylim(60,100.5)
    ax2.set_yticks([60,75,80,85,90,92.5,95,97.5,100])
    c = plt.gca().get_children()
    c_1 = ax1.get_children()
    c_2 = ax2.get_children()
    if len(df_sel) == 200:
        plt.legend([c_1[0],c_1[1],c_1[2],c_1[3],c_2[0],c_2[1],c_2[2],c_2[3]],['4LCX','5LCX','QIDACX','QIDASO4','4LCX','5LCX','QIDACX','QIDASO4'],loc='lower right')
    else:
        plt.legend([c_1[0],c_1[1],c_1[2],c_1[3],c_1[4],c_2[0],c_2[1],c_2[2],c_2[3],c_2[4]],['4LCX','5LCX','6LCX','QIDACX','QIDASO4','4LCX','5LCX','6LCX','QIDACX','QIDASO4'],loc='lower right')
    plt.xlabel("Configurations")
    title = "System:{}".format(sel)+"\n$E_{ED}$:"+" {}".format(sel_min)
    plt.title(titles[sel])
    plt.show()
    plt.savefig(sel+".pdf")
    plt.close('all')
    """
    #SEABORN TEST
    """
    fig, ax1 = plt.subplots(figsize=(10,10))
    ax1.xaxis.grid(True)
    ax1.set_ylabel("AQE%")
    ax1.set_ylim(60,100.5)
    ax1.set_yticks([60,75,80,85,90,92.5,95,97.5,100])
    ax1.set_xticks([0,2,4,6,8])
    #ax1.set_xticks([-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5])
    ax2 = ax1.twinx()
    ax2.set_ylabel("RQE%")
    ax2.set_ylim(60,100.5)
    ax2.set_yticks([60,75,80,85,90,92.5,95,97.5,100])
    if len(df_sel) == 200:
        offset= 0
        ax1=sns.violinplot(y="AQE", x=np.repeat([0-offset,2-offset,4-offset,6-offset],50),width=0.5, hue = "Config", data=df_sel_normalized,inner=None, palette="Set2",cut=0)
        ax2=sns.violinplot(y="RQE", x=np.repeat([0+offset,2+offset,4+offset,6+offset],50),width=0.5, hue="Config", data=df_sel_normalized,inner=None, palette="Pastel2",cut=0)
    else:
        offset=0
        ax1=sns.violinplot(y="AQE", x=np.repeat([0-offset,2-offset,4-offset,6-offset,8-offset],50),width=0.5, hue = "Config", data=df_sel_normalized,inner=None, palette="Set2",cut=0)
        ax2=sns.violinplot(y="RQE", x=np.repeat([0+offset,2+offset,4+offset,6+offset,8+offset],50),width=0.5, hue="Config", data=df_sel_normalized,inner=None, palette="Pastel2",cut=0)

    c = plt.gca().get_children()
ax.set_xlim(0,50)
    plt.title(titles[sel])
    plt.show()
    plt.savefig(sel+".pdf")
    plt.close('all')
    """
    """
    fig, ax1 = plt.subplots(figsize=(10,10))
    ax1.xaxis.grid(True)
    ax1.set_ylabel("AQE%")
    ax1.set_ylim(60,100.5)
    ax1.set_yticks([60,75,80,85,90,92.5,95,97.5,100])
    ax2 = ax1.twinx()
    ax2.set_ylabel("RQE%")
    ax2.set_ylim(60,100.5)
    ax2.set_yticks([60,75,80,85,90,92.5,95,97.5,100])

    if len(df_sel)==200:
        df_sel_normalized['x_1']=np.repeat([0.1,1.1,2.1,3.1],50)
        df_sel_normalized['x_2']=np.repeat([-0.1,0.9,1.9,2.9],50)

        df_x1 = pd.DataFrame({'AQE':df_sel_normalized['AQE'].astype(float),'x_1': df_sel_normalized['x_1'].astype(float),'Config':df_sel_normalized['Config']})
        df_x2 = pd.DataFrame({'RQE':df_sel_normalized['RQE'].astype(float),'x_2':df_sel_normalized['x_2'].astype(float),'Config':df_sel_normalized['Config']})
        print(df_x1.columns)
        print(len(df_x1['AQE']),len(df_x1['x_1']), len(df_x2['RQE']), len(df_x2['x_2']))
        ax1.set_xticks([0.1,1.1,2.1,3.1])
        ax2.set_xticks([-0.1,0.9,1.9,2.9])
        ax1=sns.violinplot(y="AQE", x="x_1",width=0.5, hue = "Config", data=df_x1,inner=None, palette="Set2",cut=0)
        ax2=sns.violinplot(y="RQE", x="x_2",width=0.5, hue="Config", data=df_x2,inner=None, palette="Pastel2",cut=0)
    else:
        df_sel_normalized['x_1']=np.repeat([0.1,1.1,2.1,3.1,4.1],50)
        df_sel_normalized['x_2']=np.repeat([-0.1,0.9,1.9,2.9,3.9],50)

        print(len(df_x1['AQE']),len(df_x1['x_1']), len(df_x2['RQE']), len(df_x2['x_2']))
        ax1.set_xticks([0.1,1.1,2.1,3.1,4.1])
        ax2.set_xticks([-0.1,0.9,1.9,2.9,3.9])
        df_x1 = pd.DataFrame({'AQE':df_sel_normalized['AQE'].astype(float),'x_1': df_sel_normalized['x_1'].astype(float),'Config':df_sel_normalized['Config']})
        df_x2 = pd.DataFrame({'RQE':df_sel_normalized['RQE'].astype(float),'x_2':df_sel_normalized['x_2'].astype(float),'Config':df_sel_normalized['Config']})
        ax1=sns.violinplot(y="AQE", x="x_1",width=0.5, hue = "Config", data=df_x1,inner=None, palette="Set2",cut=0)
        ax2=sns.violinplot(y="RQE", x="x_2",width=0.5, hue="Config", data=df_x2,inner=None, palette="Pastel2",cut=0)

    c = plt.gca().get_children()
    c_1 = ax1.get_children()
    c_2 = ax2.get_children()
    if len(df_sel) == 200:
        plt.legend([c_1[0],c_1[1],c_1[2],c_1[3],c_2[0],c_2[1],c_2[2],c_2[3]],['4LCX','5LCX','QIDACX','QIDASO4','4LCX','5LCX','QIDACX','QIDASO4'],loc='lower right')
    else:
        plt.legend([c_1[0],c_1[1],c_1[2],c_1[3],c_1[4],c_2[0],c_2[1],c_2[2],c_2[3],c_2[4]],['4LCX','5LCX','6LCX','QIDACX','QIDASO4','4LCX','5LCX','6LCX','QIDACX','QIDASO4'],loc='lower right')
    plt.legend(loc='lower right')
    plt.xlabel("Configurations")
    title = "System:{}".format(sel)+"\n$E_{ED}$:"+" {}".format(sel_min)
    plt.title(titles[sel])
    plt.show()
    plt.savefig(sel+".pdf")
    plt.close('all')
    """
