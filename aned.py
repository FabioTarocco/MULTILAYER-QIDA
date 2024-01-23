import numpy as np
import pandas as pd
from pprint import pprint


base_path = "./MULTIQIDA/VQE_res/"
add_path = ["3x3/{}_3x3_4LCX.pkl","3x3/{}_3x3_5LCX.pkl","3x3/{}_3x3_1ACX_1VBCX_1VLCX.pkl","3x3/{}_3x3_1ASO4_1BSO4_1LSO4.pkl",
            "2x6/{}_2x6_4LCX.pkl","2x6/{}_2x6_5LCX.pkl","2x6/{}_2x6_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","2x6/{}_2x6_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4/{}_3x4_4LCX.pkl","3x4/{}_3x4_5LCX.pkl","3x4/{}_3x4_6LCX.pkl","3x4/{}_3x4_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4/{}_3x4_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_mag/{}_3x4_mag_4LCX.pkl","3x4_mag/{}_3x4_mag_5LCX.pkl","3x4_mag/{}_3x4_mag_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_mag/{}_3x4_mag_1ASO4_1BSO4_1CSO4_1LSO4.pkl",
            "3x4_delta/{}_3x4_delta_4LCX.pkl","3x4_delta/{}_3x4_delta_5LCX.pkl","3x4_delta/{}_3x4_delta_6LCX.pkl","3x4_delta/{}_3x4_delta_1ACX_1VBCX_1VCCX_1VDCX_1VLCX.pkl","3x4_delta/{}_3x4_delta_1ASO4_1BSO4_1CSO4_1DSO4_1LSO4.pkl",
            "3x4_delta_0_1/{}_3x4_delta_0_1_4LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_5LCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ACX_1VBCX_1VCCX_1VLCX.pkl","3x4_delta_0_1/{}_3x4_delta_0_1_1ASO4_1BSO4_1CSO4_1LSO4.pkl"]
conf_layers = [1,1,3,3,1,1,5,5,1,1,1,5,5,1,1,4,4,1,1,1,5,5,1,1,4,4]
conf_ed = [1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6]
#TOLTO 3x4_6LCX un 3 in conf_ed ed un 1 in layers

NVQE = 50
systems = ["3x3","2x6","3x4","3x4_mag","3x4_delta","3x4_delta_0_1"]
e_ed = [-4.7493272585529,-6.60347247538703,-6.69168019351493,-9.50847255610575,-5.33875124065748,-4.272670379839620]
neel = [-3.0,-4.0,-4.25,-4.25,-4.25,-4.25]
systems_configs = []
zip(add_path, conf_layers, conf_ed)
for item in zip(add_path, conf_layers, conf_ed):
    systems_configs.append((base_path+item[0],item[1],e_ed[item[2]-1],neel[item[2]-1]))

#pprint(systems_configs)

out = ""


for c,l,e,n in systems_configs:
    energies = []
    aned = 0.0
    print("OPENING FILE {}".format(c))
    print("SYSTEM ED:{}".format(e))
    for i in range(1,NVQE+1):
        vqe = pd.read_pickle(c.format(i))[0]
        if l == 1:
            en_i = vqe["Energy"]
            #print(i,en_i)
            energies.append(en_i)
        if l == 2:
            en_i = vqe["Energy_12"]
            #print(i,en_i)
            energies.append(en_i)
        if l == 3:
            en_i = vqe["Energy_123"]
            #print(i,en_i)
            energies.append(en_i)
        if l == 4:
            en_i = vqe["Energy_1234"]
            #print(i,en_i)
            energies.append(en_i)
        if l == 5:
            en_i = vqe["Energy_12345"]
            #print(i,en_i)
            energies.append(en_i)
    print("Gettin Min Energy")
    e_min = np.min(energies)
    print("Absolute Normalized Energy Deviation:")
    for j in energies:
        aned += e_min - j
    o = "RESULTS for {}\n\tED: {}\tNeel:{}\n\tMean: {}\n\tStd: {}\n\tVar: {}\n\tMED(Mean Energy Deviation): {}\n\tANED (ED): {}\n\tANED (E_MIN): {}\n\tBest VQE res: {}\n\tAPQE: {}%\n\tRPQE: {}%\n\n".format(c.format(i),e,n, np.mean(energies), np.std(energies), np.var(energies), (aned/NVQE),(aned/e), (aned/e_min), e_min,(e_min/e)*100,((e_min-n)/(e -n))*100)
    print(o)
    out += o
    print("Index of mins VQE: {}".format(np.argwhere(np.array(energies)==e_min)))
with open("aned_recap.txt", 'w') as file:
    file.write(out)







