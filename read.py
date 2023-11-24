import sys
import pandas as pd
import numpy as np
from termcolor import colored
from colorama import init
import os
from pprint import pprint


print("-------------------------------------------------------------------------------------------")
print("-------------------------------------------------------------------------------------------")
n_vqe = 50
Nx = 3
Ny = 4

depth = 1
mode = 'A'
gate = 'SU4'
base_vqe = "./MULTIQIDA/VQE_res/"
f = "{}x{}/"
name = "{}_{}x{}_{}{}{}"
energies = []
time = []
for i in range(1, n_vqe+1):
    f_vqe = base_vqe + f.format(Nx,Ny) + name.format(i,Nx,Ny,depth,mode,gate) + ".pkl"
    vqe = pd.read_pickle(f_vqe)[0]
    print(colored("ITERATION {}".format(i), 'light_magenta'))
    print("VQE")
    print("\tEnergy: {}".format(vqe["Energy"]))
    energies.append(vqe["Energy"])
    time.append(vqe["Time_req"])
    print("-------------------------------------------------------------------------------------------")
print("Average Energy obtained in {} VQEs:".format(n_vqe)+colored(" \nEnergy:{} Std:{} Variance: {} ".format(np.mean(energies), np.std(energies), np.var(energies)), 'green'))
print("Min energy : {}".format(np.min(energies)))
#print(np.mean(time)/3600.0)

