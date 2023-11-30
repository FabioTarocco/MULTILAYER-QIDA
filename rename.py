import os
from os import listdir
from os.path import isfile, join

p = './MULTIQIDA/3x3/'
onlyfiles = [f for f in listdir(p) if isfile(join(p, f))]



for x in onlyfiles:
    print(x, x.replace('SU4', 'C4'))
    os.rename(x, x.replace('SU4', 'C4'))