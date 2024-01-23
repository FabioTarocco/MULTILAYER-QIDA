import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fig = plt.figure()
ax=fig.add_subplot(projection='3d')
#ax.scatter(np.random.uniform(size=1000, low=0, high=50),np.random.uniform(size=1000),np.random.uniform(size=1000),c='r',alpha=0.3)
xm = np.array([15,16,17]).T
ym = np.array([.5,.1,.7]).T
zm = np.array([.5,.9,.7]).T
atoms = pd.DataFrame({'x':xm,'y':ym,'z':zm})
bonds= [(1,2),(2,3)]
lines = []
for a1,a2 in bonds:
    x1,y1,z1 = atoms.iloc[a1-1]
    x2,y2,z2 = atoms.iloc[a2-1]
    lines.append([[x1,x2],[y1,y2],[z1,z2]])
for l in lines:
    ax.plot(l[0],l[1],l[2],c='brown',linewidth=2)
ax.scatter(xm,ym,zm,s=100,alpha=1,c='brown')

xm = np.array([.5,4,17,.9]).T
ym = np.array([.5,.1,.7,.6]).T
zm = np.array([.5,.9,.2,.3]).T
atoms = pd.DataFrame({'x':xm,'y':ym,'z':zm})
bonds= [(1,2),(2,3),(1,4),(4,3)]
lines = []
for a1,a2 in bonds:
    x1,y1,z1 = atoms.iloc[a1-1]
    x2,y2,z2 = atoms.iloc[a2-1]
    lines.append([[x1,x2],[y1,y2],[z1,z2]])
for l in lines:
    ax.plot(l[0],l[1],l[2],c='blue',linewidth=2)
ax.scatter(xm,ym,zm,s=100,alpha=1,c='blue')
ax.set_xlim(0,50)
ax.set_ylim(0,1)
ax.set_zlim(0,1)
plt.show()