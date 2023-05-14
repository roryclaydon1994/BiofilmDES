
# Standard modules
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Custom modules
import utilities as ut

# Third party modules
from tqdm import tqdm
from joblib import Parallel, delayed

# def getCellStats(file):
#     cells=ut.getCells(file)
#     n=len(cells)
#     l=np.mean([cc.length for cc in cells])
#     s=np.std([cc.length for cc in cells])
#     return n,l,s
#
# ns=[]
# avls=[]
# stds=[]
# pattern="../GeneratedOutput/SimOutput/ChainingBucklingTransition/run1/*/*.dat"
# files=glob(pattern,recursive=True)
# data=Parallel(n_jobs=12, verbose=10, backend="loky")(
#     delayed(getCellStats)(file) for file in files
#     )
# data=np.array(data)
# print(data)
# np.savetxt("../GeneratedOutput/data/average_length_cell.txt",data)

data=np.loadtxt("../GeneratedOutput/data/average_length_cell.txt")
al=np.mean(data[:,1])
print(f"{al=}")
df=pd.DataFrame(data,columns=['N','l','s'])
df=df.groupby(['N']).mean()
gdata=df.reset_index().to_numpy()

ut.setMPL()
fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
ax.scatter(gdata[:,0],gdata[:,1],fc='w',ec='k',lw=1)
ax.axhline(3.5,linestyle="--",color='g',linewidth=1.5)
ax.axhline(al,linestyle="--",color='m',linewidth=1.5)
ax.set_ylabel(r"$\langle \ell \rangle$")
ax.set_xlabel(r"$N$")
fig.savefig("../GeneratedOutput/SI/average_length_cell.pdf",
            transparent=True,dpi=300,bbox_inches='tight')
plt.show()
