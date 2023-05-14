
# Standard modules
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Custom modules
import utilities as ut
from ChainingRodShapedBacteria import ChainingRodShapedBacterium

# Third party modules
from tqdm import tqdm
from joblib import Parallel, delayed

def getChainStats(file):
    cells=ut.getCells(file)
    hp=ut.setHyperParams(file)
    print(f"lp={hp['LinkingProb']}")
    _,cs=ChainingRodShapedBacterium.computeChainDescriptors(cells)
    # except Exception as e:
    #     print(e)
    #     print(file)
    #     quit()
    cs['LinkingProb']=hp['LinkingProb']
    cs['Ncells']=len(cells)
    return cs

try:
    ndf=pd.read_csv("../GeneratedOutput/data/chain_stats.txt")
except Exception as e:
    pattern="../GeneratedOutput/SimOutput/ChainingBucklingTransition/**/*/final*.dat"
    files=glob(pattern,recursive=True)
    # files=[file for file in files if "/run1/" not in file]
    print(f"processing {len(files)} files")

    # print(getChainStats(files[len(files)//2]))
    # quit()
    data=Parallel(n_jobs=12, verbose=10, backend="loky")(
        delayed(getChainStats)(file) for file in files
        )
    df=pd.DataFrame(data)
    print(df)

    file="../GeneratedOutput/data/chain_stats.txt"
    df.to_csv(file,index=False)

ndf=ndf.groupby(['LinkingProb']).mean()
gdata=ndf.reset_index()
gdata.sort_values(by=['LinkingProb'])
print(gdata)

ut.setMPL()

fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
# ax.scatter(gdata['LinkingProb'],gdata['Ncells']/gdata['Nchains'],fc='w',ec='k',lw=1)
ax.scatter(gdata['LinkingProb'],gdata['n_c'],fc='w',ec='k',lw=1)
ax.set_ylabel(r"$\langle n_c \rangle$")
ax.set_xlabel(r"$p_\ell$")
ax.plot(gdata['LinkingProb'],1/(1-gdata['LinkingProb']),'r--')
ax.set_yscale('log')
fig.savefig(f"../GeneratedOutput/SI/chain_stat_av_number_in_chains.pdf",
            transparent=True,dpi=300,bbox_inches='tight')
plt.show()

fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
ax.scatter(gdata['LinkingProb'],gdata['persistence_length'],fc='w',ec='k',lw=1)
ax.set_ylabel(rf"$\ell_p$")
ax.set_xlabel(r"$p_\ell$")
# ax.set_yscale('log')
fig.savefig(f"../GeneratedOutput/SI/persistence_length_chains.pdf",
            transparent=True,dpi=300,bbox_inches='tight')
plt.show()
