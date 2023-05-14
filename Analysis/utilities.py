
# Standard modules
import argparse
import numpy as np
import pandas as pd
import scipy
import scipy.stats
from glob import glob
from pprint import pprint
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

# Custom modules
from visualise_biofilm import setFigSize,setClassStaticMembers
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm

import matplotlib as mpl

def setMPL():
    params={'font.family'           :'serif',
            'text.usetex'           : True,
            'axes.titlesize'        : 10,
            'axes.labelsize'        : 10,
            'xtick.labelsize'       : 10,
            'ytick.labelsize'       : 10,
            'legend.frameon'        : False,
            'lines.linewidth'       : 1.5,
            'legend.handlelength'   : 1.5,
            'legend.columnspacing'  : 0.75,
            'legend.handletextpad'  : 0.5,
            # 'legend.title_fontsize' : None,
            'legend.fontsize'       : 10,
            'font.size'             : 10,
            'legend.columnspacing'  : 1.0,
            "axes.formatter.useoffset":False,
            # "text.latex.preamble"   : [r"\usepackage[detect-all,locale=DE]{siunitx}",
            #                            r"\usepackage[T1]{fontenc}",
            #                            r"\usepackage{amsmath}",
            #                            r"\usepackage{physics}"]
            "text.latex.preamble" : "\n".join(
                 [r"\usepackage[detect-all,locale=DE]{siunitx}",
                  r"\usepackage[T1]{fontenc}",
                  r"\usepackage{amsmath}",
                  r"\usepackage{physics}"]
                  )
           }
    mpl.rcParams.update(params)

def getCells(file):
    """
        Parameters:
            file: str
                file to load cell data from

        Returns:
            cells: list of ChainingRodShapedBacterium
                cells from this data file
    """
    class_dict={ "ChainingRodShaped" : ChainingRodShapedBacterium }
    data = pd.read_csv(file,sep="\t")
    (cells,
     species_load_dict) = intialiseElementsFromData(data,class_dict)
    return cells

def setHyperParams(file):
    log_pattern=f"{os.path.split(file)[0]}/*log"
    log_file=glob(log_pattern)[0]
    hp=setClassStaticMembers(log_file)
    for key,val in hp.items():
        hp[key]=val[0]
    return hp

def getFigsize(cols=1):
    x,y=setFigSize(247)
    return cols*x,cols*y
