"""
    Script to visual a biofilm infected by a phage

    Example usage:
    python3 visualise_phage.py               \
        --sus_file data/test_bacteria_IO.dat \
        --inf_file data/test_infected_IO.dat \
        --phg_file data/test_phage_IO.dat    \
        --fig_file output/vis_phage.pdf
"""

# Add parent to path
import sys
sys.path.insert(0,'..')

# Standard libraries
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as linalg
import pandas as pd
from FindDefects import intialiseElementsFromData

# Custom modules
from FindDefects import intialiseElementsFromData
from visualise_biofilm import VisBiofilm
from visualise_biofilm import setFigSize
from phage_vis_class import Phage
from RodShapedBacteria import RodShapedBacterium
from RodShapedBacteria import setProjection       # choose plane to visualise

class VisInfectedBiofilm(VisBiofilm):
    """
        Extend concepts from visualise biofilm to handle multiple different
        species.
    """

    def __init__(self,file_name,species_load_dict):
        print(f"loading {file_name}")
        try:
            self.element_list, self.species_load_dict = intialiseElementsFromData(
                pd.read_csv(file_name,sep='\t'),
                species_load_dict
                )
        except Exception as e:
            print(e)
            raise

    def setElementColourScheme(self,colours=None):
        """
        Parameters:
            colours: list, matplotlib colour schemes

        Effect:
            sets element colours in order of list
        """
        if colours==None:
            # Do nothing - classes use their default colours
            pass
        else:
            # set element colours
            # fill in later
            pass
