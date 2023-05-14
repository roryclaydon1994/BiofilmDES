"""
    Utilities for plotting the analysis
"""

# Standard libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import axis,cm
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=10)
import numpy as np
import numpy.linalg as linalg
import scipy
import pandas as pd
import re
import os
from pprint import pprint
from functools import reduce
from shapely.geometry import Polygon,Point
from shapely.affinity import scale,rotate
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection, LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from glob import glob
from functools import reduce
from scipy.signal import find_peaks
from copy import copy
from itertools import permutations

# third party modules
from tqdm import tqdm
from skimage.morphology import skeletonize, medial_axis
from skimage.measure import EllipseModel
try:
    from sklearn.cluster import KMeans,MiniBatchKMeans,Birch,DBSCAN
    from sklearn.mixture import GaussianMixture
except Exception as e:
    print(e)

from joblib import Parallel, delayed
import multiprocessing

# Custom modules
import generalPlotting as gp
import utilities as ut
from DistributionFunctions import (getRadialDistributionFunction,
                                   getColonyDensity,
                                   computeColonyContour,
                                   getClusterContoursDescriptors,
                                   boxCount,
                                   getContourCurvature,
                                   fitEllipse,
                                   checkBuckled,
                                   getAreaFraction,
                                   filterOutliers,
                                   getCurvatureGridMethod,
                                   getCorrelationLengths,
                                   getPairCorrelationData,
                                   getWeightedPairCorrelationData,
                                   getColonyWavelength
                                   )
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from visualise_biofilm import setFigSize
from fastCellPlotting import addAllCellsToPlot
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from SphericalBacteria import SphericalBacterium

ut.setMPL()

def getMeasurementData(x):
    """
        For a given input of independent val, list of measurements, extract
        the mean and sem
    """
    y=np.array(x)
    return y[0],np.mean(y[1:]),scipy.stats.sem(y[1:])

def saveClusters(clusters,filename):
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

    with open(filename,'w') as f:
        for key,val in clusters.items():
            out=f"{key}"
            for cell in val:
                out+=f" {cell.cell_id}"
            f.write(f"{out}\n")

def loadClusters(cells,filename):
    cluster_ids={}
    with open(filename,'r') as f:
        lines=f.readlines()
        # print(lines)
        for line in lines:
            data=line.split()
            key=int(data[0])
            cluster_ids[key]=[]
            [cluster_ids[key].append(int(data[ii])) for ii in range(1,len(data))]

    cell_clusters={}
    cd={}
    for cell in cells:
        cd[f"{cell.cell_id}"]=cell

    for key,val in cluster_ids.items():
        cell_clusters[f"{key}"]=[]
        for id in val:
            cell_clusters[f"{key}"].append(cd[f"{id}"])
    return cell_clusters

class AnalysisPlotter(object):
    """
        Control plots from all the analysis.
    """

    def __init__(self, in_files,
                       analysis_dir,
                       hps,
                       class_dict={
                       "RodShaped"         : RodShapedBacterium,
                       "Spherical"         : RodShapedBacterium,
                       "ChainingRodShaped" : ChainingRodShapedBacterium
                       }
                       ):
        """
            Parameters:
                in_files: list of str
                    Data files expected with column headers
                    cell_type,cell_id,length,radius,pos_x,pos_y,pos_z,ori_x,ori_y,ori_z
                    with '\t' delimeters
                analysis_dir: str
                    The analysis directory to save the outputs in
                hps: dict
                    dictionary of the hyperparamters
                class_dict: dict
                    These define how objects will appear on the axes. The list should
                    be a dictionary of cell_type keys with value the class to use to
                    represent this.
        """
        self.in_files=in_files
        self.analysis_dir=analysis_dir
        self.hps=hps
        self.class_dict=class_dict
        self.make_output_figures=False
        self.generate_data_only=True

    def getHP(self,att,run=None):
        try:
            lp=self.hps[run][att]
        except Exception as e:
            print(e)
            try:
                lp=self.hps[att]
            except Exception as e:
                print(e)
                quit()
        return lp

    def makeFileNameFromInfile(self,in_file,dir_name,mod_stem=None,ext="dat"):
        run,repeat,stem = self.getRunRepeatStem(in_file)
        dir=self.analysis_dir
        dir=f"{dir}/{run}/{repeat}/{dir_name}/"
        if not os.path.exists(dir):
            print(f"Creating {dir}")
            os.makedirs(dir)
        if mod_stem==None:
            filename=(
                dir+
                f"{stem}_{dir_name}s.{ext}"
                )
        else:
            filename=(
                dir+
                f"{stem}_{mod_stem}.{ext}"
                )
        return filename

    @staticmethod
    def getRunRepeatStem(in_file):
        run = re.search(r"run\d+",in_file)[0]
        try:
            repeat = re.search(r"repeat\d+",in_file)[0]
        except Exception as e:
            repeat = "repeat1"
        stem=in_file.split("/")[-1].replace(".dat","")
        return run, repeat, stem

    def filterFiles(self,lps,brs,one_repeat=True):
        use_files=[]
        _,use_repeat,_ = self.getRunRepeatStem(self.in_files[0])
        for file in self.in_files:
            run,repeat,stem = self.getRunRepeatStem(file)
            if repeat != use_repeat and one_repeat:
                continue
            if not np.any(self.hps[run]['BendRig']==brs):
                continue
            elif not np.any(self.hps[run]['LinkingProb']==lps):
                continue
            use_files.append(file)
        return use_files

    def makeAnalysisDirectory(self,in_file,dname):
        """
            Based on input filename, create a new directory to save related
            analysis data to
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)
        dir=self.analysis_dir
        dir+=f"{run}/{repeat}"
        dir+=f"/{dname}/"

        if not os.path.exists(dir):
            os.makedirs(dir)

        return dir

    def getColonyContour(self,in_file,cells_):
        """
            Collect the colony contour from file, or calclate it and save to file

            Parameters:
                in_file: str
                    name of the file used to generate the contour
                cells_: list of cells
                    cells to find the contour length of
            Returns:
                coords: (N,2) float numpy array
                    contour coordinates
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)
        contour_dir=self.analysis_dir
        contour_dir+=f"{run}/{repeat}"
        contour_dir+="/contour/"
        coord_filename=(contour_dir+
                        f"{stem}_coords.txt")
        try:
            coords=np.loadtxt(coord_filename)
            print(f"Loaded contour data from {coord_filename}")
        except Exception as e:
            coords = computeColonyContour(cells_,add_links=True)
            if not os.path.exists(contour_dir):
                os.makedirs(contour_dir)
            np.savetxt(coord_filename,coords)
            print(f"Saved contour data to {coord_filename}")
        return coords

    def getAspectRatio(self,in_file,cells_):
        """
            Calculate the aspect ratio of the colony by least squares fitting an
            ellipse to the colony contour

            Parameters:
                in_file: str
                    File the aspect ratio is to be generated for
                cells_: list of cells
                    cells to find the contour length of
            Returns:
                alpha: float
                    aspect ratio
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)
        coords=self.getColonyContour(in_file,cells_)
        dir=self.analysis_dir
        dir+=f"{run}/{repeat}"
        dir+="/contour/"
        ellipse_filename=(dir+
                         f"{stem}_ellipse.pdf")

        params=fitEllipse(coords)
        xc, yc, a, b, theta = params
        alpha=max(a/b,b/a)
        if self.make_output_figures:
            fig,ax = plt.subplots(1,1)
            ax.plot(coords[:,0],coords[:,1])
            ell_patch = Ellipse((xc, yc), 2*a, 2*b, theta*180/np.pi,
                                edgecolor='red', facecolor='none'
                                )
            ax.add_patch(ell_patch)
            ax.axis("scaled")
            ax.text(xc,yc,rf"$\alpha=${alpha:3.3f}")
            fig.savefig(ellipse_filename,format="pdf",bbox_inches='tight',
                        transparent=True)
            plt.close()
        # plt.show()
        # xc, yc, a, b, theta = ell.params
        return max(a/b,b/a)

    def getPersistenceLength(self,in_file,cells_):
        """
        Get the persistence length of chains
        """
        fname=self.makeFileNameFromInfile(in_file,'persistence_length')

        # try:
        #     persistence_length = np.loadtxt(fname)
        # except Exception as e:
        #     print(e)
        chains=ChainingRodShapedBacterium.getChains(cells_)
        persistence_length=ChainingRodShapedBacterium.computePersistenceLength(chains)
            # np.savetxt(fname,persistence_length)
        return persistence_length

    def getLinkAngle(self,in_file,cells_):
        """
        Get the average linking angle between cells
        """
        fname=self.makeFileNameFromInfile(in_file,'link_angle')

        try:
            thetas = np.loadtxt(fname)
        except Exception as e:
            print(e)
            chains=ChainingRodShapedBacterium.getChains(cells_)
            thetas=np.array([ theta for chain in chains
                              for theta in ChainingRodShapedBacterium.getCurvature(chain)
                              ])
            np.savetxt(fname,thetas)
        if len(thetas):
            return np.mean(thetas)
        else:
            return 0

    def getFractalDimension(self,in_file,cells_):
        """
            Calculate the fractal dimension of the colony contour by using box
            counting.

            Parameters:
                in_file: str
                    File the aspect ratio is to be generated for
                cells_: list of cells
                    cells to find the contour length of
            Returns:
                fracdim: float
                    fractal dimension of the colony contour
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)

        dir=self.analysis_dir
        dir+=f"{run}/{repeat}"
        dir+="/fractal_dimension/"

        try:
            fractal_dim = np.loadtxt(f"{dir}/{stem}_fractal_dim.txt")
        except Exception as e:
            print(e)
            if not os.path.exists(dir):
                os.makedirs(dir)

            coords=self.getColonyContour(in_file,cells_)
            # coords = computeColonyContour(cells_)
            fractal_dim,boxes = boxCount(coords)
            np.savetxt(f"{dir}/{stem}_fractal_dim.txt",[fractal_dim])
            if self.make_output_figures==True:
                fractal_filename=(
                    f"{dir}/{stem}_fractal_dim_{fractal_dim:3.3f}.pdf"
                    )
                AnalysisPlotter.plotBoxCount(fractal_dim,boxes,coords,
                                             fractal_filename,cells_)
        return fractal_dim

    def getColonyAreaFraction(self,in_file,cells_):
        return getAreaFraction(cells_)

    def getBuckled(self,in_file,cells_):
        """
            Roughly determine if the curve has buckled

            Parameters:
                in_file: str
                    File tto check if buckled
                cells_: list of cells
                    cells to find the contour length of
            Returns:
                bukled: bool
                    True if this colony buckled
        """
        # run,repeat,stem = self.getRunRepeatStem(in_file)
        #
        # dir=self
        # dir+=f"{run}/{repeat}"
        # dir+="/buckled/"
        #
        # if not os.path.exists(dir):
        #     os.makedirs(dir)

        buckled=checkBuckled(cells_)
        print(f"{os.path.split(in_file)[-1]} === {buckled=}")
        return buckled

    def getDensity(self,in_file,cells_):
        if len(cells_)==0:
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
        coords=self.getColonyContour(in_file,cells_)
        density,_=getColonyDensity(cells_,coords=coords)
        return density

    def getArea(self,in_file,cells_):
        if len(cells_)==0:
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
        coords=self.getColonyContour(in_file,cells_)
        _,gon=getColonyDensity(cells_,coords=coords)
        return gon.area

    def getGlobalOrder(self,in_file,cells_):
        Q_file=self.makeFileNameFromInfile(in_file,"Qcell",ext='npy')
        rs_file=self.makeFileNameFromInfile(in_file,dir_name="Qcell",
                                            mod_stem="rs",ext='npy')
        try:
            Qs=np.load(Q_file)
        except Exception as e:
            N=len(cells_)
            if N==0:
                data = pd.read_csv(in_file,sep="\t")
                (cells_,
                 species_load_dict) = intialiseElementsFromData(data,self.class_dict)
            RodShapedBacterium.getNeighbourLists(cells_,r=6,find_contacts=True)
            Qs=np.zeros((N,2,2))
            rs=np.zeros((N,2))
            for ii,cc in enumerate(cells_):
                L=0
                rs[ii]=cc.rcm[:-1]
                Qs[ii]+=RodShapedBacterium.getCellQ(cc)*cc.length
                L+=cc.length
                for nn in cc.neighbours:
                    Qs[ii]+=RodShapedBacterium.getCellQ(nn)*nn.length
                    L+=nn.length
                Qs[ii]/=L

            np.save(Q_file,Qs)
            np.save(rs_file,rs)

        S=np.sqrt(np.mean(Qs[:,0,0]**2)+np.mean(Qs[:,0,1]**2))

        if self.make_output_figures==True:
            run,repeat,stem = self.getRunRepeatStem(in_file)
            if repeat=="repeat1":
                fname=self.makeFileNameFromInfile(in_file,dir_name="Qcell",
                                                  mod_stem=f"pcolor_{S=:3.3f}",ext='pdf')
                print(f"make a picture {fname=}")

                rs=np.load(rs_file)
                Ss=np.sqrt(Qs[:,0,0]**2+Qs[:,0,1]**2)
                gp.plotUnstructuredMesh(cells_,Ss,rs,fname)
        return S

    def getPairCorrelation(self,in_file,cells_,diameter=1.0):
        """

        """
        cl_file=self.makeFileNameFromInfile(in_file,"pair_correlation",ext='npy')
        rs_file=self.makeFileNameFromInfile(in_file,"pair_correlation",
                                            mod_stem='rs',ext='npy')
        phis_rs_file=self.makeFileNameFromInfile(in_file,"pair_correlation",
                                                 mod_stem='phis_rs',ext='npy')
        print(cl_file)
        try:
            c_ls=np.load(cl_file)
            rs=np.load(rs_file)
            # phis=np.load(phis_rs_file)
            # rs=np.load(phis_rs_file)
            # print(f"loaded {rs.min()=}{rs.max()=}")
        except Exception as e:
            print(e)
            N=len(cells_)
            if N==0:
                data = pd.read_csv(in_file,sep="\t")
                (cells_,
                 species_load_dict) = intialiseElementsFromData(data,self.class_dict)
                N=len(cells_)

            rs,phis=getPairCorrelationData(cells_)
            # np.save(phis_rs_file,phis)
            # np.save(phis_rs_file,rs)

            edges=np.arange(0,rs.max()+diameter,diameter)
            idx=np.digitize(x=rs,bins=edges)

            # c_ls=np.zeros((len(edges)-1,2))
            # for ii in range(len(phis)):
            #     c_ls[idx[ii]-1,0]+=phis[ii]
            #     c_ls[idx[ii]-1,1]+=1
            #
            # c_ls=c_ls[:,0]/c_ls[:,1]
            c_ls=np.zeros((len(edges)-1))
            # errors=np.zeros_like(c_ls)
            for ii in range(1,len(edges)):
                sel=np.where(idx==ii)
                c_ls[ii-1]=np.mean(phis[sel])
                # errors[ii-1]=np.std(phis[sel])
            rs=edges[:-1]+0.5*diameter
            # plt.plot(edges[:-1],c_ls)
            # plt.show()
            np.save(cl_file,c_ls)
            np.save(rs_file,rs)


            """
            Note: relegated method but can be better optimised and gives nicer
            statistics
            """

        peaks,_=find_peaks(-c_ls)
        # print(peaks)
        scal=lambda x,xi: np.exp(-x/xi)
        popt,pcov=scipy.optimize.curve_fit(
            scal,rs[:peaks[0]],c_ls[:peaks[0]],
            p0=(1)
            )
        # print(f"{popt=}")
        # print(f"{np.sqrt(np.diag(pcov))=}")
        # fig,ax=plt.subplots()
        # ax.plot(rs,c_ls)
        # ax.plot(rs,scal(rs,*popt),'--')
        # ax.set_xlabel("$r$")
        # ax.set_ylabel("$c(r)$")
        # ax.axvline(rs[peaks[0]])
        # fig_file=self.makeFileNameFromInfile(in_file,"pair_correlation",
        #                                      mod_stem='example',ext='pdf')
        # fig.savefig(fig_file,transparent=True,bbox_inches='tight')
        # plt.close()
        # quit()

        return popt[0]

    def getWeightedPairCorrelation(self,in_file,cells_,diameter=1.0):
        """

        """
        cl_file=self.makeFileNameFromInfile(in_file,"weighted_pair_correlation",ext='npy')
        rs_file=self.makeFileNameFromInfile(in_file,"weighted_pair_correlation",
                                            mod_stem='rs',ext='npy')
        print(cl_file)
        try:
            c_ls=np.load(cl_file)
            rs=np.load(rs_file)
            # print(f"loaded {rs.min()=}{rs.max()=}")
        except Exception as e:
            print(e)
            N=len(cells_)
            if N==0:
                data = pd.read_csv(in_file,sep="\t")
                (cells_,
                 species_load_dict) = intialiseElementsFromData(data,self.class_dict)
                N=len(cells_)

            rs,c_ls=getWeightedPairCorrelationData(cells_,diameter)
            np.save(cl_file,c_ls)
            np.save(rs_file,rs)

        scal=lambda x,xi: np.exp(-x/xi)
        popt,pcov=scipy.optimize.curve_fit(scal,rs,c_ls,p0=(1))

        # fig,ax=plt.subplots()
        # ax.plot(rs,c_ls)
        # ax.plot(rs,scal(rs,*popt),'--')
        # ax.set_xlabel("$r$")
        # ax.set_ylabel("$c(r)$")
        # fig_file=self.makeFileNameFromInfile(in_file,"weighted_pair_correlation",
        #                                      mod_stem='example',ext='pdf')
        # fig.savefig(fig_file,transparent=True,bbox_inches='tight')
        # plt.close()
        # quit()
        return popt[0]

    def getMeanClusterArea(self,in_file,cells_,thresh_th=3):
        clusters, cluster_name_stub = self.getCluster(in_file,cells_,
                                                      thresh_th=thresh_th)
        cluster_contours = AnalysisPlotter.getClusterContours(clusters)
        colony_contour=self.getColonyContour(in_file,cells_)

        fig_name=None
        if self.make_output_figures==True:
            run,repeat,stem = self.getRunRepeatStem(in_file)
            if repeat=="repeat1":
                fig_name=self.makeFileNameFromInfile(in_file,
                                                     dir_name=f'mean_cluster_{thresh_th}',
                                                     ext='pdf')
        mean_cluster_area=getClusterContoursDescriptors(
            clusters,
            cluster_contours,
            colony_contour,
            fig_name
            )
        return mean_cluster_area

    def getColonyCorrelations(self,in_file,cells_):
        rho_file=self.makeFileNameFromInfile(in_file,"rho",ext='npy')
        Q_file=self.makeFileNameFromInfile(in_file,"Q",ext='npy')
        C_file=self.makeFileNameFromInfile(in_file,"Correlation",ext='npy')
        return getCorrelationLengths(cells=cells_,dr=0.5,
                                     rho_file=rho_file,Q_file=Q_file,
                                     C_file=C_file)

    def getWavelength(self,in_file,cells_,nbins=10):
        return getColonyWavelength(cells_,nbins=nbins)[0]

    def getRoughness(self,in_file,cells_):
        if len(cells_)==0:
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
        ell = EllipseModel()
        coords=self.getColonyContour(in_file,cells_)
        ell.estimate(coords)
        res=ell.residuals(coords)
        return np.mean(res)

    def getCurvatures(self,in_file,cells_):
        """
            Calculate the curvature of the colony contour

            Parameters:
                in_file: str
                    File the aspect ratio is to be generated for
                cells_: list of cells
                    cells to find the contour length of
            Returns:
                max_curvature: float
                    maximum curvature of the colony contour
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)

        dir=self.analysis_dir
        dir+=f"{run}/{repeat}"
        dir+="/curvature/"

        if not os.path.exists(dir):
            os.makedirs(dir)

        fname=f"{dir}/rs_{stem}.npy"
        coords=self.getColonyContour(in_file,cells_)
        # L1,L2,rs,skel_rs=getContourCurvature(cells_,coords,fname)
        # curvs=np.array([L1,L2])
        # curvature=curvs.mean()
        # if abs(L2-L1)>0.2*np.mean([L1,L2])+1.0:
        #     print("error, curvatures are very different")
        #     print(f"{abs(L2-L1)=} {0.1*np.mean([L1,L2])+0.5=}")
        #     quit()
        L1=getCurvatureGridMethod(cells_,coords,fname)
        L2=getContourCurvature(cells_,coords,fname)
        curvature=0.5*(L1+L2)
        # np.savetxt(f"{dir}/{stem}_curv.txt",curvature)
        curvature_fname=(
            f"{dir}/{stem}_curv_{curvature:3.3f}.pdf"
            )
        if not os.path.exists(curvature_fname) and self.make_output_figures:
            rs=np.load(fname)
            skel_rs=np.load(f"{dir}/skel_{stem}.npy")
            AnalysisPlotter.plotCurvature(L1,L2,rs,skel_rs,coords,cells_,
                                          curvature_fname)
        return L1,L2

    def plotAspectRatio(self,var,show=False):
        """
            Plot the aspect ratio of the colony. This is defined by the least
            squares fitting of an ellipse to the colony

            Parameters:
                var: string
                    Name of the independent variable to plot against
                show: bool
                    Flag to switch on displaying figures generated from analysis
        """
        asps={}
        ivs={}
        # size=[]
        for in_file in self.in_files:
            run,repeat,stem = self.getRunRepeatStem(in_file)
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
            # size.append(self.getColonySize(in_file,cells_))
            asp = self.getAspectRatio(in_file,cells_)
            iv = self.getHP(var['iv'],run)
            if run not in ivs:
                if iv not in ivs.values():
                    ivs[run] = iv
                else:
                    print("This iv has already been processed")
                    quit()
            else:
                assert ivs[run]==iv, "error, ivs are different for the same run"

            if run not in asps:
                asps[run]=[ivs[run]]
            asps[run].append(asp)

        # size=np.mean(size)
        asps=np.array([getMeasurementData(data) for data in asps.values()])
        asps_data_name=f"{self.analysis_dir}/aspect_ratio.txt"
        np.savetxt(asps_data_name,asps)

        if self.generate_data_only==False:
            fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
            marker="o"
            ax.errorbar(asps[:,0],asps[:,1],asps[:,2],marker=marker,mfc="w",mec='k',
                        capsize=5,alpha=0.7)
            asps_save_name=f"{self.analysis_dir}/aspect_ratio.pdf"
            ax.set_ylabel("Aspect ratio")
            # ax.set_xlabel("Linking Probability")
            ax.set_xlabel(var['label'])
            fig.savefig(asps_save_name,format="pdf",bbox_inches='tight',
                        transparent=True)

            if show:
                plt.show()

            plt.close()

        return asps

    def plotFractalDimension(self,var,show=False):
        """
            Plot the fractal dimension of the colony contour. This uses basic
            box counting.

            Parameters:
                var : dict
                    iv: variable name, label: display label
                show: bool
                    Flag to switch on displaying figures generated from analysis
        """
        fds={} # fractal dimensions
        ivs={} # independent variables
        # size=[]
        for in_file in self.in_files:
            run,repeat,stem = self.getRunRepeatStem(in_file)
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
            # size.append(self.getColonySize(in_file,cells_))
            fd = self.getFractalDimension(in_file,cells_)

            iv = self.getHP(var['iv'],run)
            if run not in ivs:
                if iv not in ivs.values():
                    ivs[run] = iv
                else:
                    print("This iv has already been processed")
                    quit()
            else:
                assert ivs[run]==iv, "error, lps are different for the same run"

            if run not in fds:
                fds[run]=[ivs[run]]
            fds[run].append(fd)

        # size=np.mean(size)
        fds=np.array([getMeasurementData(data) for data in fds.values()])
        fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
        marker="o"
        ax.errorbar(fds[:,0],fds[:,1],fds[:,2],marker=marker,mfc="w",mec='k',
                    capsize=5,alpha=0.7)
        fds_save_name=f"{self.analysis_dir}/fractal_dimensions/fractal_dimension.pdf"
        ax.set_ylabel("$\mathcal{D}$")
        ax.set_xlabel(var['label'])
        fig.savefig(fds_save_name,format="pdf",bbox_inches='tight',
                    transparent=True)

        if show:
            plt.show()

        plt.close()

        return fds

    def getCluster(self,in_file,cells,thresh_th=3):
        """
            Collect the colony contour from file, or calclate it and save to file

            Parameters:
                in_file: str
                    name of the file used to generate the contour
                cells_: list of cells
                    cells to find the contour length of
                thresh_th: float
                    angle between cells to allow them to be in the same cluster
            Returns:
                clusters: dict
                    cluster index are keys, list of cells are values
        """
        run,repeat,stem = self.getRunRepeatStem(in_file)
        cluster_dir=self.analysis_dir
        cluster_dir+=f"{run}/{repeat}"
        cluster_dir+=f"/cluster{thresh_th}/"
        cluster_filename=(
            cluster_dir+
            f"{stem}_clusters.dat"
            )
        try:
            clusters=loadClusters(cells,cluster_filename)
            for key,cluster in clusters.items():
                cc=RodShapedBacterium.colour_fn(cluster[0].theta)
                for cell in cluster:
                    cell.colour=cc
        except Exception as e:
            print(e)
            if not os.path.exists(cluster_dir):
                os.makedirs(cluster_dir)
            clusters=RodShapedBacterium.findMicrodomains(
                cells,
                RodShapedBacterium.colour_fn,
                threshold_angle=thresh_th
                )
            saveClusters(clusters,cluster_filename)
        return clusters, f"{cluster_dir}/{stem}"

    @staticmethod
    def getClusterContours(clusters):
        cluster_contours={}
        for key,cells in tqdm(clusters.items()):
            coords = computeColonyContour(cells,add_links=True)
            cluster_contours[key]=coords
        return cluster_contours

    def plotClusterContours(self,clusters,cluster_contours,save_name):
        """
        :)
        """
        if self.make_output_figures:
            fig,ax = plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
            ax_lim=[np.inf,-np.inf]
            for key,coords in cluster_contours.items():
                ax_lim[0]=min(np.min(coords),ax_lim[0])
                ax_lim[1]=max(np.max(coords),ax_lim[1])
                p = Polygon([(v[0],v[1]) for v in coords])
                cluster = PolygonPatch(
                    p,
                    fc=clusters[key][0].colour,
                    ec='k',
                    lw=1/max(ax_lim),
                    alpha=1,
                    zorder=2
                    )
                ax.add_patch(cluster)
            ax.set_xlim(ax_lim)
            ax.set_ylim(ax_lim)
            ax.axis("scaled")
            fig.savefig(save_name,
                        format="pdf",bbox_inches='tight',
                        transparent=True)
            plt.close()

    def plotClusterAnalysis(self,var,show=False):
        cluster_analysis={}
        lps={}
        size=[]
        for in_file in self.in_files:
            print(f"{in_file=}")
            run,repeat,stem = self.getRunRepeatStem(in_file)
            ca_data_name=self.makeFileNameFromInfile(in_file,"cluster_mean")
            try:
                tca=np.loadtxt(ca_data_name)
                lp=tca[0]
                mean_cluster_area=tca[1]
                size.append(tca[2])
            except Exception as e:

                lp = self.getHP(var['iv'],run)

                data = pd.read_csv(in_file,sep="\t")
                (cells_,
                 species_load_dict) = intialiseElementsFromData(data,self.class_dict)
                tsize=self.getColonySize(in_file,cells_)
                size.append(tsize)
                clusters, cluster_name_stub = self.getCluster(in_file,cells_)

                # find the contours around each cluster
                cluster_contours = AnalysisPlotter.getClusterContours(clusters)

                # plot the cluster contours
                if self.make_output_figures:
                    self.plotClusterContours(
                        clusters,
                        cluster_contours,
                        f"{cluster_name_stub}_contour_cluster.pdf"
                        )

                # get the area of each contour
                colony_contour=self.getColonyContour(in_file,cells_)
                mean_cluster_area=getClusterContoursDescriptors(
                    clusters,
                    cluster_contours,
                    colony_contour
                    )

                np.savetxt(ca_data_name,([lp],[mean_cluster_area],[tsize]))

            if run not in lps:
                if lp not in lps.values():
                    lps[run] = lp
                else:
                    print("This lp has already been processed")
                    quit()
            else:
                assert lps[run]==lp, "error, lps are different for the same run"

            if run not in cluster_analysis:
                cluster_analysis[run]=[lps[run]]
            cluster_analysis[run].append(mean_cluster_area)

        size=np.mean(size)
        ca=np.array([getMeasurementData(data) for data in cluster_analysis.values()])
        if self.generate_data_only==False:
            fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
            marker="o"
            ax.errorbar(ca[:,0],ca[:,1],ca[:,2],marker=marker,mfc="w",mec='k',
                        capsize=5,alpha=0.7)
            ca_save_name=f"{self.analysis_dir}/{clusters}/mean_cluster_area_{size:3.3f}.pdf"
            ax.set_ylabel(r"Average domain area $\langle A\rangle$")
            ax.set_xlabel(var['label'])
            fig.savefig(ca_save_name,format="pdf",bbox_inches='tight',
                        transparent=True)

            if show:
                plt.show()

            plt.close()

        return ca

    def checkContours(self,show=False):
        """

        """
        for in_file in self.in_files:
            run, repeat, stem = self.getRunRepeatStem(in_file)
            contour_dir=self.analysis_dir
            contour_dir+=f"{run}/{repeat}"
            save_name=f"{contour_dir}/{stem}_check_contour.pdf"
            print(save_name)
            if os.path.exists(save_name):
                print(f"{save_name} already exists, skipping...")
                continue
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,
                                                            self.class_dict)
            coords = self.getColonyContour(in_file,cells_)
            density,contour_polygon=getColonyDensity(cells_,coords=coords)

            print(f"there are {len(cells_)} cells")

            fig,ax = plt.subplots(1,1,figsize=setFigSize(3*247))
            for cell in tqdm(cells_):
                cell.addElementToPlot(ax)

            ax.plot(coords[:,0],coords[:,1],'k',label="contour")
            ax.axis("scaled")
            ax.set_xlabel("$x$")
            ax.set_ylabel("$y$")
            plt.legend()
            # print(f"Saving contour fig to {save_name}")
            fig.savefig(save_name,
                        format="pdf",bbox_inches='tight',
                        transparent=True)
            if show:
                plt.show()

            plt.close()

    def plotRadialDist(self,show=True):
        """
            Plot radial density function

            Parameters:
                show: bool
                    Flag to switch on displaying figures generated from analysis
        """
        print("Update to something useful")
        quit()
        fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
        for in_file in self.in_files:
            data = pd.read_csv(in_file,sep="\t")
            cells_,species_load_dict = intialiseElementsFromData(data,self.class_dict)
            rs,probs=getRadialDistributionFunction(cells_)
            run=re.search(r"run\d+",in_file)[0]
            lp=self.getHP("LinkingProb",run)
            label=f"{lp:3.2f}"
            ax.plot(rs,probs,label=label)

        save_name=f"{self.analysis_dir}/radial_dist.pdf"
        ax.set_ylabel("Radial Probability")
        ax.set_xlabel("$r$")
        plt.legend(ncol=2,title="linking prob")
        fig.savefig(save_name,format="pdf",bbox_inches='tight',
                    transparent=True)
        if show:
            plt.show()

    def plotChainDescriptors(self,show=True):
        """
            Plot persistence length

            Parameters:
                show: bool
                    Flag to switch on displaying figures generated from analysis
        """
        p_lens=[]
        lps=[]
        for in_file in self.in_files:
            run, repeat = self.getRunRepeatStem(in_file)
            data = pd.read_csv(in_file,sep="\t")
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
            (chains,
             chain_stats
             ) = ChainingRodShapedBacterium.computeChainDescriptors(cells_)
            pprint(chain_stats)
            lp = self.getHP("LinkingProb",run)
            lps.append(lp)
            p_lens.append(chain_stats["persistence_length"])

        fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
        marker="o"
        ax.plot(lps,p_lens,marker=marker,mfc="w",mec='k')
        p_lens_save_name=f"{self.analysis_dir}/persistence_length.pdf"
        ax.set_ylabel(r"$\langle p_\ell \rangle$")
        ax.set_xlabel("Linking Probability")
        fig.savefig(p_lens_save_name,format="pdf",bbox_inches='tight',
                    transparent=True)

        if show:
            plt.show()

        plt.close()

    def plotDensity(self,var,show=True):
        """
            Plot densities

            Parameters:
                show: bool
                    Flag to switch on displaying figures generated from analysis

            Returns:
                densities: (N,M) np array
                    first column is the
        """
        densities={}
        contour_lengths={}
        lps={}
        size=[]
        for in_file in self.in_files:
            try:
                run,repeat,stem = self.getRunRepeatStem(in_file)
                data = pd.read_csv(in_file,sep="\t")
            except Exception as e:
                print(e)
                continue
            (cells_,
             species_load_dict) = intialiseElementsFromData(data,self.class_dict)
            size.append(self.getColonySize(in_file,cells_))
            coords=self.getColonyContour(in_file,cells_)
            density,contour_polygon=getColonyDensity(cells_,coords=coords)
            lp=self.getHP(var['iv'],run)
            if run not in lps:
                lps[run] = lp
            else:
                assert lps[run]==lp, "error, lps are different for the same run"

            for dd in [densities,contour_lengths]:
                if run not in dd:
                    dd[run]=[lps[run]]
            densities[run].append(density)
            contour_lengths[run].append(contour_polygon.length)

        size=np.mean(size)
        densities=np.array([getMeasurementData(data) for data in densities.values()])
        contour_lengths=np.array([getMeasurementData(data) for data in contour_lengths.values()])
        np.savetxt(f"{self.analysis_dir}/{density}/density_at_{size:3.3f}.txt",
                   densities
                   )
        np.savetxt(f"{self.analysis_dir}/{contour}/contour_lengths_at_{size:3.3f}.txt",
                   contour_lengths
                   )

        if self.generate_data_only==False:
            fig_dens,ax_dens = plt.subplots(1,1,figsize=setFigSize(247))
            fig_lens,ax_lens = plt.subplots(1,1,figsize=setFigSize(247))
            marker="o"
            ax_dens.errorbar(densities[:,0],densities[:,1],densities[:,2],
                             marker=marker,mfc="w",mec='k',
                             capsize=5,alpha=0.7)
            ax_lens.errorbar(contour_lengths[:,0],contour_lengths[:,1],
                             contour_lengths[:,2],
                             marker=marker,mfc="w",mec='k',
                             capsize=5,alpha=0.7)

            dens_save_name=f"{self.analysis_dir}/{density}/density_at_{size:3.3f}.pdf"
            lens_save_name=f"{self.analysis_dir}/{contour}/contour_lengths_at_{size:3.3f}.pdf"
            ax_dens.set_ylabel("Colony Density")
            ax_lens.set_ylabel("Colony Contour Length")

            for ax in [ax_dens,ax_lens]:
                ax.set_xlabel(var['label'])

            for fig,name in zip([fig_dens,fig_lens],[dens_save_name,lens_save_name]):
                fig.savefig(name,format="pdf",bbox_inches='tight',transparent=True)

            if show:
                plt.show()

            plt.close()

        return densities

    def plotGeneral(self,fn,var,ylabel,name,show=False):
        """
            Plot some attribute of the colony, extracted using fn

            Parameters:
                fn: function(in_file,cells_)
                    returns data for this file
                var: string
                    Name of the independent variable to plot against
                ylabel: string
                    The y axis label
                name: string or list of strings
                    name of the dependent variable
                show: bool
                    Flag to switch on displaying figures generated from analysis
        """
        # asps_save_data=f"{self.analysis_dir}/{name}/{name}.txt"
        # os.makedirs(os.path.dirname(asps_save_data),exist_ok=True)
        # try:
        #     # asps=np.loadtxt(asps_save_data)
        #     all_data=pd.read_csv(asps_save_data)
        #     asps=all_data.loc[:,[var['iv'],name,f"{name}_error"]].to_numpy()
        # except Exception as e:
        #     asps={}
        #     ivs={}
        #     # size=[]
        #     for in_file in self.in_files:
        #         run,repeat,stem = self.getRunRepeatStem(in_file)
        #         data = pd.read_csv(in_file,sep="\t")
        #         (cells_,
        #          species_load_dict) = intialiseElementsFromData(data,self.class_dict)
        #         # size.append(self.getColonySize(in_file,cells_))
        #         asp = fn(in_file,cells_)
        #         iv = self.getHP(var['iv'],run)
        #         if run not in ivs:
        #             if iv not in ivs.values():
        #                 ivs[run] = iv
        #             else:
        #                 print("This iv has already been processed")
        #                 quit()
        #         else:
        #             assert ivs[run]==iv, "error, ivs are different for the same run"
        #
        #         if run not in asps:
        #             asps[run]=[ivs[run]]
        #         asps[run].append(asp)
        #
        #     if len(self.in_files)==1:
        #         return
        #
        #     # size=np.mean(size)
        #     keys=[kk for kk in asps.keys()]
        #     asps=np.array([getMeasurementData(data) for data in asps.values()])
        #
        #
        #     headers=[key for key in self.hps[keys[0]].keys() if key!='SEED' ]
        #
        #     headers.extend([name,f'{name}_error'])
        #     all_data=[]
        #     for ii,key in enumerate(keys):
        #         all_data.append([val for key,val in self.hps[key].items() if key!='SEED'])
        #         all_data[-1].extend([*asps[ii,1:]])
        #     all_data=np.array(all_data)
        #     df=pd.DataFrame(all_data,columns=headers)
        #     df.to_csv(asps_save_data,index=False)
        print("Cores set to 12")
        df=self.getGeneralData(fn_dict={name:fn},cores=12)
        avg=df[[var['iv'],name]].copy()
        avg=avg.groupby(var['iv']).agg({name:['mean',['sem',scipy.stats.sem]]})
        avg.sort_values(by=var['iv'],inplace=True)
        # avg=avg.T.reset_index(drop=True).T.reset_index()
        avg.columns = ['_'.join(col).rstrip('_') for col in avg.columns.values]
        avg.reset_index(inplace=True)
        asps=avg.to_numpy()

        if self.generate_data_only==False:
            fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
            marker="o"
            ax.errorbar(asps[:,0],asps[:,1],asps[:,2],marker=marker,mfc="w",mec='k',
                        capsize=5,alpha=0.7)
            asps_save_name=f"{self.analysis_dir}/{name}/{name}.pdf"
            if not os.path.exists(os.path.dirname(asps_save_name)):
                os.makedirs(os.path.dirname(asps_save_name),exist_ok=True)
            ax.set_ylabel(ylabel)
            # ax.set_xlabel("Linking Probability")
            ax.set_xlabel(var['label'])
            fig.savefig(asps_save_name,format="pdf",bbox_inches='tight',
                        transparent=True)

            if show:
                plt.show()

            plt.close()

        return asps

    def getGeneralData(self,fn_dict,cores=None):
        """
            Plot some attribute of the colony, extracted using fn

            Parameters:
                fn_dict: dictionary
                    keys: names of variables
                    items: fns with inputs (in_file,cells_)
                cores: int
                    number of cpus to run on
        """
        fn_dict['size']=self.getColonySize
        name=list(fn_dict.keys())

        # add=lambda t,s: s+t
        # fname_stem=reduce(add,[f'{v}_' for v in name])[:-1]
        comp_dir=f"{self.analysis_dir}/all_data"
        os.makedirs(comp_dir,exist_ok=True)
        asps_save_data=f"{comp_dir}/general_all_data.txt"
        print(asps_save_data)
        if os.path.exists(asps_save_data):
            df_original=pd.read_csv(asps_save_data)
            name=[ nn for nn in name if nn not in df_original.columns ]

        fns = { nn:fn_dict[nn] for nn in name }
        print("Running for")
        pprint(fns)
        if name:
            run,repeat,stem = self.getRunRepeatStem(self.in_files[0])
            headers=[ key for key in self.hps[run].keys()
                      if key in name or
                         key not in ['SEED','dt','SeedHost','SeedDevice']
                         ]
            headers.extend(['run','repeat','stem'])
            headers.extend(name)

            def getParData(in_file):
                run,repeat,stem = self.getRunRepeatStem(in_file)
                data = pd.read_csv(in_file,sep="\t")
                (cells_,
                 species_load_dict) = intialiseElementsFromData(data,self.class_dict)
                asp=[]
                for key,f in fns.items():
                    print(f"Running for {key}")
                    try:
                        asp.append(f(in_file,cells_))
                    except Exception as e:
                        print(e)
                        quit()

                # asp = [f(in_file,cells_) for f in fn]
                tmp_all_data=[]
                tmp_all_data.extend([val for key,val in self.hps[run].items()
                                if key in headers])
                tmp_all_data.extend([run,repeat,stem])

                tmp_all_data.extend(asp)
                return tmp_all_data

            if cores==None:
                cores=multiprocessing.cpu_count()
            print(f"Running on {cores}/{multiprocessing.cpu_count()} cores")
            all_data=Parallel(n_jobs=cores,verbose=100)(
                delayed(getParData)(file) for file in self.in_files
                )
            all_data=np.array(all_data)
            df_new=pd.DataFrame(all_data,columns=headers)
            print(df_new)
            print("new columns")
            print(df_new.columns)
            try:
                for x in df_original.columns:
                    try:
                        df_new[x]=df_new[x].astype(df_original[x].dtypes.name)
                    except Exception as e:
                        pass
                df=pd.merge(
                    df_original,
                    df_new,
                    on=[ col for col in df_new.columns if col not in name ]
                    )
            except Exception as e:
                print(e)
                df=df_new
            df.to_csv(asps_save_data,index=False)
        else:
            df=df_original
        return df

    def phasePlot(self,fn,name,vars,lbl=None):
        if lbl==None:
            lbl=name

        data=self.getGeneralData(fn_dict={name:fn})
        nx=data[vars[0]['iv']].nunique()
        ny=data[vars[1]['iv']].nunique()
        data.sort_values(by=[ v['iv'] for v in (vars) ],inplace=True)

        avg=data.groupby('run').mean()
        avg.sort_values(by=[ v['iv'] for v in (vars) ],inplace=True)

        x=avg[vars[0]['iv']]
        y=avg[vars[1]['iv']]
        c=avg[name]

        c=c.to_numpy().reshape((nx,ny))
        xx=x.to_numpy().reshape((nx,ny))
        yy=y.to_numpy().reshape((nx,ny))
        fig,ax=plt.subplots(1,1,figsize=setFigSize(247))

        print(f"{c.min()} {c.max()}")
        if c.min()<=0:
            c=np.abs(c)
            CS=ax.pcolor(xx,yy,c,norm=mpl.colors.LogNorm(vmin=c.min(),vmax=c.max()),
                         shading='auto'
                         )
        else:
            CS=ax.pcolor(xx,yy,c,norm=mpl.colors.LogNorm(vmin=c.min(),vmax=c.max()),
                         shading='auto'
                         )
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar=fig.colorbar(CS,cax=cax)
        cbar.set_label(lbl,size=14)
        cbar.ax.tick_params(labelsize=12)
        ax.set_xlabel(vars[0]['label'])
        ax.set_ylabel(vars[1]['label'])
        ax.set_yscale('log')
        comp_dir=f"{self.analysis_dir}/{name}"
        os.makedirs(comp_dir,exist_ok=True)
        fig.savefig(f"{comp_dir}/{name}_{vars[0]['label']}_{vars[1]['label']}_phase_pcolor.pdf",
                    format="pdf",
                    bbox_inches='tight',
                    transparent=True)

    def makeCollagePlot(self,variables,full=True):
        """
        Make a scatter plot of all the param runs
        """
        brs=[]
        lps=[]
        for ii,infile in enumerate(self.in_files):
            run,repeat,stem = self.getRunRepeatStem(infile)
            brs.append(self.hps[run]['BendRig'])
            lps.append(self.hps[run]['LinkingProb'])
        brs=np.sort(np.unique(np.array(brs)))
        lps=np.sort(np.unique(np.array(lps)))

        brs=brs[brs>0.5]
        nlp=len(lps)
        nbr=len(brs)
        if nlp*nbr>40:
            brs=brs[brs>0.5]
            brs=brs[::2]
            lps=lps[lps>=0.2]
            lps=lps[::2]
            nlp=len(lps)
            nbr=len(brs)
            # while nlp*nbr>20:
            #     lps=lps[::2]
            #     nlp=len(lps)
            #     nbr=len(brs)
        if nlp*nbr>40:
            print(f"Warning! There are {nlp*nbr} to plot")
            print(brs)
            print(lps)
        use_files=self.filterFiles(lps=lps,brs=brs,one_repeat=True)

        cell_data=Parallel(n_jobs=12,verbose=100,backend='multiprocessing')(
                    delayed(self.getCellDatPar)(file) for file in use_files
                    )
        all_data=cell_data

        df=pd.DataFrame(all_data,columns=['LinkingProb','BendRig','cells'])
        df.sort_values(by=[ 'LinkingProb','BendRig' ],inplace=True)
        lps=np.sort(df['LinkingProb'].unique())
        brs=np.sort(df['BendRig'].unique())
        nlp=len(lps)
        nbr=len(brs)
        print(lps)
        print(brs)
        add=lambda t,s: s+t
        fname=f"{self.analysis_dir}/sizes_({reduce(add,[f'{v},' for v in lps])})"
        fname+=f"_({reduce(add,[f'{v},' for v in brs])}).txt"
        print(f"{fname=}")
        try:
            widths=np.loadtxt(fname.replace('sizes',"widths"))
            heights=np.loadtxt(fname.replace('sizes',"heights"))
            if len(heights)==0:
                heights=np.array([heights])
        except Exception as e:
            sizes=np.zeros((nlp,nbr,2))
            for jj,br in enumerate(brs):
                for ii,lp in enumerate(lps):
                    print(f"{lp=} {br=}")
                    rr=df.loc[(df['LinkingProb'] == lp) & (df['BendRig'] == br)]
                    if len(rr.index)==0:
                        continue
                    cells=rr['cells'].iloc[0]
                    cell_dat=np.array([ cell.rcm[:2] for cell in cells ])
                    sizes[ii,jj]=np.max(cell_dat,axis=0)-np.min(cell_dat,axis=0)

            print(sizes)
            widths=np.max(sizes[:,:,0],axis=1)
            heights=np.max(sizes[:,:,1],axis=0)
            print(widths)
            print(heights)
            np.savetxt(fname.replace('sizes',"widths"),widths)
            np.savetxt(fname.replace('sizes',"heights"),heights)

        all_cells=[]
        origin=np.array([0.5*widths[0],0.5*heights[0]])
        xticks=[]
        yticks=[]
        for jj,br in enumerate(brs):
            for ii,lp in enumerate(lps):
                print(f"{lp=} {br=}")
                rr=df.loc[(df['LinkingProb'] == lp) & (df['BendRig'] == br)]
                if len(rr.index)==0:
                    continue

                coords=np.array([np.sum(widths[:ii]) +widths[ii]*0.5  + 20*ii ,
                                 np.sum(heights[:jj])+heights[jj]*0.5 + 20*jj  ])
                xticks.append([lp,coords[0]])
                yticks.append([br,coords[1]])

                cells=rr['cells'].iloc[0]
                for cell in cells:
                    cell.rcm[:2]+=coords
                    cell.colour=RodShapedBacterium.colour_fn(cell.theta)
                all_cells.extend(cells)

        cell_dat=np.array([ cell.rcm[:2] for cell in all_cells ])
        right,up=np.max(cell_dat,axis=0)
        left,down=np.min(cell_dat,axis=0)

        fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
        ax.set_facecolor('k')
        addAllCellsToPlot(all_cells,ax,ax_rng=right-left,
                          alpha=0.8,show_id=False)
        ax.axis('scaled')
        ax.set_ylim([down-20,up+20])
        ax.set_xlim([left-20,right+20])
        ax.set_ylabel(r"$\beta$")
        ax.set_xlabel(r"$p_\ell$")

        xticks=np.array(xticks)
        ax.set_xticks(xticks[:,1])
        ax.set_xticklabels([f"{xx:2.1f}" for xx in xticks[:,0]])

        yticks=np.array(yticks)
        ax.set_yticks(yticks[:,1])
        ax.set_yticklabels([f"{yy:2.1f}" for yy in yticks[:,0]])

        print(f"{self.analysis_dir}/collage_plot.pdf")
        ext=f"_({reduce(add,[f'{v}_' for v in lps])[:-1]})"
        ext+=f"_({reduce(add,[f'{v}_' for v in brs])[:-1]})"
        if full:
            ext+="_full"
        os.makedirs(self.analysis_dir,exist_ok=True)
        fig.savefig(f"{self.analysis_dir}/collage_plot_{ext}.png",
                    format="png",bbox_inches='tight',
                    dpi=500,
                    transparent=False)
        fig.savefig(f"{self.analysis_dir}/collage_plot_{ext}.pdf",
                    format="pdf",bbox_inches='tight',
                    transparent=True)

        plt.show()
        quit()

    def clusterColonies(self,variables):
        brs=[]
        lps=[]
        for ii,infile in enumerate(self.in_files):
            run,repeat,stem = self.getRunRepeatStem(infile)
            brs.append(self.hps[run]['BendRig'])
            lps.append(self.hps[run]['LinkingProb'])
        brs=np.sort(np.unique(np.array(brs)))
        lps=np.sort(np.unique(np.array(lps)))

        fig_dir=f"{self.analysis_dir}/clusterColonies"

        os.makedirs(fig_dir,exist_ok=True)

        use_files=self.filterFiles(lps=lps,brs=brs,one_repeat=False)

        df=pd.read_csv(f"{self.analysis_dir}/all_data/general_all_data.txt")

        # Uncomment to see if some did not finish
        # grouped=df.groupby(["run","LinkingProb","BendRig"])['repeat'].sum()
        # ss=grouped[grouped!="repeat3repeat2repeat0repeat4repeat1"]
        grouped=df.groupby(["run","LinkingProb","BendRig"])['repeat'].count()
        # ss=grouped[grouped!=5]
        # print(ss)
        # print("Check the new ones have run")
        # print(df[df["run"]=='run378']['repeat'])
        # names=[
        #        # "roughness",
        #        "aspect_ratio",
        #        "density",
        #        # "S",
        #        # "curvature_L1",
        #        # "weighted_pair_correlation",
        #        "link_angle"
        #        ]
        lbls={
            "roughness": r"$\nu$",
            "aspect_ratio":r"$\alpha$",
            "density":r"$\rho$",
            # "link_angle": r"$\langle \theta^2 \rangle$",
            # "fractal_dimension": r"$\mathcal{D}$",
            "S" : "$S$",
            "curvature_L1" : r"$\mathcal{R}$",
            "weighted_pair_correlation" : r"$\tilde{\ell}_c$",
            # "KMeanLabel": "cluster label"
            }
        print(f"{df.columns=}")
        cols=['LinkingProb','BendRig']
        cols.extend(lbls.keys())
        # cols=['LinkingProb','BendRig',"roughness","aspect_ratio","density"]
        # cols=['LinkingProb','BendRig',"aspect_ratio","density","roughness","cluster_mean"]

        mean_norm=lambda x: (x-x.mean())/x.std()
        maxmin_norm=lambda x: (x-x.min())/(x.max()-x.min())

        # print(df['link_angle'])

        plotConfusion=False
        if plotConfusion:
            for nn,mm in permutations(lbls.items(),2):
                print(f"{nn=} {mm=}")
                fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
                plt.scatter(df[nn[0]],df[mm[0]],fc="None",ec='k')
                plt.xlabel(nn[1])
                plt.ylabel(mm[1])
                fig.savefig(f"{fig_dir}/{nn[0]}_{mm[0]}.pdf",
                            transparent=True,dpi=300,bbox_inches='tight')

        # df['cluster_mean']=(df['cluster_mean']-df['cluster_mean'].mean())/df['cluster_mean'].std()
        # df['cluster_mean']=(df['cluster_mean']-df['cluster_mean'].mean())/(df['cluster_mean'].std())

        for nn in lbls.keys():
            df[nn]=maxmin_norm(df[nn])
        df.sort_values(
            by=['repeat','LinkingProb','BendRig'],
            inplace=True
            )
        # df=df[df['BendRig']>0.01]
        X=df[cols].to_numpy()
        n_clusters=4
        threshold=0.5
        n_comps=3
        models={
            # f'Birch{n_clusters}' : Birch(n_clusters=n_clusters,threshold=1e-2),
            # f'DBSCAN'                      : DBSCAN(eps=0.5,min_samples=5),
            f'KMeans{n_clusters}'          : KMeans(n_clusters=n_clusters,random_state=0),
            f'MiniBatchKMeans{n_clusters}' : MiniBatchKMeans(n_clusters=n_clusters,random_state=10),
            # f'GaussianMixture{n_comps}' : GaussianMixture(n_components=n_comps),
            }

        fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
        for key,model in models.items():
            model.fit(X[:,2:])
            # assign a cluster to each example
            try:
                yhat = model.predict(X[:,2:])
            except Exception as e:
                yhat = model.fit_predict(X[:,2:])

            # print(cols[2:])
            # print(model.cluster_centers_)
            df[key]=yhat

            # # retrieve unique clusters
            clusters = np.unique(yhat)
            single_chain_label=df[
                (df['LinkingProb']==lps.max()) & (df['BendRig']>1)
                ][[key]].mode()[key][0]

            no_buckling_label=df[
                (df['LinkingProb']==lps.min()) & (df['BendRig']<0.5)
                ][[key]].mode()[key][0]

            eye=df[df['aspect_ratio']==df['aspect_ratio'].max()]
            eye_label=eye[key].iloc[0]

            tick_labels={ii:'none' for ii in range(n_clusters)}
            tick_labels[single_chain_label]='single chain'
            tick_labels[no_buckling_label]='no buckling'
            tick_labels[eye_label]='lenticular'
            for ii in range(n_clusters):
                if tick_labels[ii]=='none':
                    tick_labels[ii]='collective'
                    collective_label=ii
                    break
            # print(f"{tick_labels=}")
            lbl_dict={lbl:ii for ii,lbl in
                      enumerate([no_buckling_label,eye_label,collective_label,
                                 single_chain_label])
                     }
            def lbl_map(jj):
                return lbl_dict[jj]
            df[key]=df[key].apply(lbl_map)
            use_tick_labels={}
            for ii,jj in lbl_dict.items():
                use_tick_labels[jj]=tick_labels[ii]

            # create scatter plot for samples from each cluster
            for ii,cluster in enumerate(clusters):
                # get row indexes for samples with this cluster
                row_ix = np.where(yhat == cluster)
                # create scatter of these samples
                ax.scatter(X[row_ix,0], X[row_ix,1],label=f'cluster{ii}')


            aidan_chain_fit=lambda p: (2/(1-p))**(4) / 4e6

            addSingleChainToPhaseDiagram(ax,lps)
            ax.plot(lps,aidan_chain_fit(lps),
                    'k:',label=r'$p_\ell=1-\frac{1}{N_c}$',
                    lw=2
                    )
            ax.set_yscale('log')
            ax.set_xlabel(r'$p_\ell$')
            ax.set_ylabel(r'$\beta$')
            ax.set_ylim([0.85*brs.min(),1.15*brs.max()])
            plt.legend(frameon=True,loc='lower left')
            plt.savefig(f"{fig_dir}/cluster_plot_{key}.png")
            plt.close()
            # plt.show()

            print(df[['LinkingProb','BendRig',key]])

            # Try to find the mode again
            # cluster_filter_df=df.groupby(["LinkingProb","BendRig"])
            # cluster_filter_df=cluster_filter_df[key].apply(pd.Series.mode)
            # cluster_filter_df=cluster_filter_df.to_frame().reset_index()
            # print(cluster_filter_df[-50:])
            # print(cluster_filter_df.dtypes)
            # quit()

            df[key]=df[key].apply(lambda ii: [ii])
            cluster_filter_df=df[['LinkingProb','BendRig',key]].groupby(["LinkingProb","BendRig"])
            print(cluster_filter_df.get_group((0.90,13.13680)))
            cluster_filter_df=cluster_filter_df[key].sum()
            cluster_filter_df=cluster_filter_df.to_frame().reset_index()
            cluster_filter_df[key]=cluster_filter_df[key].apply(
                lambda ii: np.mean(scipy.stats.mode(ii)[0])
                )
            def smth(row):
                res=cluster_filter_df[
                    (cluster_filter_df['LinkingProb']==row['LinkingProb']) &
                    (cluster_filter_df['BendRig']==row['BendRig'])
                    ][key].iloc[0]
                return res
            df[key]=df.apply(smth,axis=1)
            cluster_filter_df.sort_values(
                by=['LinkingProb','BendRig'],
                inplace=True
                )
            tot_sz=0
            for i in range(5):
                sz=len(df[df['repeat']==f'repeat{i}'])
                print(f"{sz=}")
                tot_sz+=sz
            print(f"{tot_sz=}")

            # lps=np.sort(df['LinkingProb'].unique())
            # brs=np.sort(df['BendRig'].unique())
            # reps=np.sort(df['repeat'].unique())
            nlp=len(lps)
            nbr=len(brs)
            # nrp=len(reps)
            # x=df['LinkingProb']
            # y=df['BendRig']
            # c=df[key].to_numpy().reshape((nlp,nbr,nrp))
            # xx=x.to_numpy().reshape((nlp,nbr,nrp))
            # yy=y.to_numpy().reshape((nlp,nbr,nrp))
            xx=cluster_filter_df['LinkingProb'].to_numpy().reshape((nlp,nbr))
            yy=cluster_filter_df['BendRig'].to_numpy().reshape((nlp,nbr))
            c=cluster_filter_df[key].to_numpy().reshape((nlp,nbr))
            print(f"{c.shape=}")
            fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
            # cmap=cm.get_cmap('viridis',n_clusters)
            CS=ax.pcolor(xx,yy,c,
                         # norm=mpl.colors.LogNorm(vmin=c.min(),vmax=c.max()),
                         shading='auto',
                         # cmap=cmap
                         )
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar=fig.colorbar(CS,cax=cax,ticks=np.arange(n_clusters))
            cbar.set_label(rf"buckle type",size=14)
            cbar.ax.tick_params(labelsize=12)
            cbar.ax.set_yticklabels([use_tick_labels[ii] for ii in range(n_clusters)])
            addSingleChainToPhaseDiagram(ax,lps)
            addCollectiveToPhaseDiagram(ax,lps)

            ixs=np.where(c>0)
            bb=brs[np.min(ixs[1])]
            mb=lps[np.min(ixs[0])]
            addPhaseBoundaries(ax,mb=mb,bb=bb)
            ax.set_xlabel(rf"$p_\ell$")
            ax.set_ylabel(rf"$\beta$")
            ax.set_yscale('log')
            ax.legend(labelcolor='w',loc='lower left',frameon=True,
                      facecolor=cm.get_cmap('viridis')(0),
                      framealpha=0.9
                      )
            ax.set_ylim(brs.min(),brs.max())
            fig.savefig(f"{fig_dir}/pcolor_buckles_{key}.pdf",
                        format="pdf",
                        bbox_inches='tight',
                        transparent=True)
            plt.show(block=True)
            plt.close('all')
        self.makeCollagePlotLabelled(brs=brs,
                                     lps=lps,
                                     fig_dir=fig_dir,
                                     df=df,
                                     label_cols=models.keys(),
                                     type_labels=use_tick_labels)

    @staticmethod
    def evalRep(eval_dat,br,lp):
        """
        Try to guess the colony classification based on some criterion
        """
        eye=((eval_dat['aspect_ratio'].iat[0]>1.5) &
             (eval_dat['density'].iat[0]>1) &
             (eval_dat['S'].iat[0]>0.995)
             )

        col=((eval_dat['aspect_ratio'].iat[0]>1.5) &
             (eval_dat['density'].iat[0]>1) &
             (eval_dat['S'].iat[0]<=0.995)
             )

        notbuckled=((eval_dat['aspect_ratio'].iat[0]<=1.1) &
                    (eval_dat['density'].iat[0]>=1)
                    )
        assignment=None
        if br<=0.1 or lp<=0.4:
            assignment='n'
        elif br<0.3 and lp>=0.95:
            assignment='a'
        elif br>=0.3 and lp>=0.95:
            assignment='s'
        elif br>=0.5 and lp==0.9:
            assignment='s'
        elif br<0.5 and lp==0.9:
            assignment='a'
        elif br>13 and lp==0.85:
            assignment='e'
        elif br>=0.7 and br<=13 and lp==0.85:
            assignment='c'
        elif br<0.7 and lp==0.85:
            assignment='a'
        elif eye:
            assignment='e'
            print("eye")
        elif col:
            assignment='c'
            print("col")
        elif notbuckled:
            assignment='n'
            print("not buckled")
        return assignment

    def assignBuckles(self,variables,full=False):
        """
        This function asks for manual assignment of the colony morphologies,
        and then makes a phase diagram
        """
        brs=[]
        lps=[]
        for ii,infile in enumerate(self.in_files):
            run,repeat,stem = self.getRunRepeatStem(infile)
            brs.append(self.hps[run]['BendRig'])
            lps.append(self.hps[run]['LinkingProb'])
        brs=np.sort(np.unique(np.array(brs)))
        lps=np.sort(np.unique(np.array(lps)))

        if full==False:
            lps=lps[::10]
            # brs=brs[brs>1e-2]
            brs=brs[::10]
        else:
            print("Warning - running all LinkingProb and BendRig")
        # add=lambda t,s: s+t
        # list_brs=f"({reduce(add,[f'{v},' for v in brs])[:-1]})"
        # list_lps=f"({reduce(add,[f'{v},' for v in lps])[:-1]})"
        # fig_dir=f"{self.analysis_dir}/buckle_figs_lp_{list_lps}_bd_{list_brs}"
        fig_dir=f"{self.analysis_dir}/buckle_fig_all"

        os.makedirs(fig_dir,exist_ok=True)

        use_files=self.filterFiles(lps=lps,brs=brs,one_repeat=False)
        print(f"Running for {len(use_files)} files")

        # Map of buckle type to label
        equiv = {'n': "no buckling",
                 # 'a': "ambiguous",
                 'e': "lenticular",
                 # 'w': "weak collective",
                 'c': "collective",
                 's': "single chain"}

        # Filename of the output phase plot
        # fname=f"{fig_dir}/buckled_phase_plot_lp_{list_lps}_bd_{list_brs}.txt"
        other_data_df=pd.read_csv(f"{self.analysis_dir}/all_data/general_all_data.txt")
        # names=["roughness","aspect_ratio","density","S",
        #        "curvature_L1",
        #        "weighted_pair_correlation"
        #        ]
        lbls={
            "roughness": r"$\nu$",
            "aspect_ratio":r"$\alpha$",
            "density":r"$\rho$",
            "link_angle": r"$\langle \theta^2 \rangle$",
            # "fractal_dimension": r"$\mathcal{D}$",
            "S" : "$S$",
            "curvature_L1" : r"$\mathcal{R}$",
            "weighted_pair_correlation" : r"$\tilde{\ell}_c$",
            "persistence_length" : r"$\ell_p$"
            # "KMeanLabel": "cluster label"
            }
        fname=f"{fig_dir}/buckled_phase_all.txt"
        try:
            all_data=pd.read_csv(fname)
        except Exception as e:
            def makeAllFigs(file):
                run,repeat,stem = self.getRunRepeatStem(file)
                lp=self.hps[run]['LinkingProb']
                br=self.hps[run]['BendRig']
                stm=f"buckle_figs_lp_{lp}_br_{br}_rep_{repeat}"
                print(file,stm)
                fig_name=f'{fig_dir}/{stm}.png'
                if os.path.exists(fig_name):
                    # print(f"{fig_name} already exists")
                    return

                print(f"Reading {lp=} {br=}")
                data = pd.read_csv(file,sep="\t")
                cells,_ = intialiseElementsFromData(data,self.class_dict)
                for cell in cells:
                    cell.colour=RodShapedBacterium.colour_fn(cell.theta)

                cell_dat=np.array([ cell.rcm[:2] for cell in cells ])
                right,up=np.max(cell_dat,axis=0)
                left,down=np.min(cell_dat,axis=0)

                fig,ax=plt.subplots(1,1,)
                addAllCellsToPlot(cells,ax,ax_rng=right-left,
                                  alpha=0.8,show_id=False)
                eval_dat=other_data_df[
                        (other_data_df['LinkingProb']==lp) &
                        (other_data_df['BendRig']==br) &
                        (other_data_df['repeat']==repeat)
                        ]
                for ii,nn in enumerate(lbls.keys()):
                    try:
                        ax.text(0.1,0.2+ii*0.1,lbls[nn]+f"={eval_dat[nn].iat[0]:3.3f}",
                                horizontalalignment='center',
                                verticalalignment='center',
                                fontsize="x-large",
                                transform=ax.transAxes)
                    except Exception as e:
                        ax.text(0.1,0.2+ii*0.1,lbls[nn]+f"={eval_dat[nn].iat[0]}",
                                horizontalalignment='center',
                                verticalalignment='center',
                                fontsize="x-large",
                                transform=ax.transAxes)

                plt.axis('off')
                plt.show(block=False)
                ax.axis('scaled')
                ax.set_ylim([down-20,up+20])
                ax.set_xlim([left-100,right+20])
                ax.axis("off")
                fig.savefig(fig_name,bbox_inches='tight',dpi=400)
                plt.close()

            Parallel(n_jobs=-1,verbose=100)(
                        delayed(makeAllFigs)(file) for file in use_files
                        )

            all_data=other_data_df.copy(deep=True)
            all_data['BuckleType']='n'
            all_data.sort_values(by=['LinkingProb','BendRig'],inplace=True)

            """
            Find the number of runs and repeats
            Go through runs sequentially and display all repeats at once
            All repeats will be assigned the predominant buckle type
            """
            def getRunInt(file):
                run,repeat=self.getRunRepeatStem(file)[:2]
                run=run.strip('run')
                return int(run),repeat

            use_files=sorted(use_files,
                             key=getRunInt
                             )
            continue_flag=False
            for rr in reversed(range(1+getRunInt(use_files[-1])[0])):
                lp=self.hps[f'run{rr}']['LinkingProb']
                br=self.hps[f'run{rr}']['BendRig']
                if continue_flag:
                    if br==skip_br:
                        print(f"Skipping {lp=} {br=}")
                        continue
                    else:
                        continue_flag=False
                load_fname=f"{self.analysis_dir}/run{rr}/"
                reps=glob(f"{load_fname}/*")
                nreps=len(reps)
                fig,ax=plt.subplots(1,nreps,figsize=[nreps*10,10])
                for ii,rep in enumerate(reps):
                    print(f"Loading {rep}")
                    run,repeat,stem = self.getRunRepeatStem(rep)
                    lp=self.hps[run]['LinkingProb']
                    br=self.hps[run]['BendRig']

                    stm=f"buckle_figs_lp_{lp}_br_{br}_rep_{repeat}"
                    fig_name=f'{fig_dir}/{stm}.png'

                    img=plt.imread(fig_name)
                    imgplot=ax[ii].imshow(img)
                    ax[ii].axis('off')

                plt.show(block=False)
                print(equiv)
                assignment=input("Assign buckled: ")
                plt.close('all')

                all_data.loc[all_data['run']==f'run{rr}',
                             ['BuckleType']
                             ]=assignment
                if assignment=='n':
                    continue_flag=True
                    skip_br=br
            all_data.sort_values(by=['LinkingProb','BendRig','repeat'],inplace=True)
            all_data.to_csv(fname,index=False)
            # all_data["labelBuckle"] = all_data["BuckleType"].map(equiv)

            # quit()
            # for count,file in enumerate(use_files):
            #     if count%int(0.1*len(use_files))==0:
            #         print(f"=====\nprogress={100*count/len(use_files):3.2f} %\n")
            #
            #     run,repeat,stem = self.getRunRepeatStem(file)
            #     current_data=[]
            #     lp=self.hps[run]['LinkingProb']
            #     br=self.hps[run]['BendRig']
            #
            #     print(f"{lp=} {br=} {repeat=}")
            #     eval_dat=other_data_df[
            #             (other_data_df['LinkingProb']==lp) &
            #             (other_data_df['BendRig']==br) &
            #             (other_data_df['repeat']==repeat)
            #             ]
            #     if len(eval_dat)!=1:
            #         print(eval_dat)
            #         quit()
            #
            #     assignment=evalRep(eval_dat,br,lp)
            #     if assignment==None:
            #         stm=f"buckle_figs_lp_{lp}_br_{br}_rep_{repeat}"
            #         fig_name=f'{fig_dir}/{stm}.png'
            #
            #         fig,ax=plt.subplots(1,1,figsize=[10,10])
            #         img=plt.imread(fig_name)
            #         imgplot=ax.imshow(img)
            #         plt.axis('off')
            #         plt.show(block=False)
            #
            #         print(equiv)
            #         assignment=input("Assign buckled: ")
            #
            #         plt.close('all')
            #
            #     current_data.extend([ lp,br,repeat,assignment ])
            #     current_data.extend([ eval_dat[nn].iat[0] for nn in names ])
            #     all_data.append(current_data)

            # pprint(all_data)
            # columns=['LinkingProb','BendRig','repeat','BuckleType']
            # columns.extend(names)
            # print(columns)
            # df=pd.DataFrame(all_data,
            #     columns=columns
            #     )
            # df["labelBuckle"] = df["BuckleType"].map(equiv)
            # df.sort_values(by=['LinkingProb','BendRig','repeat'],inplace=True)
            # df.to_csv(fname,index=False)

        # all_data["labelBuckle"] = all_data["BuckleType"].map(equiv)
        # pprint(all_data)
        # quit()

        ndf=all_data[['LinkingProb','BendRig',"BuckleType"]].copy()
        numericequiv = { kk[1]:kk[0] for kk in enumerate(equiv.keys()) }
        print(numericequiv)
        def mapBucklesToNum(buckle):
            unq_buckles=set(buckle)
            avg=sum([numericequiv[bb] for bb in unq_buckles])/len(unq_buckles)
            # type=[int(np.floor(avg))]
            # if len(unq_buckles)>2:
            #     avg=-1
            #     print(f"{avg=} {unq_buckles=} {type=}")
            # elif len(unq_buckles)>1:
            #     type=[int(np.floor(avg)),int(np.ceil(avg))]
            #     print(f"{avg=} {unq_buckles=} {type=}")
            return avg
        ndf['numericBuckleType']=ndf['BuckleType'].apply(mapBucklesToNum)
        u_vals=ndf['numericBuckleType'].unique()
        u_vals.sort()
        new_map={ kk[1]:kk[0] for kk in enumerate(u_vals) }
        ndf['numericBuckleType']=ndf['numericBuckleType'].map(new_map)
        # all_data['numericBuckleType']=ndf['numericBuckleType']
        lblequiv = { kk[0]:kk[1] for kk in enumerate(equiv.keys()) }
        # for v in u_vals:
        #     if v not in lblequiv:
        #         lv=int(np.floor(v))
        #         lblequiv[v]=f"{lblequiv[lv+1]}/{lblequiv[lv]}"
        ylabels=[]
        for v in sorted(u_vals):
            if v not in lblequiv:
                lv=int(np.floor(v))
                lbl=f"{equiv[lblequiv[lv+1]]}/\n{equiv[lblequiv[lv]]}"
            else:
                lbl=equiv[lblequiv[v]]
            ylabels.append(lbl)

        # ndf["numericBuckleType"] = ndf["BuckleType"].map(numericequiv)
        # print(ndf)
        grouped=ndf.groupby(['LinkingProb','BendRig'])
        ndf=grouped["numericBuckleType"].max().reset_index(level=[0,1])
        # ndf["BuckleType"]=ndf["numericBuckleType"].map(lblequiv)
        # print(ndf)
        print(ndf)

        lps=np.sort(ndf['LinkingProb'].unique())
        brs=np.sort(ndf['BendRig'].unique())
        nlp=len(lps)
        nbr=len(brs)
        # print(lps)
        # print(brs)

        x=ndf['LinkingProb']
        y=ndf['BendRig']
        c=ndf["numericBuckleType"]
        xx=x.to_numpy().reshape((nlp,nbr))
        yy=y.to_numpy().reshape((nlp,nbr))
        cc=c.to_numpy().reshape((nlp,nbr))
        ixs=np.where(cc>0)
        bb=brs[np.min(ixs[1])]
        mb=lps[np.min(ixs[0])]
        fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
        CS=ax.pcolor(xx,yy,cc,
                     # norm=mpl.colors.LogNorm(vmin=c.min(),vmax=c.max()),
                     shading='auto'
                     )
        addSingleChainToPhaseDiagram(ax,lps)
        addCollectiveToPhaseDiagram(ax,lps)
        addPhaseBoundaries(ax,mb=mb,bb=bb)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar=fig.colorbar(CS,cax=cax,ticks=np.arange(len(ylabels)))
        cbar.set_label(rf"buckle type",size=14)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_yticklabels(ylabels)
        ax.set_xlabel(rf"$p_\ell$")
        ax.set_ylabel(rf"$\beta$")
        ax.set_yscale('log')
        ax.set_ylim(brs.min(),brs.max())
        ax.legend(labelcolor='w',loc='lower left',frameon=True,
                  facecolor=cm.get_cmap('viridis')(0),
                  framealpha=0.9
                  )
        fig.savefig(f"{fig_dir}/pcolor_buckles.pdf",
                    format="pdf",
                    bbox_inches='tight',
                    transparent=True)
        plt.show(block=True)
        # quit()

        use_tick_labels={ii:lbl for ii,lbl in enumerate(ylabels)}
        self.makeCollagePlotLabelled(brs=brs,
                                     lps=lps,
                                     fig_dir=fig_dir,
                                     df=ndf,
                                     label_cols=['numericBuckleType'],
                                     type_labels=use_tick_labels)

        # self.makeCollagePlotLabelled(brs=brs,
        #                              lps=lps,
        #                              fig_dir=fig_dir,
        #                              df=df,
        #                              label_cols=["BuckleType"])

    def getCellDatPar(self,file):
        run,repeat,stem = self.getRunRepeatStem(file)
        data = pd.read_csv(file,sep="\t")
        cells,_ = intialiseElementsFromData(data,self.class_dict)
        return [self.hps[run]['LinkingProb'],
                self.hps[run]['BendRig'],
                cells]

    @staticmethod
    def calcSizes(fig_dir,df,lps,brs):
        fname=f"{fig_dir}/sizes_lps_{lps[0]}_{lps[-1]}_brs_{brs[0]}_{brs[-1]}.txt"
        nlp=len(lps)
        nbr=len(brs)
        try:
            widths=np.loadtxt(fname.replace('sizes',"widths"))
            heights=np.loadtxt(fname.replace('sizes',"heights"))
            if len(heights)==0:
                heights=np.array([heights])
        except Exception as e:
            sizes=np.zeros((nlp,nbr,2))
            for jj,br in enumerate(brs):
                for ii,lp in enumerate(lps):
                    print(f"{lp=} {br=}")
                    rr=df.loc[(df['LinkingProb'] == lp) & (df['BendRig'] == br)]
                    if len(rr.index)==0:
                        continue
                    cells=rr['cells'].iloc[0]
                    cell_dat=np.array([ cell.rcm[:2] for cell in cells ])
                    sizes[ii,jj]=np.max(cell_dat,axis=0)-np.min(cell_dat,axis=0)

            print(sizes)
            widths=np.max(sizes[:,:,0],axis=1)
            heights=np.max(sizes[:,:,1],axis=0)
            print(widths)
            print(heights)
            np.savetxt(fname.replace('sizes',"widths"),widths)
            np.savetxt(fname.replace('sizes',"heights"),heights)
        return widths,heights

    def makeCollagePlotLabelled(self,brs,lps,fig_dir,df,label_cols,type_labels):
        """
            Make a collage plot and label each picture
        """
        brs=brs[brs>0.5]
        nlp=len(lps)
        nbr=len(brs)
        if nlp*nbr>40:
            brs=brs[brs>0.5]
            brs=brs[::2]
            lps=lps[lps>=0.2]
            lps=lps[::2]
            nlp=len(lps)
            nbr=len(brs)
            # while nlp*nbr>20:
            #     lps=lps[::2]
            #     nlp=len(lps)
            #     nbr=len(brs)
        if nlp*nbr>40:
            print(f"Warning! There are {nlp*nbr} to plot")
            print(brs)
            print(lps)
        use_files=self.filterFiles(lps=lps,brs=brs,one_repeat=True)

        cell_data=Parallel(n_jobs=12,verbose=100,backend='multiprocessing')(
                    delayed(self.getCellDatPar)(file) for file in use_files
                    )

        for lbl_indx,label_col in enumerate(label_cols):

            dfcell=pd.DataFrame(copy(cell_data),columns=['LinkingProb','BendRig','cells'])
            dfplot=pd.merge(df,dfcell,on=['LinkingProb','BendRig'])
            widths,heights=self.calcSizes(fig_dir,dfplot,lps,brs)

            nlabels=df[label_col].nunique()
            fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=2))
            ax.set_facecolor('k')
            norm = plt.Normalize(0, nlabels-1)
            handles=[mpl.patches.Patch(color=cm.get_cmap('viridis')(
                     norm(lbl)
                     ),
                     label=f'{type_labels[lbl]}') for lbl in range(nlabels)
                     ]
            all_cells=[]
            origin=np.array([0.5*widths[0],0.5*heights[0]])
            xticks=[]
            yticks=[]
            contours=[]
            for jj,br in enumerate(brs):
                for ii,lp in enumerate(lps):
                    print(f"{lp=} {br=}")
                    rr=dfplot.loc[(dfplot['LinkingProb'] == lp) & (dfplot['BendRig'] == br)]
                    if len(rr.index)==0:
                        continue

                    coords=np.array([np.sum(widths[:ii]) +widths[ii]*0.5  + 20*ii ,
                                     np.sum(heights[:jj])+heights[jj]*0.5 + 20*jj  ])
                    xticks.append([lp,coords[0]])
                    yticks.append([br,coords[1]])
                    label=rr[label_col].iloc[0]
                    # ax.text(*coords,f"type={label}")

                    cells=rr['cells'].iloc[0]
                    for cell in cells:
                        if lbl_indx==0:
                            cell.rcm[:2]+=coords
                        cell.colour=cm.get_cmap('viridis')(norm(label))
                        # cell.colour=RodShapedBacterium.colour_fn(cell.theta)
                        # cm.get_cmap('hsv')
                    all_cells.extend(cells)
                    tmp_contour=computeColonyContour(
                        cells,add_links=True,ret_gon=False
                        )
                    contours.append(tmp_contour)

            cell_dat=np.array([ cell.rcm[:2] for cell in all_cells ])
            # np.savetxt(
            #     f"{fig_dir}/collage_plot_{label_col}_all_cells_data.txt",
            #     cell_dat
            #     )
            # for jj,contour in enumerate(contours):
            #     np.savetxt(
            #         f"{fig_dir}/collage_plot_{label_col}_contour_{jj}.txt",
            #         cell_dat
            #         )
            right,up=np.max(cell_dat,axis=0)
            left,down=np.min(cell_dat,axis=0)

            ax_rng=right-left
            addAllCellsToPlot(all_cells,ax,ax_rng=ax_rng,
                              alpha=0.8,show_id=False,ec='w')
            ax.add_collection(
                LineCollection(contours,lw=80/ax_rng,colors='w')
                )
            ax.axis('scaled')
            ax.set_ylim([down-20,up+20])
            ax.set_xlim([left-20,right+20])
            ax.set_ylabel(r"$\beta$")
            ax.set_xlabel(r"$p_\ell$")

            xticks=np.array(xticks)
            ax.set_xticks(xticks[:,1])
            ax.set_xticklabels([f"{xx}" for xx in xticks[:,0]])

            yticks=np.array(yticks)
            ax.set_yticks(yticks[:,1])
            ax.set_yticklabels([f"{yy:2.1f}" for yy in yticks[:,0]])

            ax.legend(
                handles=handles,
                loc='lower left',
                frameon=True,
                framealpha=0.7
                )
            os.makedirs(fig_dir,exist_ok=True)
            fig.savefig(f"{fig_dir}/collage_plot_{label_col}.png",
                        format="png",bbox_inches='tight',
                        dpi=400)
            # fig.savefig(f"{self.analysis_dir}/collage_plot_{label_col}.pdf",
            #             format="pdf",bbox_inches='tight',
            #             transparent=True)
            plt.show(block=False)
        plt.show()
        plt.close('all')

    def findCollectiveBucklingPoint(self,df,variables,keys):
        nx=df[variables[0]['iv']].nunique()
        ny=df[variables[1]['iv']].nunique()
        df.sort_values(by=[ v['iv'] for v in (variables) ],inplace=True)
        avg=df.groupby('run').mean()
        avg.sort_values(by=[ v['iv'] for v in (variables) ],inplace=True)
        avg=avg.reset_index()
        avg=avg[[variables[0]['iv'],variables[1]['iv'],*keys]]
        brs=df[variables[1]['iv']].unique()
        for br in brs:
            tmp_df=avg[avg['BendRig']==br]
            print(tmp_df)
            plt.title(f"{br=}")
            for key in keys:
                plt.plot(tmp_df['LinkingProb'],np.abs(tmp_df[key]),label=key)
            plt.plot(tmp_df['LinkingProb'],(1/(1-tmp_df['LinkingProb'])),'--',
                     label="av chain length")
            plt.legend()
            plt.show()
        quit()

    @staticmethod
    def plotEvolutionQuantity(quantity,analysis_dir,qname,label=None,poster=False):
        """
            Plot the time evolution of a quantity for several colony sizes

            Parameters:
                quantity: dict
                    fields: lps are the linking probs, qname is the mean
                            and err_qname is the sem
                    analysis_dir: str
                        path to where the results will be collated
                    qname: str
                        The quantity name to be investigated
                    label: str (optional)
                        label for legend/axis
                    poster: bool (optional)
                        set fig size and font size for poster
        """
        if not label:
            label=qname.capitalize()

        if poster:
            figsize=[10,10/1.61803398875]
            fontsize=30
            pp_params={
                'axes.titlesize'        : fontsize,
                'axes.labelsize'        : fontsize,
                'xtick.labelsize'       : fontsize,
                'ytick.labelsize'       : fontsize,
                'legend.fontsize'       : fontsize,
                'font.size'             : fontsize,
            }
            mpl.rcParams.update(pp_params)
            lw=2.5
            markersize=12

        fig,ax = plt.subplots(1,1,figsize=figsize)
        for key,qs_at_ss in quantity.items():
            print(key,qs_at_ss)
            lps=qs_at_ss['lps']
            qs=qs_at_ss[f'{qname}']
            errs=qs_at_ss[f'err_{qname}']
            try:
                if np.any(np.isnan(errs)):
                    ax.plot(lps,qs,
                            marker="o",mfc="w",mec='k',
                            label=f"{float(key.strip('s_')):3.0f}",
                            lw=lw,
                            markersize=markersize,
                            alpha=0.7)
                else:
                    ax.errorbar(
                            lps,qs,
                            yerr=errs,
                            marker="o",mfc="w",mec='k',
                            label=f"@ {float(key.strip('s_')):3.0f}",
                            capsize=5,alpha=0.7)
            except Exception as e:
                print(e)
                ax.plot(lps,qs,
                        marker="o",mfc="w",mec='k',
                        label=f"@ {key:3.0f}")
        ax.set_ylabel(f"Colony {label}")
        if poster:
            ax.set_xlabel("$p_\ell$")
        else:
            ax.set_xlabel("Linking Probability")
        plt.legend(title="colony size")

        if poster:
            fig.savefig(f"{analysis_dir}/poster_{qname}.pdf",
                        format="pdf",bbox_inches='tight',
                        transparent=True)
        else:
            fig.savefig(f"{analysis_dir}/{qname}.pdf",
                        format="pdf",bbox_inches='tight',
                        transparent=True)
        plt.show()

    @staticmethod
    def getColonySize(file,cells):
        return reduce(lambda a,b: a+b,
                      [cell.length+cell.diameter for cell in cells]
                      )

    @staticmethod
    def plotBoxCount(fractal_dim,boxes,coords,fractal_filename,cells):
        contourgon = Polygon([(v[0],v[1]) for v in coords])

        fig,ax = plt.subplots(1,1,figsize=setFigSize(8*247))
        ax.plot(*contourgon.exterior.xy,'m',lw=0.7,zorder=3)
        for box in boxes:
            ax.plot(*box.exterior.xy,'k',alpha=0.5,lw=0.7,zorder=1)

        #----
        # springs = [ spring for cell in cells
        #             for spring in cell.getSpringGons(radiusfactor=0.7)
        #             ]
        # ax.add_collection(
        #     LineCollection(springs,colors='r',zorder=4)
        #     )
        # cellgons=[ cell.getPolygon(0.7*cell.radius) for cell in cells ]
        # patches = [ PolygonPatch(polygon=gon,alpha=0.4) for gon in cellgons ]
        # ax.add_collection(PatchCollection(patches,match_original=True,zorder=1))
        #----
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.axis("scaled")
        # plt.show()
        cntr=contourgon.centroid
        ax.text(cntr.x,cntr.y,rf"$\mathcal{{D}}=${fractal_dim:3.3f}",
                horizontalalignment='center',verticalalignment='center'
                )
        fig.savefig(fractal_filename,format="pdf",bbox_inches='tight',
                    transparent=True)
        plt.close()

    @staticmethod
    def plotCurvature(L1,L2,rs,skel_rs,coords,cells,fname):

        # find the centre of the colony and the rough radius
        params=fitEllipse(coords)
        xc, yc, a, b, theta = params
        circ=Point([xc,yc]).buffer(1)
        ellipse=rotate(scale(circ,0.5*a,0.5*b),theta,use_radians=True)

        fig,ax=plt.subplots()
        RodShapedBacterium.sus_vis_radius_factor=0.7
        for cell in cells:
            cell.addElementToPlot(ax,zorder=1)
        for rr in rs:
            p=rr[:2]
            r=rr[2]
            ax.plot(*p,'o',mec='k',mfc='w',alpha=0.8,zorder=3,
                    markersize=20/(max(a,b)),
                    mew=5/(max(a,b)))
            ax.plot(*Point(p).buffer(r).exterior.xy,c=cm.magma(r),alpha=0.8,
                    zorder=3,lw=5/(max(a,b)))
            # print(f"{p} is {r} from boundary")
        ax.axis('scaled')
        ax.plot(coords[:,0],coords[:,1],label="contour",alpha=0.8)
        ax.set_xlim(np.min(coords[:,0]),np.max(coords[:,0]))
        ax.set_ylim(np.min(coords[:,1]),np.max(coords[:,1]))
        fig.savefig(fname,format='pdf',bbox_inches='tight',transparent=True)
        # plt.show()
        plt.close()

        dr=0.5
        density_field, density_grid = RodShapedBacterium.computeDensity(
            cells,
            dr=dr,
            fname=fname.replace('rs','rho')
            )
        density_field[density_field>0.5] =1
        density_field[density_field<=0.5]=0
        # skeleton=skeletonize(density_field)
        # Compute the medial axis (skeleton) and the distance transform
        skel, distance = medial_axis(density_field,return_distance=True)
        dist_on_skel = distance * skel * dr
        density_field[density_field<=0.05]=np.nan
        fig,ax=plt.subplots()
        for cell in cells:
            cell.addElementToPlot(ax,zorder=1)
        # cs=ax.contourf(*density_grid,density_field,
        #                cmap=cm.PuBu_r,
        #                zorder=0)
        dist_on_skel[dist_on_skel<1e-5]=np.nan
        CS=ax.pcolor(*density_grid,dist_on_skel,cmap=cm.magma,
                     zorder=4,
                     shading='auto')
        # prs=rs[ np.logical_and(rs[:,2]>0.95*curvature,rs[:,2]<1.05*curvature) ]
        # ax.plot(prs[:,0],prs[:,1],'o',mec='k',mfc='None',alpha=0.1,zorder=2)

        dists=[]
        for rr in skel_rs:
            if Point(rr[:2]).within(ellipse) and rr[2]>6:
                dists.append(rr)

        for rr in dists:
            p=rr[:2]
            r=rr[2]
            ax.plot(*p,'o',mec='k',mfc='w',alpha=0.8,zorder=3,
                    markersize=12/(max(a,b)),
                    mew=5/(max(a,b)))
            ax.plot(*Point(p).buffer(r).exterior.xy,c=cm.magma(r),alpha=0.8,
                    zorder=3,lw=5/(max(a,b)))
        ax.axis('scaled')
        ax.plot(coords[:,0],coords[:,1],label="contour",alpha=0.8)
        ax.set_xlim(np.min(coords[:,0]),np.max(coords[:,0]))
        ax.set_ylim(np.min(coords[:,1]),np.max(coords[:,1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar=fig.colorbar(CS,cax=cax)
        cbar.set_label("$r$ (\si{\micro\meter})",size=14)
        cbar.ax.tick_params(labelsize=12)
        fig.savefig(fname.replace('.pdf','_dist_on_skel.pdf'),format='pdf',
                    bbox_inches='tight',transparent=True)
        plt.close()
        # plt.show()
        plt.close()

        fig,ax=plt.subplots()
        data=ax.hist(rs[:,2],bins=20,density=True,alpha=0.5,label='random maxes')
        # ax.hist(skel_rs[:,2],bins=20,density=True,
        #         alpha=0.5,label='medial dist')
        ax.set_xlabel('$\max(r)$')
        ax.set_ylabel('$\mathcal{P}(r)$')
        for cc in [[L1,'random maxes','m'],[L2,'medial dist','k']]:
            curvature,lbl,color=cc
            ax.axvline(curvature,linestyle='--',color=color,zorder=2,label=lbl)
        plt.legend()
        fig.savefig(fname.replace('.pdf','hist.pdf'),format='pdf',
                    bbox_inches='tight',transparent=True)
        plt.close()
        # plt.show()

    def getMaxOverlap(self,in_file,cells):
        """
            Run overlap analysis and save histogram
        """
        dir=self.makeAnalysisDirectory(in_file,'overlaps')
        overlap_file=f"{dir}/overlaps.txt"
        if not os.path.exists(overlap_file):
            overlaps=RodShapedBacterium.getOverlaps(cells)
            np.savetxt(overlap_file,overlaps)
            if self.make_output_figures:
                fig,ax=plt.subplots(1,1,figsize=setFigSize(247))
                ax.hist(overlaps,density=True)
                ax.set_xlabel("$h$")
                ax.set_ylabel("$\mathcal{P}(h)$")
                fig.savefig(f"{dir}/overlap_hist.pdf",format="pdf",
                            bbox_inches='tight',transparent=True
                            )
        else:
            overlaps=np.loadtxt(overlap_file)
        return np.max(overlaps)

    def plotLengthScale(self,mean_cluster_area,characteristic_radius,mod='',
                        label=None):
        print(f"{mean_cluster_area[0]=}")
        # or_L = np.sqrt(mean_cluster_area[0,1])
        # err_or_L = 0.5*mean_cluster_area[0,2]/or_L

        or_L = mean_cluster_area[0,1]
        err_or_L = mean_cluster_area[0,2]

        phi=or_L/characteristic_radius[:,1]
        errors=np.abs(phi)*np.sqrt(
            (characteristic_radius[:,2]/characteristic_radius[:,1])**2 +
            (err_or_L/or_L)**2
            )

        _dir=f"{self.analysis_dir}/LengthScales"
        os.makedirs(_dir,exist_ok=True)
        np.savetxt(f"{_dir}/competition_length_scales{mod}.txt",
                   [characteristic_radius[:,0],phi,errors]
                   )
        if self.generate_data_only==False:
            fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
            ax.errorbar(characteristic_radius[:,0],
                        phi,
                        errors,
                        marker='o',
                        mfc='w',
                        alpha=0.7,
                        capsize=5,
                        mec='k'
                        )
            if label==None:
                label=r"$\mathcal{K}\sqrt{\langle A(0) \rangle}$"
            ax.set_ylabel(label)
            ax.set_xlabel(r"$p_\ell$")
            save_name=f"{_dir}/competition_length_scales{mod}.pdf"
            fig.savefig(save_name,format="pdf",bbox_inches='tight',transparent=True)
            plt.show()

def addSingleChainToPhaseDiagram(ax,lps):
    """
    Add the critical buckling of a single chain to the plot. Changed this to
    p_l = 0.5 as this would better match the way in which I assigned clusters.
    """
    # single_chain_fit=lambda p: ((2*(1+np.log(0.01)/np.log(p)))**(1/0.277))/4e6
    single_chain_fit=lambda p: (3.5**2)*2e-4*(((1+np.log(0.5)/np.log(p)))**(4))
    ax.plot(lps,single_chain_fit(lps),'w-.',
        # label=r'$p_\ell^{N_c-1}=0.01$'
        label="Single chain buckling"
        )

def addCollectiveToPhaseDiagram(ax,lps):
    kappa=0.1
    ax.plot(lps,kappa/lps,'w:',lw=2,
            # label=r"$p_\ell\sim \frac{1}{\beta}$"
            label=r"Collective buckling"
            )

    # Rc=50
    # alpha=1e-4
    # kappa=(np.pi*Rc**2/8) * np.sqrt( (2e-4*4.22*alpha)/(2*3.75) )
    # print(f"{kappa=}")
    # quit()
    # ax.plot(lps,kappa/lps,'y:',lw=2,label=r"$p_\ell=\frac{\pi R_c^2}{8\beta}\sqrt{\frac{\mu K_{F0}\alpha}{2 \langle \ell \rangle}}$")
    #
    # kappa=np.pi**3*np.sqrt((4.22**3)/(2*2e-4*3.75**3))
    # ax.plot(lps,kappa/lps,'g-.',lw=2,label=r"$p_\ell=\frac{\pi^3}{\beta}\sqrt{\frac{K_{F0}^3}{2\mu l^3}}$")

def addPhaseBoundaries(ax,mb=0.25,bb=0.1):
    ax.axvline(mb,linestyle='--',color='w',lw=2,label='Microdomains boundary')
    ax.axhline(bb,linestyle=(0, (3, 1, 1, 1, 1, 1)),color='w',lw=2,label='Buckling boundary')
