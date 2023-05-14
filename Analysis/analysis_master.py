"""
    Run through all generated data and obtain quantitative information about
    colony formation.

    At the moment, different simulation runs are assumed to be labelled run###
    within the provided directory.

    Otherwise, any directory inside the GeneratedOutput folder may be accessed,
    but at present not sweeped.
"""

# Add parent to path
import sys

sys.path.insert(0, '..')

# Standard libraries
import numpy as np
import numpy.linalg as linalg
import pandas as pd
import argparse
import os
from glob import glob
import subprocess
import multiprocessing as mp
from joblib import Parallel, delayed
from copy import copy
import re
from pprint import pprint
from functools import reduce

# Third party
from hidden_prints import HiddenPrints
from tqdm import tqdm

# Custom modules
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm
from visualise_biofilm import setFigSize,setClassStaticMembers
from phage_vis_class import Phage
from infected_vis_class import Infected
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from AG43RodShapedBacteria import AG43RodShapedBacterium
from SphericalBacteria import SphericalBacterium
from RodShapedBacteria import setProjection  # choose plane to visualise
from PlotAnalysisResults import AnalysisPlotter,loadClusters,saveClusters
import utilities as ut
ut.setMPL()

# Remove later
import matplotlib.pyplot as plt

var_map={
    'LinkingProb'    : { "label" : r'$p_\ell$'},
    'Kappa'          : { "label" : r'$\kappa$'},
    'RodGrowthRate'  : { "label" : r'$g$' },
    'RodAspectRatio' : { "label" : r'$\ell_d$' },
    'adhKapp'        : { "label" : r'$\kappa$' },
    'adhDist'        : { "label" : r'$\delta$' }
    }

def setUpAP(dirs,analysis_dir,stem=None):
    """
    Better would be find all files then filter for stem
    """
    analysis_dir=f"{analysis_dir}/"
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    use_files=[]
    for dir in dirs:
        # print(f"{dir}/**/{stem}*")
        files=glob(f"{dir}/**/final*")
        if files:
            use_files.extend(files)

    if not use_files:
        print("Error, found no files!")
        quit()

    # process files sequentially
    hps={}
    for ii,file in enumerate(use_files):
        dir=os.path.split(file)[0]
        log_file=glob(f"{dir}/*log")[0]
        hp=pd.read_csv(log_file,sep="\t").to_dict('list')
        if ii==0:
            keys=[kk for kk in hp.keys()]
        for key,val in hp.items():
            hp[key]=val[0]
        run,repeat,tstem=AnalysisPlotter.getRunRepeatStem(file)
        hps[run]=hp

    use_files=[]
    for dir in dirs:
        print(dir)
        if stem!=None:
            files=glob(f"{dir}/**/{stem}*.dat")
        else:
            files=glob(f"{dir}/**/*.dat")
        if files:
            use_files.extend(files)

    AP=AnalysisPlotter(in_files=use_files,
                       analysis_dir=analysis_dir,
                       hps=hps,class_dict=class_dict)
    return AP

# Default class dictionary
class_dict={
   "RodShaped"         : RodShapedBacterium,
   "Spherical"         : RodShapedBacterium,
   "ChainingRodShaped" : ChainingRodShapedBacterium,
   "AG43RodShaped"     : AG43RodShapedBacterium
   }

def getHyperParams(dirs,stem):
    """
    Parameters:
        dirs: list of str
            list of directories to process
        stem: str
            file stem to load from
    Returns:
        hps: dict of dict
            dictionary of the hyperparameters with keys matching run id
    """
    hps={}
    for dir in dirs:
        repeats = glob(f"{dir}/repeat*")
        if not repeats:
            repeats=[dir]

        # Get log file
        try:
            log_file = f"{repeats[0]}/{stem}.log"
            hp=pd.read_csv(log_file,sep="\t").to_dict('list')
            for key,val in hp.items():
                hp[key]=val[0]
            # if "Rd" in hp:
            #     RodShapedBacterium.depletion_radius=hp['Rd'][0]
            run,_,_ = AnalysisPlotter.getRunRepeatStem(dir)
            hps[run]=hp
        except Exception as e:
            print(e)
            quit()
    return hps

def runComparison(in_files,
                  analysis_dir,
                  hps,
                  var,
                  class_dict={
                      "RodShaped"         : RodShapedBacterium,
                      "Spherical"         : RodShapedBacterium,
                      "ChainingRodShaped" : ChainingRodShapedBacterium
                      },
                  show=False):
    """
        Compare quantitative desciptors of the input colony data at the same
        colony size

        These include:
            Radial density function

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
            show: bool
                Flag to switch on displaying figures generated from analysis
    """
    if not os.path.exists(analysis_dir):
        print(f"Creating {analysis_dir}")
        os.makedirs(analysis_dir)
    else:
        print(f"{analysis_dir} already exists")

    if len(in_files)>1:
        in_files.sort(key=lambda x: int(re.search(r"(?<=run)(\d+)",x).group(1)))

    AP=AnalysisPlotter(in_files=in_files,
                       analysis_dir=analysis_dir,
                       hps=hps,class_dict=class_dict)

    AP.make_output_figures=True
    AP.generate_data_only=False

    results={}

    # try:
    # aps = AP.plotAspectRatio(var,False)

    # except Exception as e:
    #     print(e)

    # This is could be changed to two plot general functions in future
    # densities = AP.plotDensity(var,show)
    # results["densities"]=densities
    # print("Completed density analysis")
    show=True
    # wavelength=lambda a,b: AP.getWavelength(a,b)[0]
    # wavlen = AP.plotGeneral(fn=wavelength,
    #                         var=var,
    #                         ylabel="$\lambda$",
    #                         name="global_wavelength",
    #                         show=show)
    # results["wavelength"]=wavlen
    # print("Completed persistence length analysis")
    #
    # perlen = AP.plotGeneral(fn=AP.getPersistenceLength,
    #                         var=var,
    #                         ylabel="$\ell_p$",
    #                         name="persistence_length",
    #                         show=show)
    # results["persistence_length"]=perlen
    # print("Completed persistence length analysis")
    #
    # thetas = AP.plotGeneral(fn=AP.getLinkAngle,
    #                         var=var,
    #                         ylabel=r"$\langle \theta^2 \rangle$",
    #                         name="link_angle",
    #                         show=show)
    # results["link_angle"]=thetas
    # print("Completed link angle analysis")
    # quit()

    # pc = AP.plotGeneral(fn=AP.getPairCorrelation,var=var,
    #                    ylabel="$\ell_c$",
    #                    name="pair_correlation",
    #                    show=show)
    # results["pair_correlation"]=pc
    # print("Completed pair correlation analysis")
    #
    # wpc = AP.plotGeneral(fn=AP.getWeightedPairCorrelation,var=var,
    #                    ylabel=r"$\tilde{\ell}_c$",
    #                    name="weighted_pair_correlation",
    #                    show=show)
    # results["weighted_pair_correlation"]=wpc
    # print("Completed pair correlation analysis")
    #
    # getL1=lambda a,b: AP.getCurvatures(a,b)[0]
    # max_radius_curvatures_L1 = AP.plotGeneral(fn=getL1,var=var,
    #                                        ylabel="$\mathcal{R}$",
    #                                        name="curvature_L1",
    #                                        show=show)
    # results["max_radius_curvature_L1"]=max_radius_curvatures_L1
    # print("Completed curvature analysis L1")

    # getL2=lambda a,b: AP.getCurvatures(a,b)[1]
    # max_radius_curvatures_L2 = AP.plotGeneral(fn=getL2,var=var,
    #                                        ylabel="$\mathcal{R}$",
    #                                        name="curvature_L2",
    #                                        show=show)
    # results["max_radius_curvature_L2"]=max_radius_curvatures_L2
    # print("Completed curvature analysis L2")
    #
    # if len(in_files)>1:
    #     AP.plotLengthScale(wpc,max_radius_curvatures_L1,mod='_L1_wpc',
    #                        label=r"$\tilde{\ell}_c/\mathcal{R}$")
    #     # AP.plotLengthScale(cluster_analysis,max_radius_curvatures_L1,'_L1')
    #     # AP.plotLengthScale(cluster_analysis,max_radius_curvatures_L2,'_L2')
    #
    # quit()
    #
    thresh_th=3
    mean_cluster_3=lambda a,b: AP.getMeanClusterArea(a,b,thresh_th)
    # mca = AP.plotGeneral(fn=mean_cluster_3,var=var,
    #                      ylabel=r"$\langle A_3 \rangle$",
    #                      name=f"cluster_mean_3",
    #                      show=show)
    # S = AP.plotGeneral(fn=AP.getGlobalOrder,var=var,
    #                    ylabel="$\mathcal{S}$",
    #                    name="S",
    #                    show=show)
    # results["S"]=S
    # print("Completed order analysis")

    # mca = AP.plotGeneral(fn=mean_cluster_3,var=var,
    #                      ylabel=r"$\langle A_5 \rangle$",
    #                      name=f"cluster_mean_5",
    #                      show=show)
    # mca = AP.plotGeneral(fn=mean_cluster_3,var=var,
    #                      ylabel=r"$\sqrt{\langle A_5 \rangle}$",
    #                      name="sqrt_cluster_mean_5",
    #                      show=show)

    thresh_th=3
    mean_cluster_3=lambda a,b: np.sqrt(AP.getMeanClusterArea(a,b,thresh_th))
    mca = AP.plotGeneral(fn=mean_cluster_3,var=var,
                         ylabel=r"$\sqrt{\langle A_3 \rangle}$",
                         name="sqrt_cluster_mean_3_mod",
                         show=show)

    thresh_th=5
    mean_cluster_5=lambda a,b: np.sqrt(AP.getMeanClusterArea(a,b,thresh_th))
    mca = AP.plotGeneral(fn=mean_cluster_5,var=var,
                         ylabel=r"$\sqrt{\langle A_5 \rangle}$",
                         name="sqrt_cluster_mean_5_mod",
                         show=show)

    roughness = AP.plotGeneral(fn=AP.getRoughness,var=var,
                               ylabel=r"$\nu$",
                               name="roughness",
                               show=show)

    # cluster_analysis = AP.plotClusterAnalysis(var,show)
    # results["clusters"]=cluster_analysis
    # print("Completed cluster analysis")

    # buckles = AP.plotGeneral(fn=AP.getBuckled,var=var,
    #                          ylabel="buckled",
    #                          name="buckles",
    #                          show=show)
    # results["buckles"]=buckles
    # print("Completed buckled analysis")

    asps = AP.plotGeneral(fn=AP.getAspectRatio,var=var,
                          ylabel=r"$\alpha$",
                          name="aspect_ratio",
                          show=show)
    results["aspect_ratios"]=asps
    print("Completed aspect ratio analysis")

    density = AP.plotGeneral(fn=AP.getDensity,var=var,
                             ylabel=r"$\rho$",
                             name="density",
                             show=show)
    results["density"]=density
    print("Completed fractal dimension analysis")

    # fractal_dims = AP.plotGeneral(fn=AP.getFractalDimension,var=var,
    #                               ylabel="$\mathcal{D}$",
    #                               name="fractal_dimension",
    #                               show=show)
    # results["fractal_dims"]=fractal_dims
    # print("Completed fractal dimension analysis")

    max_overlaps = AP.plotGeneral(fn=AP.getMaxOverlap,var=var,
                                  ylabel="$h$",
                                  name="overlaps",
                                  show=show)
    results["max_overlaps"]=max_overlaps
    print("Completed overlap analysis")


    # except Exception as e:
    #     print(e)

    # try:

    # except Exception as e:
    #     print(e)

    # try:
    # AP.plotChainDescriptors(show)
    # except Exception as e:
    #     print(e)
    # AP.checkContours(show)

    # AP.plotRadialDist(show)
    quit()
    return results

def getFileSize(fname,class_dict=class_dict):
    data = pd.read_csv(fname,sep="\t")
    (cells,
     species_load_dict) = intialiseElementsFromData(data,class_dict)
    return AnalysisPlotter.getColonySize(cells)

def getFilesAtSizes(dir,stem,sizes):
    """
        For files in dir, find all which are closest to the set of sizes

        Parameters:
            dir: str
                path to files
            stem: str
                stem of filenames
            sizes: (N,) floats, np array
                list of sizes to find

        Returns:
            files: list of str
                those files which were closest

    """
    selected_files=[]
    files=glob(f"{dir}/{stem}*.dat")
    files.sort(
        key=lambda x: int(re.search(fr"(?<={args.stem}\_)(\d+)",x).group(1))
        )
    file_size=[]
    for file in tqdm(files):
        file_size.append([file,getFileSize(file)])

    for ss in tqdm(sizes):
        file=min(file_size, key=lambda x: abs(x[1]-ss) )
        if (abs(file[1]-ss)/ss)<1e-1:
            selected_files.append(file[0])
        else:
            selected_files.append(None)
            print(f"Following file is out of tolerance of {ss}")
            print(file)
    return selected_files

def getAllFilesAtSizes(sizes,dirs,analysis_dir,stem):
    """
        Find all the files which are at the

        Parameters
            sizes: array of float
                sizes to determine the file names at
            dirs: list of str
                directories to search for files in
            analysis_dir: str
                directory to load/save files at the correct size to
            stem: str
                file stems to load from

        Returns:
            selected_files: list of str
                all files at the required sizes (within 10%)
    """
    sf_fname=f"{analysis_dir}/selected_files.csv"
    selected_files={}
    if os.path.exists(sf_fname):
        sf=pd.read_csv(sf_fname,dtype=str)
        selected_files=sf.to_dict('list')
        print(f"Loaded following files from {sf_fname}")

    print(selected_files.keys())
    for ss in sizes:
        if f's_{ss}' not in selected_files.keys():
            selected_files[f's_{ss}']=[]
            print(ss)
    prev_files=[vv for sl in selected_files.values() for vv in sl]
    for dir in dirs:
        repeats = glob(f"{dir}/repeat*")
        if not repeats:
            repeats=[dir]
        for repeat in repeats:
            already_processed=len([vv for vv in prev_files if repeat in str(vv)])
            if already_processed>len(sizes):
                print(f"Error! There are {already_processed}/{len(sizes)}")
                quit()
            elif already_processed<len(sizes):
                print(f"Processing {repeat}")
                files=getFilesAtSizes(repeat,stem,sizes)
                for file,ss in zip(files,sizes):
                    selected_files[f's_{ss}'].append(file)
                df=pd.DataFrame.from_dict(selected_files)
                df.to_csv(sf_fname,index=False)
                print(f"Updated {sf_fname} with")
                pprint(files)
    return selected_files

def analyseEvolutionAtSizes(args,sizes=None):
    """
    Parameters:
        args: namespace of input parameters
            arguments parsed from main
    """
    # Data files
    dirs=glob(f"{args.base_path}/{args.data_dir}/**/run*",recursive=True)

    max_size=37500 # Size of colony no simulation will exceed

    # Compare outputs at these colony sizes to evaluate the colony evolution
    if sizes==None:
        sizes=np.array([*np.arange(max_size/4,max_size,max_size/4), max_size])

    # Save the files to this directory
    analysis_dir=f"{args.base_path}/{args.analysis_dir}/"
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    # each size will be a key corresponding to another dictionary with the
    # lps and densities at that size - same for aspect ratios aps
    densities={}
    aps={}

    # Attempt to load precalculated densities at all sizes we need
    sizes_to_be_calculated=[]
    for ss in sizes:
        densities_file=f"{analysis_dir}/densities_{ss:3.3f}.csv"
        if os.path.exists(densities_file):
            dens=np.genfromtxt(densities_file,
                               delimiter='\t')
            densities[f's_{ss}']={}
            densities[f's_{ss}']['lps']=dens[:,0]
            densities[f's_{ss}']['densities']=dens[:,1]
            densities[f's_{ss}']['err_densities']=dens[:,2]
        else:
            sizes_to_be_calculated.append(ss)

        aspect_file=f"{analysis_dir}/aps_{ss:3.3f}.csv"
        if os.path.exists(aspect_file):
            ap=np.genfromtxt(aspect_file,
                               delimiter='\t')
            aps[f's_{ss}']={}
            aps[f's_{ss}']['lps']=ap[:,0]
            aps[f's_{ss}']['aps']=ap[:,1]
            aps[f's_{ss}']['err_aps']=ap[:,2]
        else:
            if ss not in sizes_to_be_calculated:
                sizes_to_be_calculated.append(ss)

    # Everything went as planned and all densities we need exist
    if not sizes_to_be_calculated:
        AnalysisPlotter.plotEvolutionQuantity(densities,analysis_dir,
                                              qname="densities",
                                              poster=True)
        AnalysisPlotter.plotEvolutionQuantity(aps,analysis_dir,qname="aps",
                                              label="Aspect Ratios",
                                              poster=True)
        return

    print("Need to calculate densities for the following sizes")
    pprint(sizes_to_be_calculated)

    # Get the hyperparameters
    hps=getHyperParams(dirs,stem=args.stem)

    # load files for which the sizes have already been calculated
    selected_files=getAllFilesAtSizes(
        sizes,
        dirs,
        analysis_dir,
        args.stem
        )

    for ss in sizes_to_be_calculated:
        try:
            files_at_current_size=selected_files[ss]
        except Exception as e:
            files_at_current_size=selected_files[f's_{ss}']

        # for ff in files_at_current_size:
        #     curr_size=getFileSize(ff)
        #     if np.allclose(curr_size,ss)

        tmp_aps,tmp_dens = runComparison(
                              files_at_current_size,
                              analysis_dir,
                              hps,
                              class_dict={
                                  "RodShaped"         : RodShapedBacterium,
                                  "Spherical"         : RodShapedBacterium,
                                  "ChainingRodShaped" : ChainingRodShapedBacterium
                                  },
                              show=False)
        densities[f's_{ss}']={}
        densities[f's_{ss}']['lps']=tmp_dens[:,0]
        densities[f's_{ss}']['densities']=tmp_dens[:,1]
        densities[f's_{ss}']['err_densities']=tmp_dens[:,2]
        np.savetxt(f"{analysis_dir}/densities_{ss:3.3f}.csv",tmp_dens,
                    delimiter='\t')
        aps[f's_{ss}']={}
        aps[f's_{ss}']['lps']=tmp_aps[:,0]
        aps[f's_{ss}']['aps']=tmp_aps[:,1]
        aps[f's_{ss}']['err_aps']=tmp_aps[:,2]
        np.savetxt(f"{analysis_dir}/aps_{ss:3.3f}.csv",tmp_aps,
                    delimiter='\t')
    pprint(densities)
    pprint(aps)

    AnalysisPlotter.plotEvolutionQuantity(densities,analysis_dir,
                                          qname="densities",
                                          poster=True)
    AnalysisPlotter.plotEvolutionQuantity(aps,analysis_dir,qname="aps",
                                          label="Aspect Ratios",
                                          poster=True)

def analyseQuantityEvolution(args,var=None):
    """
    Find the evolution of the energy for now

    Parameters:
        args: namespace of input parameters
            arguments parsed from main
    """
    bend_only=False
    # Data files
    # dirs=glob(f"{args.base_path}/{args.data_dir}/**/repeat*/*",recursive=True)
    print("This will only work for a single repeat at the moment")
    log_file = f"{args.base_path}/{args.data_dir}/{args.stem}.log"
    setClassStaticMembers(log_file)
    hp=pd.read_csv(log_file,sep="\t").to_dict('list')
    lp=hp['LinkingProb'][0]
    analysis_dir=f"{args.base_path}/{args.analysis_dir}/Energy/"

    if args.end_index==-1:
        end=len(glob(f"{args.base_path}/{args.data_dir}/{args.stem}_*.dat"))-1
    else:
        end=args.end_index
    par_args_list = [ f"{args.base_path}/{args.data_dir}/{args.stem}_{ii:05d}.dat"
                      for ii in range(end + 1)
                      ]
    try:
        energies=np.loadtxt(f"{analysis_dir}/energy_{bend_only=}.txt")
    except Exception as e:
        print(e)
        class_dict={
           "RodShaped"         : RodShapedBacterium,
           "Spherical"         : RodShapedBacterium,
           "ChainingRodShaped" : ChainingRodShapedBacterium
           }
        contact_energies=[]
        spring_energies=[]
        sizes=[]
        for file in par_args_list:
            cells, species_load_dict = intialiseElementsFromData(
                pd.read_csv(file,sep='\t'),
                class_dict
                )
            N=len(cells)
            hertzian_energy=RodShapedBacterium.getContactEnergies(cells)
            spring_energy=ChainingRodShapedBacterium.getSpringEnergies(cells,bend_only)
            contact_energies.append(np.sum(hertzian_energy)/N)
            spring_energies.append(np.sum(spring_energy)/N)
            sizes.append(getFileSize(file))

        # Save the files to this directory
        if not os.path.exists(analysis_dir):
            os.makedirs(analysis_dir)

        energies=np.array([sizes,contact_energies,spring_energies]).T
        print(energies)
        np.savetxt(f"{analysis_dir}/energy_{bend_only=}.txt",energies)
    print(energies.shape)

    if bend_only:
        sp_label=r"$E_{\beta}$"
    else:
        sp_label=r"$E_{\beta}+E_{\kappa}$"
    fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
    ax.plot(energies[:,0],energies[:,1]/np.max(energies[:,1]),label="$E_{Y}$")
    ax.plot(energies[:,0],energies[:,2]/np.max(energies[:,2]),label=sp_label)
    ax.set_ylabel("$E/E_{max}$")
    ax.set_xlabel("Colony size")
    # ax.set_yscale("log")
    plt.legend()
    fig.savefig(f"{analysis_dir}/energy_{lp=}_{bend_only}.pdf",format="pdf",bbox_inches='tight',
                transparent=True)
    plt.show()

    fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
    indces=energies[:,1]>0
    ax.plot(energies[indces,0],energies[indces,2]/energies[indces,1])

    ax.set_ylabel(sp_label+"$/E_{Y}$")
    ax.set_xlabel("Colony size")
    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.savefig(f"{analysis_dir}/energy_{lp=}_{bend_only}_ratio.pdf",format="pdf",bbox_inches='tight',
                transparent=True)
    plt.show()
    # Parallel(n_jobs=6, verbose=10, backend="multiprocessing")(
    #     delayed(find_energies)(**par_args)
    #     for par_args in par_args_list
    #     )

    # Save the files to this directory
    # analysis_dir=f"{args.base_path}/{args.analysis_dir}/Comparison/"
    # if not os.path.exists(analysis_dir):
    #     os.makedirs(analysis_dir)
    #
    # qnty={}

    quit()

def analyseEvolutionAll(dirs,analysis_dir,args,var_map=var_map):
    """
        Make plots of the evolution of all quantities as a function of colony
        size
    """
    AP=setUpAP(dirs,analysis_dir)
    fn_dict={
        'fractal_dimension' : AP.getFractalDimension,
        'aspect_ratio'      : AP.getAspectRatio,
        'density'           : AP.getDensity,
        'roughness'         : AP.getRoughness,
        'S'                 : AP.getGlobalOrder
        }
    df=AP.getGeneralData(fn_dict,cores=args.max_cores)
    quit()
    labels=[r'$\mathcal{D}$',r'$\alpha$',r'$\rho$']
    pprint(df)
    base_state={
        'RodAspectRatio' : 4.0,
        'RodGrowthRate'  : 0.0002,
        'adhKapp'        : 0.01,
        'adhDist'        : 0.1
        }
    base_case = df.loc[
        (df['RodAspectRatio'] == base_state['RodAspectRatio']) &
        (df['RodGrowthRate']  == base_state['RodGrowthRate'])  &
        (df['adhKapp']        == base_state['adhKapp'])        &
        (df['adhDist']        == base_state['adhDist'])
        ]
    df=df.loc[
        df['RodAspectRatio'].isin(np.array([3,4,5],dtype=int))
        ]
    """
        Plot evolutions of all quantities
    """

    print(base_case)
    keys=[ kk for kk in base_state.keys() ]
    for key,val in base_state.items():
        other_keys=[kk for kk in keys if kk != key]
        select=pd.DataFrame.copy(
            df.loc[df[key]!=base_state[key]]
            )
        select=select.append(base_case)
        print(f"{key} has {select[key].nunique()} unique entries")
        avg=select.groupby([key,'stem']).mean()
        avg.sort_values(by=[key,'stem'],inplace=True)
        pprint(avg)
        plt.close()
        figsize=setFigSize(247)
        fig,ax=plt.subplots(1,len(names),figsize=[len(names)*figsize[0],figsize[1]])
        for iv,new_df in avg.groupby(level=0):
            for ii in range(len(names)):
                ax[ii].plot(new_df.loc[:,'size'],
                            new_df.loc[:,names[ii]],
                            label=f'{key}={iv}')
                ax[ii].set_xlabel("size")
                ax[ii].set_ylabel(labels[ii])
                ax[ii].set_yscale('log')
        plt.legend()

        # add=lambda t,s: s+t
        # fname_stem=reduce(add,[f'{v}_' for v in names])[:-1]
        comp_dir=f"{analysis_dir}/all_data/figures"
        os.makedirs(comp_dir,exist_ok=True)
        fig_name=f"{comp_dir}/evolution_{key}.pdf"

        fig.savefig(fig_name,
                    format="pdf",
                    bbox_inches='tight',
                    transparent=True)
        plt.show()

    """
        Plot final state of all quantities
    """
    keys=[ kk for kk in base_state.keys() ]
    for key,val in base_state.items():
        other_keys=[kk for kk in keys if kk != key]
        select=pd.DataFrame.copy(
            df.loc[ ( df[key]!=base_state[key] ) &
                    ( df['stem'].str.contains('final') ) ]
            )
        select.append(base_case)
        print(f"{key} has {select[key].nunique()} unique entries")
        avg=select.groupby('run').mean()
        avg.sort_values(by=key,inplace=True)
        pprint(avg)
        plt.close()
        fig,ax=plt.subplots(1,len(names))
        for ii in range(len(names)):
            ax[ii].plot(avg.loc[:,key],avg.loc[:,names[ii]])
            ax[ii].set_xlabel(var_map[key]['label'])
            ax[ii].set_ylabel(labels[ii])
        plt.show()
    quit()

    """
        Phase plots of all quantities
    """

def analyseSingleTime(args,var,dirs):
    in_files=[]
    hps={}
    for dir in dirs:
        if len(dirs)>1:
            run=re.search(r"run\d+",dir)[0]
        else:
            run=""
        # print(f"processing {dir}")
        repeats = glob(f"{dir}/repeat*")
        if not repeats:
            repeats=[dir]
        num_repeats=len(repeats)

        try:
            log_file = f"{repeats[0]}/{args.stem}.log"
            hp=pd.read_csv(log_file,sep="\t").to_dict('list')
            for key,val in hp.items():
                hp[key]=val[0]
            # if "Rd" in hp:
            #     RodShapedBacterium.depletion_radius=hp['Rd'][0]
            hps[run]=hp
        except Exception as e:
            print(e)
            quit()

        if args.phage:
            class_dict["Phage"]=Phage

        for repeat in repeats:
            if args.start_index==-1:
                start_index=len(glob(f"{repeat}/{args.stem}_*.dat"))
                in_file = f"{repeat}/final_{start_index:05d}.dat"
            else:
                start_index=args.start_index
                in_file = f"{repeat}/{args.stem}_{start_index:05d}.dat"

            if os.path.exists(in_file):
                print(f"exists {in_file=}")
                in_files.append(in_file)
            else:
                print(f"{in_file} not found, skipping...")
                continue

            rep_num = re.search(r"repeat\d+",repeat)
            if rep_num:
                analysis_dir=f"{args.base_path}/{args.analysis_dir}/{run}/{rep_num[0]}"
            else:
                analysis_dir=f"{args.base_path}/{args.analysis_dir}/{run}"
            if args.compare_only==False:
                print(f"Reading from output file {in_file}")
                print(f"Writing to analysis directory {analysis_dir}")
                # try:
                runComparison(
                    in_files=[in_file],
                    analysis_dir=analysis_dir,
                    hps=hp,
                    var=var,
                    class_dict=class_dict,
                    show=True
                    )
                # except Exception as e:
                #     print(e)

    if len(dirs)>1:
        print("=== Begin Comparisons ===")
        print(f"Reading from output files:")
        pprint(in_files)
        analysis_dir=f"{args.base_path}/{args.analysis_dir}/"
        print(f"Writing to analysis directory {analysis_dir}")
        runComparison(
            in_files=in_files,
            hps=hps,
            var=var,
            analysis_dir=analysis_dir,
            class_dict=class_dict,
            show=False
            )

def makeBucklingPhaseDiagram(dirs,variables,analysis_dir,args):
    AP=setUpAP(dirs,analysis_dir,stem='final')
    # AP.makeCollagePlot(variables=variables)
    # quit()

    # getLPrll=lambda a,b: AP.getColonyCorrelations(a,b)[0]
    # getLPerp=lambda a,b: AP.getColonyCorrelations(a,b)[1]

    # fns=[AP.getFractalDimension,
    #      AP.getAspectRatio,
    #      AP.getDensity,
    #      AP.getRoughness,
    #      getLPrll,
    #      getLPerp]
    # names=['fractal_dimension','aspect_ratio','density',
    #        'roughness',"c_par","c_perp"]
    # labels=[r'$\mathcal{D}$',r'$\alpha$',r'$\rho$',r'$\nu$',
    #         r'\ell^\par_c',r'\ell^\perp_c']
    # df=AP.getGeneralData(fns,names,cores=args.max_cores)

    mean_cluster_3=lambda a,b: AP.getMeanClusterArea(a,b,thresh_th=3)
    mean_cluster_5=lambda a,b: AP.getMeanClusterArea(a,b,thresh_th=5)

    getL1=lambda a,b: AP.getCurvatures(a,b)[0]
    getL2=lambda a,b: AP.getCurvatures(a,b)[1]

    wavelength=lambda a,b: AP.getWavelength(a,b)[0]

    labels=[r'$\langle \theta^2 \rangle$',r'$\mathcal{D}$',r'$\alpha$',r'$\rho$',
            r'$\nu$',"$S$",r"$\mathcal{R}$",r"$\mathcal{R}$",r"$\tilde{\ell}_c$",
            r"$\ell_c$",r"$\langle A_3 \rangle$",r"$\langle A_5 \rangle$",
            r"$\ell_p$"]
    fn_dict={
        'link_angle'                 : AP.getLinkAngle,
        'fractal_dimension'          : AP.getFractalDimension,
        'aspect_ratio'               : AP.getAspectRatio,
        'density'                    : AP.getDensity,
        'roughness'                  : AP.getRoughness,
        'S'                          : AP.getGlobalOrder,
        'curvature_L1'               : getL1,
        'curvature_L2'               : getL2,
        # 'c_par'             : getLPrll,
        'weighted_pair_correlation'  : AP.getWeightedPairCorrelation,
        'pair_correlation'           : AP.getPairCorrelation,
        'cluster_mean_3'             : mean_cluster_3,
        'cluster_mean_5'             : mean_cluster_5,
        'persistence_length'         : AP.getPersistenceLength,
        # 'wavelength'                 : wavelength
        }



    df=AP.getGeneralData(fn_dict,cores=args.max_cores)

    # AP.findCollectiveBucklingPoint(df,variables=variables,
    #                                keys=['persistence_length',
    #                                      'pair_correlation']
    #                                )
    for f,n,l in zip(fn_dict.values(),fn_dict.keys(),labels):
        # local_dict={n:f}
        break
        if n=="cluster_mean_5":
            AP.phasePlot(f,n,variables,lbl=l)
            plt.show()
            plt.close()

    AP.clusterColonies(variables=variables)
    AP.assignBuckles(variables=variables,full=True)
    quit()

    AP.makeCollagePlot(variables=variables)

    quit()
    # AP.assignBuckles(variables=variables)

    quit()

    names=[]
    name="roughness"
    AP.phasePlot(AP.getRoughness,name,variables)
    plt.show()
    plt.close()

    # name="area_fraction"
    # AP.phasePlot(AP.getColonyAreaFraction,name,variables)
    # plt.show()
    # plt.close()

    name="density"
    AP.phasePlot(AP.getDensity,name,variables,lbl=r"$\rho$")
    plt.show()
    plt.close()

    name="aspect_ratio"
    AP.phasePlot(AP.getAspectRatio,name,variables,lbl=r"$\alpha$")
    plt.show()
    plt.close()

    name="fractal_dimension"
    AP.phasePlot(AP.getFractalDimension,name,variables,lbl=r"$\mathcal{D}$")
    plt.show()
    plt.close()

    # name="curvature_L1"
    # getL1=lambda a,b: AP.getCurvatures(a,b)[0]
    # AP.phasePlot(getL1,name,variables)
    # plt.show()
    # plt.close()

    quit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='visualise_phage.py',
        usage='%(prog)s [options] path',
        description='Visualisation script for infected biofilms',
        epilog='Author: Rory Claydon'
    )
    parser.add_argument("-S", "--start-index", dest="start_index",
                        type=int,
                        required=False,
                        help="load data from"
                             + " SimOutput/output_{start_index:05d}.dat")
    parser.add_argument("-E", "--end-index",dest="end_index",
                        type=int,
                        required=False,
                        default=None,
                        help="load data up until"
                             + " SimOutput/output_{start_index:05d}.dat"
                             + " if end=-1 use all in directory")
    parser.add_argument("-P","--phage",action="store_true",required=False)
    parser.add_argument("--step", type=int, required=False,
                        default=1,
                        help="skip step number outputs")
    parser.add_argument("-I", "--in-file", dest="in_file",
                        type=str,
                        required=False,
                        default=None,
                        help="output data filename if different from default")
    parser.add_argument("--max-cores", type=int, default=mp.cpu_count(),
                        required=False, action="store", dest="max_cores",
                        help="maximum number of cores to use to produce output")
    parser.add_argument("-BD", "--base-dir", type=str, required=False, default="../GeneratedOutput",
                        action="store", dest="base_path",
                        help="default directory to look for inputs")
    parser.add_argument("-DD", "--data-dir", type=str, required=False, default="SimOutput",
                        action="store", dest="data_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("-AD", "--analysis-dir", type=str, required=False, default="AnalysisResults",
                        action="store", dest="analysis_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("--stem", type=str, required=False, default="biofilm",
                        action="store", dest="stem",
                        help="default stem for input file names i.e. path/{stem}_%05d.dat")
    parser.add_argument("-C","--compare-only",action="store_true",required=False,dest="compare_only")
    parser.add_argument("--phase-diagram",action="store_true",
                        required=False,dest="phase_diagram")
    args = parser.parse_args()

    if not args.start_index and args.phase_diagram:
        variables=[
            { 'iv': "LinkingProb", 'label' : "$p_\ell$"},
            { 'iv': "BendRig", 'label' : r"$\beta$"}
            ]
        dirs=glob(f"{args.base_path}/{args.data_dir}/**/run*",recursive=True)
        makeBucklingPhaseDiagram(dirs,variables,f"{args.base_path}/{args.analysis_dir}/",args)
        analyseEvolutionAll(dirs,f"{args.base_path}/{args.analysis_dir}/",args)
        quit()
        analyseQuantityEvolution(args)
        analyseEvolutionAtSizes(args)

    if args.in_file:
        print("Specific in file not supported yet")
        quit()

    dirs=glob(f"{args.base_path}/{args.data_dir}/**/run*",recursive=True)
    if not dirs:
        dirs=[f"{args.base_path}/{args.data_dir}"]
    # print("Found the following directories:")
    # pprint(dirs)

    if args.phase_diagram:
        variables=[
            { 'iv': "LinkingProb", 'label' : "$p_\ell$"},
            { 'iv': "BendRig", 'label' : r"$\beta$"}
            ]
        makeBucklingPhaseDiagram(dirs,variables,f"{args.base_path}/{args.analysis_dir}/")
        quit()

    if args.end_index is None:
        # var={ 'iv': "Kappa", 'label' : "$\kappa$" }
        # var={ 'iv': "RodGrowthRate", 'label' : "$g$" }
        var={ 'iv': "LinkingProb", 'label' : "$p_\ell$" }
        analyseSingleTime(args,var,dirs)
    else:
        print("not sure what to do here")
        quit()
        analyseEvolutionAtSizes()
