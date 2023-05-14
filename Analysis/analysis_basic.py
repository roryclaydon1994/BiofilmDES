"""

"""

import argparse
import multiprocessing as mp
from glob import glob
from analysis_master import analyseSingleTime, makeBucklingPhaseDiagram

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

    args.start_index=-1
    args.end_index=-1
    args.compare_only=True
    dirs=glob(f"{args.base_path}/{args.data_dir}/**/run*",recursive=True)
    variables=[
        { 'iv': "LinkingProb", 'label' : "$p_\ell$"},
        { 'iv': "BendRig", 'label' : r"$\beta$"}
        ]
    # makeBucklingPhaseDiagram(dirs,variables,f"{args.base_path}/{args.analysis_dir}/",args)
    analyseSingleTime(args,variables[0],dirs)
