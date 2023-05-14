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
from visualise_biofilm import setFigSize
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm

params = {'font.family': 'serif',
          'text.usetex': True,
          'axes.titlesize': 10,
          'axes.labelsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'legend.frameon': False,
          'lines.linewidth': 1.5,
          'legend.handlelength': 1.5,
          'legend.columnspacing': 0.75,
          'legend.handletextpad': 0.5,
          # 'legend.title_fontsize' : None,
          'legend.fontsize': 10,
          'font.size': 10,
          'legend.columnspacing': 1.0,
          "axes.formatter.useoffset": False,
          # "text.latex.preamble"   : r"\usepackage[detect-all,locale=DE]{siunitx},\usepackage[T1]{fontenc},\usepackage{amsmath}"
          }
mpl.rcParams.update(params)


class AnyObject(object):
    pass


class ColouredBacteria(object):
    def __init__(self, cell, colour):
        self.cell = cell
        self.cell.colour = colour

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        return self.cell.legend_artist(legend, orig_handle, fontsize, handlebox)


def orderData(x, y, z):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    inds = np.argsort(x)
    return x[inds], y[inds], z[inds]


def getStatistics(dir, quantity='lengths', ret_hps=False):
    """
        Get the buckling probability, avg critical length and sem of lengths
        Parameters:
            dir: str
                directory to look in for sim files
        Returns:
            lp: float
                linking probability
            buckling_prob: float
                chance of buckling
            av_l: float
                average critical length
            sem_l: float
                standard error of the mean for the critical length
    """
    hps = getHyperParams(dir)
    try:
        lengths = np.loadtxt(f"{dir}/{quantity}.dat")
    except Exception as e:
        try:
            lengths = np.loadtxt(f"{dir}/{quantity}.dat", skiprows=1)
        except Exception as e:
            print(e)
            raise

    lps = hps["LinkingProb"]
    probs = hps["NBuckles"] / hps["NTrials"]
    av_l = np.nan_to_num(lengths.mean())
    sem_l = np.nan_to_num(scipy.stats.sem(lengths))
    if ret_hps == True:
        return (lps, probs, av_l, sem_l), hps
    else:
        return (lps, probs, av_l, sem_l)


def getHyperParams(dir):
    """
        Parameters:
            dir: str
                directory to look in for log file
        Returns:
            hps: dict
                dictionary of the hyperparamters
    """
    try:
        hps = pd.read_csv(f"{dir}/criticalLength.log", sep="\t").to_dict("list")
    except Exception as e:
        try:
            hps = pd.read_csv(f"{dir}/LcDist.log", sep="\t").to_dict("list")
        except Exception as e:
            print(e)
            quit()

    for key in hps.keys():
        hps[key] = hps[key][0]

    return hps


def getCells(file):
    """
        Parameters:
            file: str
                file to load cell data from

        Returns:
            cells: list of ChainingRodShapedBacterium
                cells from this data file
    """
    class_dict = {"ChainingRodShaped": ChainingRodShapedBacterium}
    data = pd.read_csv(file, sep="\t")
    (cells,
     species_load_dict) = intialiseElementsFromData(data, class_dict)
    return cells


def averageNumberBeforeBuckling(dir):
    """
        Parameters:
            dir: str
                directory to search in for data files
        Returns:
            lp: float
                linking probability
            av_N: float
                average number of bacteria in the chain before it buckles
            sem_N: float
                standard error of the mean of the above
    """
    hps = getHyperParams(dir)

    Ns = []
    for file in glob(f"{dir}/buckled_true*"):
        cells = getCells(file)
        Ns.append(len(cells))

    return (hps["LinkingProb"],
            np.nan_to_num(np.mean(Ns)),
            np.nan_to_num(scipy.stats.sem(Ns)))


def getProbDistBuckling(data_dir, analysis_dir, lc, avN, show=True):
    """
        Plot the probability of buckling as a function of linking prob

        Parameters:
            data_dir: str
                directory to look in for trial dirs
            analysis_dir: str
                directory to save data to
            lc: float
                critical chain length determined from simulations
            avN: float
                average number of bacteria in a chain which buckles
            show: bool
                if true show plot
    """
    lps = []
    bps = []
    err_bps = []
    for dd in glob(f"{data_dir}/*"):
        try:
            lp, bp, av_l, sem_l = getStatistics(dd)
            lps.append(lp)
            bps.append(bp)
            hps = getHyperParams(dd)
            err_bps.append(1 / np.sqrt(hps["NTrials"]))
        except Exception as e:
            print(e)
            print(f"Failed for {dd}")

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    marker = "o"
    lps, bps, err_bps = orderData(lps, bps, err_bps)
    ax.errorbar(lps, bps, err_bps, marker=marker, mfc="w", mec='k',
                capsize=5, alpha=0.7, label="simulation", zorder=1)
    ps = np.linspace(lps.min(), lps.max())
    # ax.plot(ps,ps**(-1+lc/4.5),"r:",
    #         label=r"$p_\ell^{(\langle L \rangle_\mathrm{sim}/\langle\ell\rangle)-1}$",
    #         zorder=2
    #         )
    ax.plot(ps, ps ** (-1 + avN), "k--",
            label=r"$p_\ell^{N_c-1}$",
            zorder=3
            )
    lp_critical = 0.5 ** (1 / (-1 + avN))
    print(f"{lp_critical=}")
    ax.axvline(lp_critical, linestyle=":", color='k', zorder=4,
               label=r"$p_\ell^c$")  # =10^{-\frac{2}{(\langle L \rangle_\mathrm{sim}/\langle\ell\rangle)-1}}$")
    ax.axhline(0.5, linestyle=":", color='m', label=r"$\mathcal{P}_{b}^{c}$")
    save_name = f"{analysis_dir}/buckling_prob_vs_linking_prob.pdf"
    ax.set_ylabel("$\mathcal{P}_b$")
    ax.set_xlabel("$p_\ell$")
    plt.legend(ncol=2, facecolor='w', framealpha=0.7, frameon=False, edgecolor=None,
               loc='upper left',
               # handletextpad=0.6,borderpad=0.2,
               # columnspacing=0.5,handlelength=1.0
               )

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()
    print(f"Saving to {save_name}")


def getProbDistChainingLength(data_dir, analysis_dir, show=True):
    """
        Plot the distribution of crticial length as a function of linking prob

        Parameters:
            data_dir: str
                directory to look in for trial dirs
            analysis_dir: str
                directory to save data to
            show: bool
                if true show plot
    """
    lps = []
    avls = []
    semls = []
    for dd in glob(f"{data_dir}/*"):
        try:
            lp, bp, av_l, sem_l = getStatistics(dd)
            lps.append(lp)
            avls.append(av_l)
            semls.append(sem_l)
        except Exception as e:
            print(e)
            print(f"Failed for {dd}")

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    marker = "o"
    lps, avls, semls = orderData(lps, avls, semls)
    ax.errorbar(lps, avls, semls, marker=marker, mfc="w", mec='k',
                capsize=5, alpha=0.7, zorder=1)
    save_name = f"{analysis_dir}/crticial_length_vs_linking_prob.pdf"
    ax.set_ylabel("$L$")
    ax.set_xlabel("$p_\ell$")

    avlc = avls[avls > 0].mean()
    ax.axhline(avlc, linestyle='--', label=rf"$L_c={avlc:3.1f}$", color='k', zorder=2)
    ax.legend()

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()
    print(f"saved to {save_name}")
    return avlc


def getProbDistChainingTime(data_dir, analysis_dir, show=True):
    """
        Plot the distribution of crticial time as a function of linking prob

        Parameters:
            data_dir: str
                directory to look in for trial dirs
            analysis_dir: str
                directory to save data to
            show: bool
                if true show plot
    """
    lps = []
    avts = []
    semts = []
    gs = []
    for dd in glob(f"{data_dir}/*"):
        try:
            (lp, bp, av_t, sem_t), hps = getStatistics(dd, quantity='times', ret_hps=True)
            lps.append(lp)
            g = hps['RodGrowthRate']
            avts.append(g * av_t)
            semts.append(g * sem_t)
        except Exception as e:
            print(e)
            print(f"Failed for {dd}")

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    marker = "o"
    lps, avts, semts = orderData(lps, avts, semts)
    ax.errorbar(lps, avts, semts, marker=marker, mfc="w", mec='k',
                capsize=5, alpha=0.7, zorder=1)
    save_name = f"{analysis_dir}/crticial_time_vs_linking_prob.pdf"
    ax.set_ylabel(r"$\tau$")
    ax.set_xlabel("$p_\ell$")

    avtc = avts[avts > 0].mean()
    ax.axhline(avtc, linestyle='--', label=rf"$\tau_c={avtc:3.1f}$", color='k', zorder=2)
    ax.legend()

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()
    print(f"saved to {save_name}")
    return avtc


def getProbDistChainingNumber(data_dir, analysis_dir, show=True):
    """
        Plot the distribution of number of bacteria as a function of linking prob

        Parameters:
            data_dir: str
                directory to look in for trial dirs
            analysis_dir: str
                directory to save data to
            show: bool
                if true show plot
    """
    lps = []
    avNs = []
    semNs = []
    for dd in glob(f"{data_dir}/*"):
        try:
            lp, av_N, sem_N = averageNumberBeforeBuckling(dd)
            lps.append(lp)
            avNs.append(av_N)
            semNs.append(sem_N)
        except Exception as e:
            print(e)
            print(f"Failed for {dd}")

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247 // 2))
    marker = "o"
    lps, avNs, semNs = orderData(lps, avNs, semNs)
    ax.errorbar(lps, avNs, semNs, marker=marker, mfc="w", mec='k',
                capsize=5, alpha=0.7, zorder=1)
    avNc = avNs[avNs > 0].mean()
    ax.axhline(avNc, linestyle='--', label=rf"$N_c={avNc:3.1f}$", color='k', zorder=2)
    save_name = f"{analysis_dir}/number_cells_pre_buckle_vs_linking_prob.pdf"
    ax.set_ylabel("$N$")
    ax.set_xlabel("$p_\ell$")
    ax.legend()

    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()
    print(f"saved to {save_name}")
    return avNc


def plotCellsComparison(lp, buckled_cells, non_buckled_cells, analysis_dir, show=False):
    """
        Compare a buckling chain to a non_buckling chain in the same run
    """
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(2 * 247))

    ax_lim = 0
    for cells in [buckled_cells, non_buckled_cells]:
        try:
            clim = np.max(VisInfectedBiofilm.getAxisLimits(cells)) + 5
        except Exception as e:
            clim = 0
        ax_lim = max(ax_lim, clim)

    for cell in buckled_cells:
        cell.rcm[1] -= 2.5
        cell.addElementToPlot(ax)

    red = (200 / 255, 30 / 255, 30 / 255, 1)
    blue = (138 / 255, 43 / 255, 226 / 255, 1)
    chains = ChainingRodShapedBacterium.getChains(non_buckled_cells)
    for chain, cc in zip(chains, [red, blue]):
        chain_colour = cc
        for cell in chain:
            cell.rcm[1] += 2.5
            cell.colour = cc
            cell.addElementToPlot(ax)

    ax.axis("scaled")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_xlim([-0.5 * ax_lim, 0.5 * ax_lim])
    ax.set_ylim([-20, 20])
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axis('off')

    # Create a nice legend
    legend_labels = []
    legend_elements = []
    handler_map = {}

    cell = ChainingRodShapedBacterium(0, 0, 0, 0, 0, 0, 1, 0, 0)
    legend_elements.append(cell)
    handler_map[ChainingRodShapedBacterium] = cell
    legend_labels.append("buckled chain")

    cell = ChainingRodShapedBacterium(0, 0, 0, 0, 0, 0, 1, 0, 0)
    cell2 = ColouredBacteria(cell, red)
    legend_elements.append(cell2)
    handler_map[ColouredBacteria] = cell2
    legend_labels.append("non buckled chain")

    ax.legend(handles=legend_elements,
              labels=legend_labels,
              handler_map=handler_map,
              loc='upper right',
              framealpha=0,
              frameon=False)

    save_name = f"{analysis_dir}/check_buckling_condition_lp_{lp:3.3f}.pdf"
    fig.savefig(save_name,
                format="pdf", bbox_inches='tight',
                transparent=True)

    if show:
        plt.show()
    plt.close()
    print(f"saved to {save_name}")


def getBucklingSubSample(data_dir, analysis_dir):
    """
        For each trial, visualise a file denoted as buckled and one as not

        Parameters:
            data_dir: str
                directory to look in for trial dirs
            analysis_dir: str
                directory to save data to
    """
    for dd in glob(f"{data_dir}/*"):

        try:
            buckled_file = np.random.choice(glob(f"{dd}/**/buckled_true*", recursive=True))
            print(buckled_file)
            buckled_cells = getCells(buckled_file)
        except Exception as e:
            buckled_cells = []

        try:
            non_buckled_file = np.random.choice(glob(f"{dd}/buckled_false*"))
            print(non_buckled_file)
            non_buckled_cells = getCells(non_buckled_file)
        except Exception as e:
            non_buckled_cells = []
        try:
            hps = getHyperParams(dd)
            plotCellsComparison(hps["LinkingProb"], buckled_cells,
                                non_buckled_cells, analysis_dir)
        except Exception as e:
            print(e)


def getData(data, dir):
    """
    """
    new_data = []
    log = f"{dir}/LcDist.log"
    ldf = pd.read_csv(log, sep="\t")
    ar = ldf["RodAspectRatio"][0]
    new_data.append(ar)

    dat = f"{dir}/critical_dists.dat"
    ddf = pd.read_csv(dat, sep="\t")
    ddf = ddf.to_numpy()
    mns = np.mean(ddf, axis=0)
    # sems = scipy.stats.sem(ddf,axis=0)
    sems = np.std(ddf, axis=0)
    [new_data.append(ll) for rr in zip(mns, sems) for ll in rr]
    data.append(new_data)


def plotDistributions(data_dir, analysis_dir, lc=None, tc=None):
    """
    """
    # Extract data
    data = []
    for dir in glob(f"{data_dir}/*"):
        try:
            getData(data, dir)
        except Exception as e:
            print(e)
    data = np.array(data)

    inds = np.argsort(data[:, 0])
    data = data[inds]
    ar, ls, els, ts, ets, ns, ens = [data[:, ii] for ii in range(data.shape[1])]
    l0 = 0.5 * (ar - 1)
    av_l = 0.5 * (l0 + ar)
    av_t = np.log(3) * av_l
    exp_ls = (1 + 0.5 * (ar + l0)) * ns
    geff = 1 / (1 + 0.5 * ar)
    h_ar = np.linspace(0, ar[-1])

    # fig,ax = plt.subplots(1,1,figsize=setFigSize(247))
    # marker="o"
    # f = lambda x,c: l0*np.exp(c*x)
    # popt,_ = scipy.optimize.curve_fit(f,ts,ls,sigma=els)
    # # print(popt)
    # ax.errorbar(ts,ls,yerr=els,xerr=ets,label="sim")
    #
    # h_ts=np.linspace(0,ts[-1])
    # ax.plot(h_ts,f(h_ts,*popt),"--",
    #         label=fr"$ L =\ell_0\exp {{{popt[0]:3.3f}}}\tau_c$"
    #         )
    # ax.set_ylabel("$\ell_{c}$")
    # ax.set_xlabel(fr"$\tau$")
    # ax.set_xlim([0,None])
    # ax.legend()
    # save_name=f"{analysis_dir}/lcdists_l(t).pdf"
    # if not os.path.exists(analysis_dir):
    #     os.makedirs(analysis_dir)
    # fig.savefig(save_name,format="pdf",bbox_inches='tight',
    #             transparent=True)
    #
    # plt.show()
    # plt.close()

    for key in ["lcs", "taus", "nums"]:
        fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
        marker = "o"

        # Plot and save data
        if key == "lcs":
            continue
            f = lambda x, m, c: m * x + c
            popt, pconv = scipy.optimize.curve_fit(f, ar, ls, sigma=els)
            print(popt)
            errors = np.sqrt(np.diag(pconv))
            print(popt)
            print(errors)
            ax.errorbar(ar, ls, els, marker=marker, mfc="w", mec='k',
                        capsize=5, alpha=0.7, label="sim")
            # ax.plot(ar,exp_ls,":",
            #         label=r"$(1+\langle\ell\rangle)N_c$")
            ax.plot(h_ar, f(h_ar, *popt), "--",
                    label=fr"$L=m\alpha+c$"
                    )
            # ax.plot(ar,,"--")
            ax.set_ylabel("$L$")
            annotation = (f"$m={popt[0]:3.3f}\pm{errors[0]:3.3f}$\n" +
                          f"$c={popt[1]:3.3f}\pm{errors[1]:3.3f}$\n" +
                          r"$\tilde{\ell}_c=" + f"{lc:3.3f}$"
                          )
            ax.text(0.05, 0.95, annotation,
                    horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes)
            if lc:
                ax.axhline(lc, linestyle="--", color="k", label=r"$L$")
        elif key == "taus":
            ts /= 4
            ets /= np.sqrt(4)
            f = lambda x, m, c: m * x + c
            popt, pconv = scipy.optimize.curve_fit(f, ar, ts, sigma=ets)
            errors = np.sqrt(np.diag(pconv))
            print(popt)
            print(errors)

            f = lambda x, a: a + 0.25 * (x+1) * np.log((x+1))
            popt, pconv = scipy.optimize.curve_fit(f, ar, ts, sigma=ets)
            errors = np.sqrt(np.diag(pconv))
            print("power law", popt)
            print("errors", errors)
            # ax.errorbar(ar,ts,ets,marker=marker,mfc="w",mec='k',
            #             capsize=5,alpha=0.7,label="sim")
            ax.plot(ar, ts, marker=marker, mfc="w", mec='k', markevery=0.2,
                    alpha=0.7, label="discrete")
            ax.fill_between(ar, ts + ets, ts - ets, color='gray', alpha=0.2)
            ax.plot(ar, f(ar, *popt), "-.", alpha=0.7, color='r',
                    label=fr"fit"
                    )
            ax.set_ylabel(r"$t_{c}$ (\si{\hour})")
            # ax.text(0.05,0.95,f"$m={popt[0]:3.3f}\pm{errors[0]:3.3f}$\n$c={popt[1]:3.3f}\pm{errors[1]:3.3f}$\n"+r"$\tilde{\tau}_c="+f"{tc:3.3f}$",
            #         horizontalalignment='left',
            #         verticalalignment='top', transform=ax.transAxes)
            if tc:
                # ax.axhline(tc, linestyle="--", color="k", label=r"continuum pred.")
                ax.plot(ar, tc * (ar + 1) / (5 + 1), linestyle="--", color="k", label=r"continuum")
        elif key == "nums":
            continue
            f = lambda x, n0, a: n0 * x ** a
            popt, _ = scipy.optimize.curve_fit(f, ar, ns, sigma=ens)
            print(popt)
            ax.errorbar(ar, ns, ens, marker=marker, mfc="w", mec='k',
                        capsize=5, alpha=0.7, label="sim")
            ax.plot(h_ar[1:], f(h_ar[1:], *popt), "--",
                    label=fr"$N_c={popt[0]:3.3f}\alpha^{{{popt[1]:3.3f}}}$"
                    )
            # ax.plot(ar,1+2*np.floor(ts/av_t),":")
            ax.set_ylabel("$N_{c}$")

        ax.set_xlabel(r"$\ell_{div}$ (\si{\um})")
        # ax.set_xlim([0,None])
        ax.legend(loc="best")
        save_name = f"{analysis_dir}/lcdists_{key}_thesis.pdf"
        if not os.path.exists(analysis_dir):
            os.makedirs(analysis_dir)
        print(save_name)
        fig.savefig(save_name, format="pdf", bbox_inches='tight',
                    transparent=True)

        plt.show()
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='findCriticalLength.py',
        usage='%(prog)s [options] path',
        description='Find probability distribution of buckling lengths',
        epilog='Author: Rory Claydon'
    )
    parser.add_argument("-BD", "--base-dir", type=str,
                        required=False, default="../GeneratedOutput",
                        action="store", dest="base_path",
                        help="default directory to look for inputs")
    parser.add_argument("-DD", "--data-dir", type=str, required=False,
                        default="SimOutput",
                        action="store", dest="data_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("-AD", "--analysis-dir", type=str, required=False,
                        default="AnalysisResults",
                        action="store", dest="analysis_dir",
                        help="default directory name in base for inputs")
    parser.add_argument("-RD", "--run-dir", type=str, required=False,
                        default="ChainingCriticalLengths_Aspect_5_Thesis",
                        # default="ChainingCriticalLengths_Aspect_5",
                        action="store", dest="run_dir",
                        help="default directory name in data dir for inputs")
    parser.add_argument("--distribution-dir", type=str, required=False,
                        # default="LcDist",
                        default="LcDistFinal",
                        action="store", dest="dist_dir",
                        help="top directory for the distribution sweeps")

    args = parser.parse_args()

    # Set up distribution directories
    dist_dir = f"{args.base_path}/{args.data_dir}/{args.dist_dir}/"
    dist_analysis_dir = f"{args.base_path}/{args.analysis_dir}/{args.dist_dir}/"

    # getBucklingSubSample(dist_dir,dist_analysis_dir)

    dir = f"{args.base_path}/{args.data_dir}/{args.run_dir}/"
    analysis_dir = f"{args.base_path}/{args.analysis_dir}/{args.run_dir}/"
    tc = getProbDistChainingTime(dir, analysis_dir, show=True)
    lc = getProbDistChainingLength(dir, analysis_dir, show=True)
    N = getProbDistChainingNumber(dir, analysis_dir, show=True)
    print(f"{lc=} {N=}")
    getProbDistBuckling(dir, analysis_dir, lc=lc, avN=N, show=True)

    # getBucklingSubSample(dir,analysis_dir)

    tc = 3.178
    lc = 4 * tc
    # plotDistributions(dist_dir,dist_analysis_dir,lc=lc,tc=tc)
