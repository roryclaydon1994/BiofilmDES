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
import utilities as ut
from visualise_biofilm import setFigSize
from ChainingRodShapedBacteria import ChainingRodShapedBacterium
from FindDefects import intialiseElementsFromData
from visualise_infected_biofilm import VisInfectedBiofilm

ut.setMPL()

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))


def getScalingLabel(iv, dv):
    scaling_lookup = {
        'nums': {
            'RodGrowthRate': "{-1/4}",
            'RodAspectRatio': "{-1/2}",
            'BendRig': "{1/4}",
            "Kappa": "0"
        },
        'lengths': {
            'RodGrowthRate': "{-1/4}",
            'RodAspectRatio': "{1/2}",
            'BendRig': "{1/4}",
            "Kappa": "0"
        },
        'taus': {
            'RodGrowthRate': "{-5/4}",
            'RodAspectRatio': "{1/2}",
            'BendRig': "{1/4}",
            "Kappa": "0"
        }
    }
    return scaling_lookup[dv][iv]


def getAidanNc(l, mu, B, xi=200):
    """
    Aidan scaling argument
    """
    # a = B/(mu*xi*l**2)
    a = B / (mu * xi * l ** 3)
    return a ** (0.25)


def getAidanLc(l, mu, B, xi=200):
    """
    Aidan scaling argument
    """
    return getAidanNc(l, mu, B, xi=200) * l


def getAidanTc(l, mu, B, xi=200):
    """
    Aidan scaling argument
    """
    return getAidanNc(l, mu, B, xi=200) * l / mu
    # return getAidanNc(l,mu,B,xi=200)*l/mu


def getTheoreticalNc(l, mu, B, xi=200, pole=1):
    """
    Takes dimensional units
    """
    if pole == 1:
        f_star = (B / (1 + l)) * (2 * l - 1 + 1 / l)  # Considering one pole
    elif pole == 2:
        f_star = (B / (1 + 0.5 * l)) * (l - 1 + 2 / l)
    else:
        f_star = (B / (1 + l)) * (2 * l + 2 + 1 / l)  # Considering one pole
        print(f"{pole=} for full bending force")
    return np.sqrt(f_star / (mu * l * xi))


def getTheoreticalLc(l, mu, B, pole=1):
    return getTheoreticalNc(l, mu, B, pole=1) * l


def getTheoreticalTc(l, mu, B, pole=1):
    return getTheoreticalNc(l, mu, B, pole=1) * l


def getFitFn(ivargs, dvargs, ls, mus, Bs):
    """
        Return the function with the parameter which is varying as input
    """
    if ivargs['ivar'] == 'RodAspectRatio':
        fn = lambda x: dvargs['g'](l=x, mu=mus, B=Bs)
    elif ivargs['ivar'] == 'Kappa':
        fn = lambda x: dvargs['g'](l=ls, mu=mus, B=Bs)
    elif ivargs['ivar'] == 'BendRig':
        fn = lambda x: dvargs['g'](l=ls, mu=mus, B=x)
    elif ivargs['ivar'] == 'RodGrowthRate':
        fn = lambda x: dvargs['g'](l=ls, mu=x, B=Bs)
    else:
        print(":(")
        quit()
    return fn


def getBestFit(x, y, fn, first_only=False):
    fns = [lambda x, a: a * fn(x),
           lambda x, a: fn(x) + a
           ]
    res_min = np.inf

    for ii, _fn in enumerate(fns):
        popt, pcov = scipy.optimize.curve_fit(_fn, x, y)
        if np.max(np.sqrt(np.diag(pcov))) == np.inf:
            continue
        data = _fn(x, *popt)
        res = np.mean((y - data) ** 2)
        if res <= res_min:
            best_popt = popt
            best_pcov = pcov
            best_fn = _fn
            best_ii = ii
        if first_only:
            break

    print(f"Best fn at index {best_ii}")
    return best_fn, popt, pcov


def plotParameter(_dir, ivargs, dvargs, save_dir, show=False):
    if dvargs['dvar'] == "taus":
        return 0, 0, 0, 0

    ivars = []
    dvars = []
    errors = []
    mus = []
    Bs = []
    ls = []
    for dd in os.listdir(_dir):
        if dd == "L_2" or dd == "L_2.01508" or dd == "L_2.03015":
            # print("here")
            continue
        log_file = f"{_dir}/{dd}/LcDist.log"
        dat_file = f"{_dir}/{dd}/critical_dists.dat"
        log_df = pd.read_csv(log_file, sep='\t')
        log_df["BendRig"] *= 4e6  # Convert to dimensional units
        log_df["Kappa"] *= 4e6  # Convert to dimensional units
        log_df["RodAspectRatio"] = 1 + 0.75 * log_df["RodAspectRatio"] - 0.25
        dat_df = pd.read_csv(dat_file, sep='\t')
        ivars.append(log_df.at[0, ivargs['ivar']])
        mus.append(log_df.at[0, "RodGrowthRate"])
        Bs.append(log_df.at[0, "BendRig"])
        ls.append(log_df.at[0, "RodAspectRatio"])
        dvars.append(dat_df[dvargs['dvar']].mean())
        errors.append(scipy.stats.sem(dat_df[dvargs['dvar']]))
        if dvargs['dvar'] == 'taus':
            # print(f"scaling tau with mu: {mus[-1]}")
            dvars[-1] /= mus[-1]
            if _dir == grwth_dir:
                dvars[-1] /= 2e-4
                # print(f"scaled time to {dvars[-1]}")
            errors[-1] /= np.sqrt(mus[-1])

    idx = np.argsort(ivars)
    ivars = np.array(ivars)[idx]
    dvars = np.array(dvars)[idx]
    Bs = np.array(Bs)[idx]
    ls = np.array(ls)[idx]
    mus = np.array(mus)[idx]
    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(1))
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    marker = "o"
    # ax.errorbar(ivars,dvars,errors,ls='-',mfc='None',mec='k',capsize=5,
    #             marker=marker,
    #             markevery=(0.0,0.2),
    #             errorevery=int(0.1*len(dvars)),
    #             zorder=1,
    #             label="sim")
    ax.plot(ivars, dvars, ls='-', mfc='None', mec='k',
            marker=marker,
            markevery=(0.0, 0.2),
            zorder=1,
            label="sim")
    ax.fill_between(ivars, dvars + errors, dvars - errors, color="gray", alpha=0.2)
    # style=['--',':','-.']
    # for ii,pole in enumerate(range(1,4)):
    #     ax.plot(ivars,dvargs['f'](l=ls,mu=mus,B=Bs,pole=pole),style[ii],label=f'{pole=}')

    fn = getFitFn(ivargs, dvargs, ls, mus, Bs)
    fit_fn, popt, pcov = getBestFit(ivars, dvars, fn, first_only=True)

    print(dvargs['dvar'], ivargs['ivar'], popt, np.sqrt(np.diag(pcov)))

    # lbl=rf"$dv \sim iv^{{{getScalingLabel(ivargs['ivar'],dvargs['dvar'])}}}$"
    lbl = rf"theory"
    lbl = lbl.replace('dv', dvargs['symbol'])
    lbl = lbl.replace('iv', ivargs['symbol'])
    sim_dat = fit_fn(ivars, *popt)
    ax.plot(ivars, sim_dat, "g:",
            marker='x',
            mec='k',
            mfc='None',
            markevery=(0.2 / 3, 0.2),
            label=lbl,
            alpha=0.8,
            zorder=2
            )

    # lbl=rf"$dv^\mathrm{{sim}} \sim iv^{{{b:3.3f}}}$"
    # lbl=rf"$dv^\mathrm{{sim}} \sim iv^{{{b:2.2f}}}$"
    # lbl=rf"$dv^\mathrm{{sim}}$"
    lbl = rf"fit"
    lbl = lbl.replace('dv', dvargs['symbol'])
    lbl = lbl.replace('iv', ivargs['symbol'])

    # if ivargs['ivar'] == "RodAspectRatio":
    #     power_law = lambda x,a,b,c: a*x**b+c
    #     popt,pcov=scipy.optimize.curve_fit(power_law,ivars,dvars)
    #     a,b,c=popt
    #     th_dat = power_law(ivars,*popt)
    #     bl = th_dat > 0
    #     ax.plot(ivars[bl],th_dat[bl],"m:",
    #             label=lbl,
    #             alpha=0.8,
    #             zorder=3
    #             )
    #     print(f"{dvargs['dvar']} {ivargs['ivar']} a: {a} b: {b} c: {c}")
    # else:
    power_law = lambda x, a, b: a * x ** b
    popt, pcov = scipy.optimize.curve_fit(power_law, ivars, dvars)
    a, b = popt
    aerr, berr = np.sqrt(np.diag(pcov))
    th_dat = power_law(ivars, *popt)
    bl = th_dat >= 0
    if ivargs['ivar'] == 'Kappa':
        bl = th_dat > 0
    ax.plot(ivars[bl], th_dat[bl], "k--",
            markevery=(2 * 0.2 / 3, 0.2),
            marker='d',
            mec='k',
            mfc='None',
            label=lbl,
            alpha=1.0,
            zorder=3
            )
    # try:
    #     print(f"{dvargs['dvar']} {ivargs['ivar']} a: {a} b: {b} c: {c}")
    # except Exception as e:
    #     print(f"{dvargs['dvar']} {ivargs['ivar']} a: {a} b: {b}")

    ax.set_xlabel(ivargs['label'])
    ax.set_ylabel(dvargs['label'])
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
    ax.yaxis.set_minor_locator(ticker.FixedLocator([-1]))
    if ivargs['ivar'] == 'RodAspectRatio' and dvargs['dvar'] == 'nums':
        plt.legend()
    ax.set_ylim([0.925 * dvars.min(), 1.075 * dvars.max()])
    # print(ax.get_yticklabels())
    fig.savefig(
        f"{save_dir}/{dvargs['dvar']}_{ivargs['ivar']}_new_scaling.pdf",
        bbox_inches="tight",
        transparent=True
    )
    # if show:
    #     plt.show(block=True)
    # plt.close()

    return a, b, aerr, berr


if __name__ == "__main__":
    kappa_dir = "../GeneratedOutput/SimOutput/ScalingExpKappa"
    grwth_dir = "../GeneratedOutput/SimOutput/ScalingExpGrwthRateCorrected"
    bendy_dir = "../GeneratedOutput/SimOutput/ScalingExpBendRig"
    lengt_dir = "../GeneratedOutput/SimOutput/LcDistFinal"
    save_dir = "../GeneratedOutput/AnalysisResults/ScalingExp"

    kappa_var = {"ivar": 'Kappa', "symbol": r"\kappa",
                 "label": r'$K$ (\si{\newton})', 'dir': kappa_dir}
    grwth_var = {"ivar": 'RodGrowthRate',
                 "symbol": r'\mu_0',
                 "label": '$\mu_0$ (\si{\micro\metre\per\hour})', 'dir': grwth_dir}
    bendy_var = {"ivar": 'BendRig', "symbol": r'B',
                 "label": r'$B$ (\si{\newton\metre^2})', 'dir': bendy_dir}
    lengt_var = {"ivar": 'RodAspectRatio',
                 "symbol": r'\langle \ell \rangle',
                 'label': r'$\langle \ell \rangle + d_0$ (\si{\micro\metre})', 'dir': lengt_dir}
    vargs = [kappa_var, grwth_var, bendy_var, lengt_var]

    num_var = {"dvar": 'nums', 'symbol': "N_c", 'label': "$N_c$",
               "f": getTheoreticalNc,
               "g": getAidanNc
               }
    len_var = {"dvar": 'lengths', 'symbol': 'L_c',
               'label': "$L_c$ (\si{\micro\metre})",
               "f": getTheoreticalLc,
               "g": getAidanLc
               }
    tau_var = {"dvar": 'taus', 'symbol': r'\tau_c',
               'label': r"$T_c$ (\si{\hour})",
               "f": getTheoreticalTc,
               "g": getAidanTc
               }
    dargs = [num_var, len_var, tau_var]

    scalings = []
    for dv in dargs:
        for iv in vargs:
            a, b, aerr, berr = plotParameter(iv['dir'], iv, dv, save_dir=save_dir, show=True)
            scalings.append([dv['dvar'], iv['ivar'], a, aerr, b, berr])
    plt.show()
    with open(f"{save_dir}/scaling_power_laws.csv", "w") as f:
        f.write("iv\tdv\ta\taerr\tb\tberr\n")
        for scale in scalings:
            f.write(f"{scale[0]}\t{scale[1]}\t{scale[2]:2.3f}\t{scale[3]:2.3f}\t{scale[4]:2.3f}\t{scale[5]:2.3f}\n")
