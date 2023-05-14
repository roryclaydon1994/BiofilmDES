"""
Frozen mode analysis of compressible, growing filament
"""
import os

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from joblib import Parallel, delayed
# import multiprocessing as mp
from timeit import default_timer as timer

# Custom modules
import findBucklingEigenvalues as fbe
import utilities as ut

ut.setMPL()


class HPs:
    d = 1  # diameter (micron)
    g = 0.75 / np.log(2)  # growth rate (per hr)
    Y0 = 4e6  # Young's modulus (Pa)
    l0 = 2  # initial length (microns)
    xi = 200  # friction coefficient (Pa hr)
    k0 = (Y0 * d ** 2) / (g * xi * l0 ** 2)  # compression modulus DES value
    beta = (Y0 * d ** 2) / (g * xi * l0 ** 4)  # bending modulus DES value


def getTc(kappa, a, N, t_max):
    print(f"\n--------------------  {kappa=:2.2e}_ --------------------")
    if np.isclose(kappa, 0.0):
        t_c, t_err = 0, 0
    else:
        t_c, t_err = fbe.findTCrit(kappa=kappa, beta=HPs.beta, a=a, N=N, t_max=t_max, bc='Fily')
    return t_c, t_err


def runSweep(save_path='kappa_sweep_critical_ts.txt'):
    # t_max=np.log(k)
    t_max = 10
    N = 501
    a, da = np.linspace(0, 1, N, retstep=True, endpoint=True)
    # _chis = np.linspace(0, HPs.beta / HPs.k0, num=51, endpoint=True)
    # _chis = np.logspace(0, 0, num=7, endpoint=True) * HPs.beta / HPs.k0
    kappas = HPs.k0 * np.logspace(0, -5, num=51, endpoint=True)
    kappas = kappas[::-1]

    def _getTc(c):
        return getTc(c, a, N, t_max)

    data = Parallel(n_jobs=1, verbose=10)(delayed(_getTc)(k) for k in kappas)
    data = np.array(data)
    t_cs, t_errs = data[:, 0], data[:, 1]
    _chis = HPs.beta / kappas
    np.savetxt(save_path, np.transpose([_chis, t_cs, t_errs]))
    return _chis, t_cs, t_errs


def plotSweepResults(chis, ts, errs, fname="kappa_sweep.pdf", savedir='./', show=True):
    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(1))
    ax.set_xscale('log')
    chis = 1.0 / chis
    p1 = ax.plot(chis, ts, marker='o', mfc="w", mec='k', alpha=0.7,
                 label=rf"$\tau_c(\chi)$")
    ax.fill_between(chis, ts - errs, ts + errs, color='gray', alpha=0.2)

    # Small kappa regime
    f = lambda x, a, b: a + b * np.log(x)
    popt, pcov = scipy.optimize.curve_fit(f, chis[:18], ts[:18])
    print(popt, np.sqrt(np.diag(pcov)), f"{0.5/HPs.g=}")
    ax.plot(chis[:18], f(chis[:18], *popt), alpha=0.8, color='k', ls='--')

    # transition regime
    # f = lambda x, a, b: a + b * np.log(x)
    # popt, pcov = scipy.optimize.curve_fit(f, chis[20:33], ts[20:33])
    # print(popt, np.sqrt(np.diag(pcov)))
    # ax.plot(chis[20:33], f(chis[20:33], *popt), alpha=0.8, color='k', ls='-.')

    ax.set_xlabel(r"$\kappa/\beta$")
    ax.set_ylabel(r"$\tau_c$")
    # previous tc 3.503881292257138
    t_inextensible = ts[-8].mean()
    ax.axhline(t_inextensible, xmin=0.6, linestyle='--', alpha=0.8, color='k')
    ax.axvline(np.sqrt(1/HPs.beta), alpha=0.6, linestyle=':', color='m')
    # ax.axhline(0.25*np.log(HPs.beta))
    # ax.plot(chis, 0.5*np.log(chis*HPs.beta))

    fig.savefig(f"{savedir}/{fname}", format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def getData(save_path, force_recalculate=False):
    if force_recalculate:
        chis, tcrits, terrs = runSweep(save_path)
    else:
        try:
            data = np.loadtxt(save_path)
            chis = data[:, 0]
            tcrits = data[:, 1]
            terrs = data[:, 2]
        except FileNotFoundError as e:
            print(e)
            chis, tcrits, terrs = runSweep(save_path)
        except Exception as e:
            print(e)
            quit()

    return chis, tcrits, terrs


def checkLambdaAnalytic():
    m = 101
    a = np.linspace(0, 1, m)
    kappas = HPs.k0 * np.logspace(-4, -3, num=1, endpoint=True)
    Ms = [51, 11]
    for k, M in zip(kappas, Ms):
        ll = fbe.getLam0(a=a, k=k, N=m, t_max=4)
        fbe.plotLambda(ll=ll, a=a, kappa=k, t_max=4, M=M)


if __name__ == "__main__":
    # checkLambdaAnalytic()
    # quit()
    data_filename = 'kappa_sweep_critical_ts.txt'
    result_directory = f"../GeneratedOutput/frozen_mode_extensible/"
    os.makedirs(result_directory, exist_ok=True)
    chis, ts, es = getData(f"{result_directory}/{data_filename}".replace('//', '/'), force_recalculate=False)
    plotSweepResults(chis, ts, es, savedir=result_directory)
