"""
    Old parameters
    d=1 # diameter (micron)
    g=6 # growth rate (per hr)
    Y0=4e6 # Young's modulus (Pa)
    Y=2*Y0
    l0=1 # initial length (microns)
    xi=200 # friction coefficient (Pa hr)
"""
# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.stats
from scipy import sparse
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from pprint import pprint
from scipy.special import erf
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, curve_fit
import matplotlib as mpl
import os
from tqdm import tqdm
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                   mark_inset)
# from joblib import Parallel, delayed
# import multiprocessing as mp
from timeit import default_timer as timer
# import gc

# Custom modules
from visualise_biofilm import setFigSize
from findCriticalLengthProb import plotDistributions
from findBucklingEigenvalues import makeD1, makeD2, makeD3, makeD4
import utilities as ut

ut.setMPL()


def getLambda(s, t):
    return 0.5 * np.exp(2 * t) * s * (s - 1)


def plotLambda(ll, s):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.plot(s, ll(s, 0), label=r"$t={0}$")
    ax.plot(s, ll(s, eta_max / 2), label=r"$t=\frac{1}{4}\ln(\beta)$")
    ax.plot(s, ll(s, eta_max), label=r"$t=\frac{1}{2}\ln(\beta)$")
    plt.legend()
    save_name = f"{savedir}/check_lambda.pdf"
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    plt.close()


def setUpMatrix(s, eta, beta=1):
    N = len(s)
    et_conv = 1 / beta

    a4 = -np.exp(-4 * eta / beta)
    A4 = diags([a4], offsets=[0], format='csr', shape=(N, N))

    a2 = 0.5 * s * (s - 1) * et_conv
    A2 = diags([a2], offsets=[0], format='csr', shape=(N, N))

    a1 = 2 * s * et_conv
    A1 = diags([a1], offsets=[0], format='csr', shape=(N, N))

    L = A4.dot(D4) + A2.dot(D2) + A1.dot(D1)
    return L


def solveEigenProblemAtT(s, eta, beta):
    L = setUpMatrix(s, eta, beta=beta)
    # cond=np.linalg.cond(L.toarray())
    # print(f"{cond}")
    # quit()
    # try:
    #     w,v=eigs(L,k=2,which="LR")
    # except Exception as e:
    #     w=0
    #     print(e)
    #     print(f"Failed at {t}")
    #     quit()

    w, v = np.linalg.eig(L.toarray())
    w = scipy.linalg.eigvals(L.toarray())
    inds = np.argpartition(np.real(w), -2)[-2:]
    inds = inds[np.argsort(np.real(w[inds]))][::-1]  # not necessarily ordered
    w = w[inds]
    v = v[inds]
    return w


def solveEigenProblem(s, eta_max, beta):
    min_ws = []
    etas = np.linspace(0, eta_max, 30)

    for eta in etas:
        w = solveEigenProblemAtT(s, eta=eta, beta=beta)
        min_ws.append(w)
    return min_ws, etas


def findTCrit(s, eta_max, beta):
    # f=lambda tc: np.imag(solveEigenProblemAtT(s,tc,beta=beta)[0])+1e-7
    # sol=root_scalar(f,bracket=[0,eta_max],x0=0.25*eta_max)
    # print(sol)

    f = lambda etac: np.real(solveEigenProblemAtT(s, etac, beta=beta)[1])
    sol = root_scalar(f, bracket=[0, eta_max], x0=0.25 * eta_max)
    print(sol)

    return sol.root, 0.0


def sweepFn(s, num=25, dir=""):
    betas_ = np.logspace(-4, 4, num, endpoint=True)

    eta_crits = []
    eta_errs = []
    betas = []  # successful betas

    for bb in tqdm(betas_):
        try:
            eta_max = max(4 * bb, np.log(bb) * bb)
            eta_crit, eta_err = findTCrit(s, eta_max, bb)
            plotAllEigs(eta_crit, beta=bb, dir=dir, eta_max=eta_max)
            betas.append(bb)
            eta_crits.append(eta_crit)
            eta_errs.append(eta_err)
            print(f"{bb=}")
        except Exception as e:
            print(e)

    return np.array(betas), np.array(eta_crits), np.array(eta_errs)


def plotEigenvalues(ts, mws, beta=1, eta_crit=None, save_name="bucklingEigenvaluesPython.pdf", show=True):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.plot(ts, np.real(mws[:, 0]), "b-", marker='o', mfc='w', label="$\Re(\omega_1)$",
            markevery=(0.0, 0.15), alpha=0.7)
    ax.plot(ts, np.imag(mws[:, 0]), "--", marker='d', mfc='w', label="$\Im(\omega_1)$",
            markevery=(0.05, 0.15), alpha=0.7)
    ax.plot(ts, np.real(mws[:, 1]), "m-", marker='o', label="$\Re(\omega_2)$",
            markevery=(0.0, 0.17), alpha=0.7)
    ax.plot(ts, np.imag(mws[:, 1]), ":", marker='d', label="$\Im(\omega_2)$",
            markevery=(0.05, 0.17), alpha=0.7)
    ax.set_xlim([ts[0], ts[-1]])
    ax.set_ylim([-10, 10])

    if eta_crit:
        # ax.axvline(eta_crit,linestyle=':',color='k',label=r"$\tau_c$")
        ax.axvline(eta_crit, linestyle=':', color='k', ymax=0.5)

        txt = ax.text(0.63, 0.25, rf"$\tau_c={eta_crit:3.3f}$",
                      horizontalalignment='center',
                      verticalalignment='center', transform=ax.transAxes)
    plt.legend(loc="upper left", borderpad=0.2, borderaxespad=0.2, handletextpad=0.5,
               ncol=2, columnspacing=1)

    ax.set_xlabel(r"$\tau_0$")
    ax.set_ylabel(r"$\omega$")

    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def plotAllEigs(etacrit, beta, dir, eta_max):
    minw, ts = solveEigenProblem(s, eta_max, beta)
    minw = np.array(minw)
    plotEigenvalues(ts, minw, etacrit,
                    save_name=f"{dir}/bucklingEigenvaluesPython_{bc}_{N=}_tc_{beta=:6.6f}.pdf",
                    show=False
                    )
    plotLogEigenvalues(ts, minw, etacrit,
                       save_name=f"{dir}/bucklingLogEigenvaluesPython_{bc}_{N=}_tc_{beta=:6.6f}.pdf",
                       show=False
                       )


def plotSweepResults(betas, etacrits, etaerrs=None, save_name="tcs.pdf", show=True):
    tcrits = np.copy(etacrits) / betas
    terrs = np.copy(etaerrs) / betas
    # tcrits /= g
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    f = lambda x, a, b: a + b * np.log(x)
    # betas *= 1e-12 * g * xi * l0 ** 4
    popt, pcov = curve_fit(f, betas, tcrits)
    print(popt, np.sqrt(np.diag(pcov)))

    af = f"{popt[0]:2.2f}"
    bf = f"{popt[1]:1.2e}"
    print(f"{bf=}")
    p1 = ax.plot(betas, tcrits, marker='o', mfc="w", mec='k', alpha=0.7,
                 label=rf"Simulation")
    p2 = ax.plot(betas, f(betas, *popt), "--", zorder=3, alpha=0.7, color='k',
                 label=rf"$t_c = p + q\ln(B)$"
                 # label=(f"Fit ${{{af}}}"+r"\log{\qty(\num{"+
                 #         bf+
                 #         "}" + r"B)}$"
                 #         )#+f"\n$a={af}$,\n$b={bf}$"
                 )
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\tau_c$")
    plt.legend()

    save_name = f"{save_name}"
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()
    print(f"Expected buckling should be {f(beta,*popt)=}")
    quit()

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.plot(betas, np.exp(tcrits))

    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$\ell_c/\ell_0$")

    save_name = save_name.replace("tcs_", "ells_")
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def plotLogEigenvalues(ts, mws, eta_crit=None, save_name="bucklingEigenvaluesPythonLog.pdf", show=True):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.semilogy(ts, np.abs(np.real(mws[:, 0])), "b-", label="$|\Re(\omega_1)|$")
    ax.semilogy(ts, np.abs(np.imag(mws[:, 0])), "--", label="$|\Im(\omega_1)|$")
    ax.semilogy(ts, np.abs(np.real(mws[:, 1])), "m-", label="$|\Re(\omega_2)|$")
    ax.semilogy(ts, np.abs(np.imag(mws[:, 1])), ":", label="$|\Im(\omega_2)|$")

    if eta_crit:
        ax.axvline(eta_crit, linestyle=':', color='k', label=r"$\eta_c$")
    plt.legend(loc="best")

    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("$\omega$")

    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


if __name__ == "__main__":

    N = 101
    s, ds = np.linspace(0, 1, N, retstep=True, endpoint=True)

    plot_single = True

    savedir = "../GeneratedOutput/frozen_mode/inextensible_correct_parameters2/"
    if not os.path.exists(f"{savedir}"):
        os.makedirs(savedir)

    bc = "Fily"
    D1 = makeD1(ds, N)
    D2 = makeD2(ds, N)
    D3 = makeD3(ds, N, bc=bc)
    D4 = makeD4(ds, N, bc=bc)
    # plotLambda(getLambda,s)

    # d=1 # diameter (micron)
    # g=6 # growth rate (per hr)
    # Y0=4e6 # Young's modulus (Pa)
    # Y=2*Y0
    # l0=1 # initial length (microns)
    # xi=200 # friction coefficient (Pa hr)
    d = 1  # diameter (micron)
    g = 0.75 / np.log(2)  # growth rate (per hr)
    print(f"{g=}")
    Y0 = 4e6  # Young's modulus (Pa)
    l0 = 2  # initial length (microns)
    xi = 200  # friction coefficient (Pa hr)
    B = Y0 * d ** 4  # bending modulus in Pa micron^4
    beta = B / (g * xi * l0 ** 4)  # Bending modulus
    eta_max = beta * np.log(beta)  # This is time, not sure why I used this symbol again

    if plot_single:

        print(f"{beta=} {eta_max=}")

        minw, etas = solveEigenProblem(s, eta_max, beta)
        minw = np.array(minw)
        print(f"The maximum real parts are: {max(np.real(minw[:, 0]))} and {max(np.real(minw[:, 1]))}")
        eta_crit, eta_err = findTCrit(s, eta_max, beta)

        plotEigenvalues(etas / beta, minw, beta=beta, eta_crit=eta_crit / beta,
                        save_name=f"{savedir}/bucklingEigenvaluesPython_{N=}_{bc}_etac_{beta=:6.6f}.pdf"
                        )
        plotLogEigenvalues(etas, minw, eta_crit=eta_crit,
                           save_name=f"{savedir}/bucklingLogEigenvaluesPython_{N=}_{bc}_tc_{beta=:6.6f}.pdf"
                           )
        # use_data_dir='LcDistHighRes'
        use_data_dir = 'LcDistFinal/'
        dist_dir = f"../GeneratedOutput/SimOutput/{use_data_dir}/"
        dist_analysis_dir = f"{savedir}/{use_data_dir}/"
        tcrit = eta_crit / beta
        tcrit /= g
        print(f"{tcrit=}")
        plotDistributions(dist_dir, dist_analysis_dir, tc=tcrit, lc=l0 * np.exp(tcrit))
    else:
        savedir += "/sweep/"
        if not os.path.exists(f"{savedir}"):
            os.makedirs(savedir)
        tc_data_save_name = f"etas_data_{N=}_{bc}.txt"
        fpath = f"{savedir}/{tc_data_save_name}"
        if not os.path.exists(fpath):
            betas, etacrits, etaerrs = sweepFn(s, dir=savedir)
            np.savetxt(fpath, np.transpose([betas, etacrits, etaerrs]))
        else:
            data = np.loadtxt(fpath)
            betas = data[:, 0]
            etacrits = data[:, 1]
            etaerrs = data[:, 2]

        tc_fig_save_name = f"{savedir}/tcs_{N=}_{bc}.pdf"
        plotSweepResults(betas, etacrits, etaerrs, save_name=tc_fig_save_name,
                         show=True)
