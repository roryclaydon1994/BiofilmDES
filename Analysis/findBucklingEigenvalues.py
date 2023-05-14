r"""
    Determine the fastest growing eigenmode of the linearised 2d equation:
    K\lambda_0'\vartheta_1' + K\qty(\lambda_0-1)\vartheta_1''-
    B\qty(\frac{\vartheta_1''}{\lambda_0})''+
    K\lambda_0\vartheta_1'+
    \xi\qty(\omega\lambda_0\vartheta_1+g\alpha\vartheta_1')=0

    \vartheta'(\alpha=\{-1,1\},t)=0 \qq{and} \vartheta''(\alpha=\{-1,1\},t)=0\,,

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
import utilities as ut

ut.setMPL()

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


def pt(t, k=1, L=1):
    bn = lambda n: k * (np.pi * n / L) ** 2
    p_0 = lambda n: np.exp(-t) / (n * np.pi)
    p_1 = lambda n: np.nan_to_num(
        np.sqrt(8 * np.pi * bn(n)) * np.exp(np.exp(-2 * t) * bn(n) / (2))
    )
    erfs = lambda n: erf(np.sqrt(bn(n) / (2))) - erf(np.exp(-t) * np.sqrt(bn(n) / (2)))
    p_2 = lambda n: 4 * np.exp((np.exp(-2 * t) - 1) * bn(n) / (2))
    p = lambda m: p_0(m) * (p_1(m) * erfs(m) + p_2(m) - 4 * np.exp(t))
    return p


def lam0(alpha, t, k=1, L=1, M=11):
    p = pt(t, k=k, L=L)
    Ms = np.arange(1, M, 2)
    ll = p(Ms[:, np.newaxis]) * np.sin(Ms[:, np.newaxis] * alpha[np.newaxis, :] * np.pi / L)
    return 1 + np.sum(ll, axis=0)


def dlam0_da(alpha, t, g=0.4, beta=1, k=1, L=1, M=8):
    p = pt(t, g=g, beta=beta, k=k, L=L)
    Ms = np.arange(1, M, 2)
    aa = Ms[:, np.newaxis] * np.pi / L
    ll = p(Ms[:, np.newaxis]) * aa * np.cos(aa * alpha[np.newaxis, :])
    return np.sum(ll, axis=0)


def dlam02_da2(alpha, t, g=0.4, beta=1, k=1, L=1, M=8):
    p = pt(t, g=g, beta=beta, k=k, L=L)
    Ms = np.arange(1, M, 2)
    aa = Ms[:, np.newaxis] * np.pi / L
    ll = -p(Ms[:, np.newaxis]) * (aa ** 2) * np.sin(aa * alpha[np.newaxis, :])
    return np.sum(ll, axis=0)


def checkDefinitions(k, ll, N, t_max):
    alpha, da = np.linspace(0, 1, N, retstep=True, endpoint=True)
    D2 = makeD2(dr=da, N=N)
    t = t_max

    l = ll(t)
    lp = np.gradient(l, da, edge_order=2)
    lpp = np.gradient(lp, da, edge_order=2)

    plt.plot(alpha, lam0(alpha, t, M=20), label="$l$")
    # plt.plot(alpha,dlam0_da(alpha,t,M=4),label="$l'$")
    # plt.plot(alpha,dlam02_da2(alpha,t,M=4),label="$l''$")

    plt.plot(alpha, l, label="$nl$")
    plt.plot(alpha, lp, label="$nl'$")
    plt.plot(alpha, lpp, label="$nl''$")
    plt.legend()
    plt.show()

    plt.spy(D2)
    plt.show()
    plt.close()


def checkDerivatives():
    f = a ** 5

    plt.plot(a, 5 * a ** 4, label=r"$f'$")
    plt.plot(a, D1.dot(f), "--", label=r"$\tilde{f}'$")
    plt.legend()
    plt.show()
    plt.close()

    plt.plot(a, 20 * a ** 3, label=r"$f''$")
    plt.plot(a, D2.dot(f), "--", label=r"$\tilde{f}''$")
    plt.plot(a, d2dx(f, da), "--", label=r"$\hat{f}''$")
    plt.legend()
    plt.show()
    plt.close()

    plt.plot(a, 60 * a ** 2, label=r"$f'''$")
    plt.plot(a, D3.dot(f), "--", label=r"$\tilde{f}'''$")
    plt.legend()
    plt.show()
    plt.close()

    plt.plot(a, 120 * a, label=r"$f''''$")
    plt.plot(a, D4.dot(f), "--", label=r"$\tilde{f}''''$")
    plt.legend()
    plt.show()
    plt.close()


def makeD1(dr, N):
    ldd = np.zeros(N - 2)
    ld = -np.ones(N - 1)
    dd = np.zeros(N)
    ud = np.ones(N - 1)
    udd = np.zeros(N - 2)

    # top row
    dd[0] = -15 / 4
    ud[0] = 4
    udd[0] = -0.25

    # bottom row
    dd[-1] = 15 / 4
    ld[-1] = -4
    ldd[-1] = 0.25

    return diags([ldd, ld, dd, ud, udd],
                 offsets=[-2, -1, 0, 1, 2],
                 shape=(N, N)) * (0.5 / dr)


def makeD2(dr, N):
    ldd = np.zeros(N - 2)
    ld = np.ones(N - 1)
    dd = np.full(N, -2.0)
    ud = np.ones(N - 1)
    udd = np.zeros(N - 2)

    # top row
    dd[0] = 1.75
    ud[0] = -2
    udd[0] = 0.25

    # bottom row
    ldd[-1] = 0.25
    ld[-1] = -2
    dd[-1] = 1.75
    return diags([ldd, ld, dd, ud, udd],
                 format="csr",
                 offsets=[-2, -1, 0, 1, 2],
                 shape=(N, N)) * (1 / dr ** 2)


def makeD3(dr, N, bc='Fily'):
    if bc == "Fily":
        # Make arrays
        ldd = -np.ones(N - 2)
        ld = 2 * np.ones(N - 1)
        dd = np.zeros(N)
        ud = -2 * np.ones(N - 1)
        udd = np.ones(N - 2)

        dd[0] = -45 / 2
        dd[1] = 3
        dd[-1] = 45 / 2
        dd[-2] = -3

        ld[0] = -7 / 4
        ld[-2] = 9 / 4
        ld[-1] = -24

        ud[0] = 24
        ud[1] = -9 / 4
        ud[-1] = 7 / 4

        ldd[-1] = 3 / 2

        udd[0] = -3 / 2

        D3 = diags([ldd, ld, dd, ud, udd],
                   offsets=[-2, -1, 0, 1, 2],
                   format="csr",
                   shape=(N, N)) * (0.5 / dr ** 3)

    elif bc == "Shelley":

        # Make arrays
        ldd = -np.ones(N - 2)
        ld = 2 * np.ones(N - 1)
        dd = np.zeros(N)
        ud = -2 * np.ones(N - 1)
        udd = np.ones(N - 2)

        # row 1:
        dd[0] = 0
        ud[0] = 0
        udd[0] = 0

        # row 2:
        dd[1] = 0
        ud[1] = 0
        udd[1] = 0
        ld[0] = 0

        # row N-2:
        dd[-2] = 0
        ld[-2] = 0
        ldd[-2] = 0
        ud[-1] = 0

        # row N-1:
        dd[-1] = 0
        ld[-1] = 0
        ldd[-1] = 0

        D3 = diags([ldd, ld, dd, ud, udd],
                   offsets=[-2, -1, 0, 1, 2],
                   format="csr",
                   shape=(N, N)) * 0.5

        row = np.array([0, 0, 0, 0, 0,
                        1, 1, 1, 1, 1,
                        N - 2, N - 2, N - 2, N - 2, N - 2,
                        N - 1, N - 1, N - 1, N - 1, N - 1])

        col = np.array([0, 1, 2, 3, 4,
                        0, 1, 2, 3, 4,
                        N - 1, N - 2, N - 3, N - 4, N - 5,
                        N - 1, N - 2, N - 3, N - 4, N - 5])

        data = np.array([-5 / 2, 9, -12, 7, -3 / 2,
                         -3 / 2, 5, -6, 3, -1 / 2,
                         3 / 2, -5, 6, -3, 1 / 2,
                         5 / 2, -9, 12, -7, 3 / 2])

        D3bcs = scipy.sparse.csr_matrix((data, (row, col)), shape=(N, N))

        D3 += D3bcs
        D3 /= dr ** 3

    else:
        print("Please check bcs are either Fily or Shelley. Exiting...")
        quit()

    D3.eliminate_zeros()
    return D3


def makeD4(dr, N, bc="Fily"):
    # Make arrays

    if bc == "Fily":

        ldd = np.ones(N - 2)
        ld = -4 * np.ones(N - 1)
        dd = 6 * np.ones(N)
        ud = -4 * np.ones(N - 1)
        udd = np.ones(N - 2)

        dd[0] = 21
        dd[1] = 3
        dd[-2] = 3
        dd[-1] = 21

        ld[0] = -1 / 4
        ld[-2] = -15 / 4
        ld[-1] = -24

        ud[-1] = -1 / 4
        ud[1] = -15 / 4
        ud[0] = -24

        ldd[-1] = 3

        udd[0] = 3

        D4 = diags([ldd, ld, dd, ud, udd],
                   offsets=[-2, -1, 0, 1, 2],
                   format="csr",
                   shape=(N, N)) * (1 / dr ** 4)

    elif bc == "Shelley":

        # Make arrays
        ldd = np.ones(N - 2)
        ld = -4 * np.ones(N - 1)
        dd = 6 * np.ones(N)
        ud = -4 * np.ones(N - 1)
        udd = np.ones(N - 2)

        # row 1:
        dd[0] = 0
        ud[0] = 0
        udd[0] = 0

        # row 2:
        dd[1] = 0
        ud[1] = 0
        udd[1] = 0
        ld[0] = 0

        # row N-2:
        dd[-2] = 0
        ld[-2] = 0
        ldd[-2] = 0
        ud[-1] = 0

        # row N-1:
        dd[-1] = 0
        ld[-1] = 0
        ldd[-1] = 0

        D4 = diags([ldd, ld, dd, ud, udd],
                   offsets=[-2, -1, 0, 1, 2],
                   format="csr",
                   shape=(N, N))

        row = np.array([0, 0, 0, 0, 0, 0,
                        1, 1, 1, 1, 1, 1,
                        N - 2, N - 2, N - 2, N - 2, N - 2, N - 2,
                        N - 1, N - 1, N - 1, N - 1, N - 1, N - 1])

        col = np.array([0, 1, 2, 3, 4, 5,
                        0, 1, 2, 3, 4, 5,
                        N - 1, N - 2, N - 3, N - 4, N - 5, N - 6,
                        N - 1, N - 2, N - 3, N - 4, N - 5, N - 6])

        data = np.array([3, -14, 26, -24, 11, -2,
                         2, -9, 16, -14, 6, -1,
                         2, -9, 16, -14, 6, -1,
                         3, -14, 26, -24, 11, -2
                         ])

        D4bcs = scipy.sparse.csr_matrix((data, (row, col)), shape=(N, N))

        D4 += D4bcs

        D4 /= dr ** 4

    else:
        print("Please check bcs are either Fily or Shelley. Exiting...")
        quit()

    D4.eliminate_zeros()
    # pprint(D4.toarray())
    return D4


def getDs(da, N, bc="Fily"):
    D1 = makeD1(da, N)
    D2 = makeD2(da, N)
    D3 = makeD3(da, N, bc=bc)
    D4 = makeD4(da, N, bc=bc)
    return D1, D2, D3, D4


def setUpMatrix(ll, a, t, chi=1, bc='Fily'):
    N = len(a)
    da = a[1] - a[0]

    l = ll(t)
    lp = ddx(l, da)
    lpp = d2dx(l, da)
    # lp=np.gradient(l,da,edge_order=2)
    # lpp=np.gradient(lp,da,edge_order=2)

    a4 = -(chi / l ** 2) * np.exp(-4 * t)
    A4 = diags([a4], offsets=[0], format='csr', shape=(N, N))

    a3 = ((2 * chi * lp) / l ** 3) * np.exp(-4 * t)
    A3 = diags([a3], offsets=[0], format='csr', shape=(N, N))

    a2 = (
                 (chi / l) * (lpp / l ** 2 - (2 * lp ** 2) / l ** 3) * np.exp(-2 * t)
                 + (l - 1) / l
         ) * np.exp(-2 * t)
    A2 = diags([a2], offsets=[0], format='csr', shape=(N, N))

    a1 = ((2 * lp) / l) * np.exp(-2 * t)
    A1 = diags([a1], offsets=[0], format='csr', shape=(N, N))

    D1, D2, D3, D4 = getDs(da=da, N=N, bc=bc)
    L = A4.dot(D4) + A3.dot(D3) + A2.dot(D2) + A1.dot(D1)
    return L


def solveEigenProblemAtT(ll, a, t, kappa, beta, bc='Fily'):
    _chi = beta / kappa
    L = setUpMatrix(ll, a, t, chi=_chi, bc=bc)
    w, v = np.linalg.eig(L.toarray())
    inds = np.argpartition(np.real(w), -2)[-2:]
    inds = inds[np.argsort(np.real(w[inds]))][::-1]  # not necessarily ordered
    w = w[inds]
    v = v[inds]
    return w


def solveEigenProblem(ll, a, t_max, chi, bc='Fily'):
    min_ws = []
    ts = np.linspace(0, t_max, 200)

    for t in ts:
        w = solveEigenProblemAtT(ll, a, t, chi=chi, bc=bc)
        min_ws.append(w)
    return min_ws, ts


def findTCritRand(ll, a, t_max, chi, bc='Fily', num_attempts=10):
    f = lambda tc: np.real(solveEigenProblemAtT(ll, a, tc, chi=chi, bc=bc)[1])
    roots = []
    guesses = np.random.random(size=num_attempts) * t_max
    for guess in guesses:
        sol = root_scalar(f, bracket=[0, t_max], x0=guess)
        roots.append(sol.root)
    print(f"{chi=} {np.mean(roots)=} {np.std(roots)=}")
    return np.mean(roots), np.std(roots)


def findTCrit(kappa, beta, N, a, t_max, bc='Fily'):
    t0 = timer()
    ll = getLam0(a=a, k=kappa, N=N, t_max=t_max, method="BDF")
    print(f"Took {timer() - t0}s to find lambda")

    def getFirstPositive(tc):
        max_ev = np.real(solveEigenProblemAtT(ll=ll, a=a, t=tc, kappa=kappa, beta=beta, bc=bc)[1])
        # print(max_ev)
        return kappa*max_ev

    try:
        sol = root_scalar(getFirstPositive, bracket=[0, t_max], x0=0.75 * t_max)
    except Exception as e:
        print(e)
        ts = np.linspace(0, t_max)
        minw = []
        for t in ts:
            w = solveEigenProblemAtT(ll=ll, a=a, t=t, kappa=kappa, beta=beta, bc=bc)
            minw.append(w)
        plotEigenvalues(ts, np.array(minw))
        quit()

    # print(sol)
    print(f"{kappa=} chi={beta / kappa} tau_c: {sol.root=}")
    return sol.root, 0.0


def plotEigenvalues(ts, mws, t_crit=None, fname="bucklingEigenvaluesPython.pdf", show=True):
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

    if t_crit:
        # ax.axvline(t_crit,linestyle=':',color='k',label=r"$\tau_c$")
        ax.axvline(t_crit, linestyle=':', color='k', ymax=0.5)

        txt = ax.text(0.63, 0.25, rf"$\tau_c={t_crit:3.3f}$",
                      horizontalalignment='center',
                      verticalalignment='center', transform=ax.transAxes)
    plt.legend(loc="upper left", borderpad=0.2, borderaxespad=0.2, handletextpad=0.5,
               ncol=2, columnspacing=1)

    ax.set_xlabel(r"$\tau_0$")
    ax.set_ylabel(r"$\omega$")

    fig.savefig(fname, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def plotLogEigenvalues(ts, mws, tcrit=None, fname="bucklingEigenvaluesPythonLog.pdf", show=True):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.semilogy(ts, np.abs(np.real(mws[:, 0])), "b-", label="$|\Re(\omega_1)|$")
    ax.semilogy(ts, np.abs(np.imag(mws[:, 0])), "--", label="$|\Im(\omega_1)|$")
    ax.semilogy(ts, np.abs(np.real(mws[:, 1])), "m-", label="$|\Re(\omega_2)|$")
    ax.semilogy(ts, np.abs(np.imag(mws[:, 1])), ":", label="$|\Im(\omega_2)|$")

    # axins = inset_axes(ax, width="100%", height="100%",
    #                bbox_to_anchor=(.35, .15, .4, .4),
    #                bbox_transform=ax.transAxes, loc=3)
    #
    # axins.loglog(ts[1:],np.abs(np.real(mws[:,0])[1:]),
    #              "b-",label="$|\Re(\omega_1)|$")
    # axins.loglog(ts[1:],np.abs(np.real(mws[:,1])[1:]),
    #              "m-",label="$|\Re(\omega_2)|$")
    # axins.set_xscale("log")
    # axins.set_yscale("log")
    # axins.set_yticks([],minor=True)
    # axins.set_yticks([2.4,3.0,3.6])
    # axins.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    if tcrit:
        ax.axvline(tcrit, linestyle=':', color='k', label=r"$\tau_c$")
    plt.legend(loc="best")

    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel("$\omega/\kappa$")

    analysis_dir = f"../GeneratedOutput/{savedir}"
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    save_name = f"{analysis_dir}/{fname}"
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def plotAllEigs(tcrit, chi):
    minw, ts = solveEigenProblem(ll, a, t_max, chi)
    minw = np.array(minw)
    plotEigenvalues(ts, minw, tcrit,
                    fname=f"../GeneratedOutput/{savedir}/bucklingEigenvaluesPython_{bc}_{N=}_tc_{k=:6.6f}_{chi=:6.6f}.pdf",
                    show=False
                    )
    plotLogEigenvalues(ts, minw, tcrit,
                       fname=f"bucklingLogEigenvaluesPython_{bc}_{N=}_tc_{k=:6.6f}_{chi=:6.6f}.pdf",
                       show=False
                       )


def ddx(y, dr):
    dydx = np.zeros(len(y))

    # Interior derivatives
    dydx[1:-1] = y[2:] - y[:-2]

    # Edge derivatives
    dydx[0] = -3 * y[0] + 4 * y[1] - y[2]
    dydx[-1] = 3 * y[-1] - 4 * y[-2] + y[-3]

    return dydx * 0.5 / dr


def d2dx(y, dr):
    d2ydx = np.zeros(len(y))

    # Interior derivatives
    d2ydx[1:-1] = (y[:-2] - 2 * y[1:-1] + y[2:])

    # Edge derivatives
    d2ydx[0] = 2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]
    d2ydx[-1] = 2 * y[-1] - 5 * y[-2] + 4 * y[-3] - y[-4]

    return d2ydx / dr ** 2


def dl_dt(t, y, dr, k, J):
    dydt = k * np.exp(-2 * t) * d2dx(y, dr) - y
    dydt[0] = 0
    dydt[-1] = 0
    return dydt


def getJac(t, y, dr, k, J):
    jac = J * (k * np.exp(-2 * t) / dr ** 2)

    return jac


def makeJ(N):
    ld = np.ones(N - 1)
    dd = np.full(N, -2) - 1
    ud = np.ones(N - 1)
    dd[0] = 0
    ud[0] = 0
    ld[-1] = 0
    dd[-1] = 0

    J = diags([ld, dd, ud],
              offsets=[-1, 0, 1],
              shape=(N, N))
    return J


def getLam0(a, k=1, N=501, t_max=7, method='BDF'):
    da = a[1] - a[0]
    N = len(a)
    ic = np.ones(N)
    J = makeJ(N)
    sol = solve_ivp(dl_dt,
                    t_span=[0, t_max],
                    y0=ic,
                    jac=getJac,
                    atol=1e-8,
                    rtol=1e-6,
                    dense_output=True,
                    args=(da, k, J),
                    method=method)

    return sol.sol


def sweepFn(ll, a, t_max, num=31):
    chis = np.logspace(-8, 0, num, base=10, endpoint=True)

    # ll,a=getLam0(k=k,N=N,t_max=t_max)

    tcrits = []
    terrs = []

    for cc in tqdm(chis):
        tcrit, t_err = findTCrit(ll, a, t_max, chi=cc)
        plotAllEigs(tcrit, chi=cc)
        tcrits.append(tcrit)
        terrs.append(t_err)

    return chis, np.array(tcrits), np.array(terrs)


def plotSweepResults(chis, tcrits, terrs=None, fname="tcs.pdf", show=True):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    # try:
    #     p1=ax.errorbar(chis,tcrits,terrs,marker='o',mfc="w",mec='k',
    #                 capsize=5,alpha=0.7)
    # except:
    p1 = ax.plot(l0 ** 2 * chis, tcrits, marker='o', mfc="w", mec='k', alpha=0.7,
                 label=rf"$\tau_c(\chi)$")

    # axins = inset_axes(ax, width="100%", height="100%",
    #                bbox_to_anchor=(.35, .15, .55, .55),
    #                bbox_transform=ax.transAxes, loc=3)
    #
    # ip1=axins.plot(chis,tcrits)
    # f = lambda x,a,b: a*x**b
    # popt, pcov = curve_fit(f,chis,tcrits,
    #                        p0=[0.5,2.4],bounds=(0,4))
    # print(popt)
    # ip2=axins.plot(chis,f(chis,*popt),"--")
    # axins.set_xscale("log")
    # axins.set_yscale("log")
    # axins.set_yticks([],minor=True)
    # axins.set_yticks([2.4,3.0,3.6])
    # axins.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    # af=f"{popt[0]:2.2f}"
    # bf=f"{popt[1]:2.2f}"
    # p2=ax.plot(chis,f(chis,*popt),"--",zorder=3,
    #            label=r"${a}\chi^{b}$,"+f"\n$a={af}$,\n$b={bf}$")
    ax.set_xlabel("$B/K$ (\si{\metre^2})")
    ax.set_ylabel(r"$\tau_c$")
    # plt.legend()

    analysis_dir = f"../GeneratedOutput/{savedir}"
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    save_name = f"{analysis_dir}/{fname}"
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()

    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.plot(chis, np.exp(tcrits))

    ax.set_xlabel("$\chi$")
    ax.set_ylabel(r"$\ell_c/\ell_0$")

    analysis_dir = f"../GeneratedOutput/{savedir}/"
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)

    save_name = analysis_dir + fname.replace("tcs_", "ells_")
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    if show:
        plt.show()
    plt.close()


def calc_lambda_analytic(a, t, k, M=11):
    def calc_pnt(n):
        try:
            with np.errstate(invalid='raise'):
                kn = (np.pi * n) ** 2 * k
                emt = np.exp(-t)
                large_number = 0.5 * kn * emt ** 2
                le = emt * np.exp(large_number)

                p0 = - 4 / (n * np.pi)
                p1 = -p0 * emt * np.exp(0.5 * kn * (emt ** 2 - 1.0))
                erf1 = erf(np.sqrt(0.5 * kn))
                p2 = (2.0 / n) * le * np.sqrt(2.0 * kn / np.pi) * (erf1 - erf(emt * np.sqrt(0.5 * kn)))
            return p0 + p1 + p2
        except OverflowError as e:
            print(e)
            return 0.0

    lam = 1 + np.sum([calc_pnt(n) * np.sin(np.pi * a * n) for n in range(1, M, 2)], axis=0)
    return lam

def calc_lambda_analytic_low_k(a, t, k, M=11):
    def calc_pnt(n):
        p0 = lambda n: -(4/(n*np.pi)) * (1 - np.exp(-t))
        p1 = lambda n: (2 * np.pi * n) * np.exp(3*t) * (np.exp(t)-1)**2 * k
        return p0(n) + p1(n)

    lam = 1 + np.sum([calc_pnt(n) * np.sin(np.pi * a * n) for n in range(1, M, 2)], axis=0)
    return lam

def plotLambda(ll, a, t_max=7, M=11, kappa=None):
    fig, ax = plt.subplots(1, 1, figsize=setFigSize(247))
    ax.plot(a, ll(0), label=r"$t={0}$")
    ax.plot(a, ll(t_max / 2), label=r"$t=\frac{1}{2}\ln(\kappa)$")
    ax.plot(a, ll(t_max), label=r"$t=\ln(\kappa)$")
    if kappa is not None:
        # ax.plot(a, lam0(a, t_max, kappa, L=1.0, M=M), ':o', alpha=0.5, label='analytic_old')
        ax.plot(a, calc_lambda_analytic(a, t_max, kappa, M=M), '--', alpha=0.7, label='analytic_new')
        # ax.plot(a, calc_lambda_analytic_low_k(a, t_max, kappa, M=101), '--', alpha=0.7, label='analytic_low')
    plt.legend()
    plt.title(f"{kappa=}")
    save_name = "../GeneratedOutput/check_lambda.pdf"
    fig.savefig(save_name, format="pdf", bbox_inches='tight',
                transparent=True)
    plt.show()
    plt.close()


if __name__ == "__main__":

    d = 1  # diameter (micron)
    g = 0.75 / np.log(2)  # growth rate (per hr)
    Y0 = 4e6  # Young's modulus (Pa)
    Y = Y0
    l0 = 2  # initial length (microns)
    xi = 200  # friction coefficient (Pa hr)
    chi = (d / l0) ** 2  # beta / kappa
    k = (Y0 * d ** 2) / (g * xi * l0 ** 2)  # increased k to account for the spring k=(Y*d**2)/(g*xi*l0**2)
    # t_max=np.log(k)
    t_max = 7
    N = 201
    a, da = np.linspace(0, 1, N, retstep=True, endpoint=True)
    savedir = "frozen_mode_extensible/"
    # D1 = makeD1(da, N)
    # D2 = makeD2(da, N)
    # D3 = makeD3(da, N, bc=bc)
    # D4 = makeD4(da, N, bc=bc)

    print(f"{k=:6.6f} {chi=:6.6f} {t_max=:3.3f}")

    # checkDerivatives()
    # quit()

    t0 = timer()
    ll = getLam0(a=a, k=k, N=N, t_max=t_max)
    print(f"Took {timer() - t0:3.3f}s to find ll")

    plotLambda(ll)

    minw, ts = solveEigenProblem(ll, a, t_max, chi, bc="Fily")
    minw = np.array(minw)
    print(f"The maximum real parts are: {max(np.real(minw[:, 0]))} and {max(np.real(minw[:, 1]))}")
    tcrit, terr = findTCrit(ll, a, t_max, chi)

    plotEigenvalues(ts, minw, tcrit,
                    fname=f"../GeneratedOutput/{savedir}/bucklingEigenvaluesPython_"
                          f"{N=}_{bc}_tc_{k=:6.6f}_{chi=:6.6f}.pdf")

    dist_dir = f"../GeneratedOutput/SimOutput/LcDist/"
    dist_analysis_dir = f"../GeneratedOutput/figures/LcDist/"
    tc_data_fname = f"tcs_data_{N=}_{bc}_{k=:6.6f}_{chi=:6.6f}.txt"
    if not os.path.exists(f"../GeneratedOutput/{savedir}/{tc_data_fname}"):
        chis, tcrits, terrs = sweepFn(ll, a, t_max, num=31)
        np.savetxt(f"../GeneratedOutput/{savedir}/{tc_data_fname}", np.transpose([chis, tcrits, terrs]))
    else:
        data = np.loadtxt(f"../GeneratedOutput/{savedir}/{tc_data_fname}")
        chis = data[:, 0]
        tcrits = data[:, 1]
        terrs = data[:, 2]

    plotDistributions(dist_dir, dist_analysis_dir, tc=tcrit / g, lc=l0 * np.exp(tcrit / g))

    tc_fig_fname = f"tcs_{N=}_{bc}_{k=:6.6f}_{chi=:6.6f}.pdf"
    plotSweepResults(chis, tcrits, terrs=None, fname=tc_fig_fname, show=False)

    plotEigenvalues(ts, minw, tcrit,
                    fname=f"../GeneratedOutput/{savedir}/bucklingEigenvaluesPython_{N=}_{bc}_tc_{k=:6.6f}_{chi=:6.6f}.pdf")
    plotLogEigenvalues(ts, minw, tcrit, fname=f"bucklingLogEigenvaluesPython_{N=}_{bc}_tc_{k=:6.6f}_{chi=:6.6f}.pdf")
