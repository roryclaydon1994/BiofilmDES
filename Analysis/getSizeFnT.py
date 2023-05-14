# Standard modules
import re
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pandas as pd
import os
import scipy.optimize as scopt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, zoomed_inset_axes

# Third party

# Custom modules
import utilities as ut

ut.setMPL()
save_dir = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/figures/sizes_t/"

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1, 1))


def sort_by_time(_file):
    """
    Extract the time index from the file name
    :param _file:
    :return:
    """
    res = re.findall("(?<=_)\d+", _file)
    try:
        return int(res[0])
    except IndexError:
        return np.inf
    except Exception as _e:
        print(_e)
        quit()


def getGrowthRates(_dir):
    """
    :param _dir:
    :return:
    """
    _dirs = glob(f"{_dir}/repeat*")
    growth_rates = pd.Series(name='grwth_rate', dtype=np.float64)
    avgs = []
    all_data = []
    for c, tmp_dir in enumerate(_dirs):
        print(tmp_dir)
        sort_files = sorted(os.listdir(tmp_dir), key=sort_by_time)
        for r, file in enumerate(sort_files):
            if 'dat' not in file:
                continue
            # if r % 20 != 0:
            #     continue
            if "200" not in file:
                continue
            # if r < 100:
            #     continue

            try:
                dat = pd.read_csv(f"{tmp_dir}/{file}", sep='\t')
                if len(dat) > 100:
                    avgs.append(dat['grwth_rate'].mean())
                    all_data.append(dat['grwth_rate'].to_numpy() / 0.5e-4)
                    growth_rates = pd.concat([growth_rates, dat['grwth_rate']], ignore_index=True)
            except Exception as e:
                print(f"{e} for {file}")
                quit()
    return growth_rates, np.array(avgs), all_data


def getData(_dirs):
    """

    :param _dirs:
    :return:
    """

    cols = len(_dirs)
    rows = max([len(os.listdir(_dir)) for _dir in _dirs]) - 1
    # lengths = np.zeros((rows, cols))
    nums = np.zeros((rows, cols))

    for c, _dir in enumerate(_dirs):
        print(_dir)
        sort_files = sorted(os.listdir(_dir), key=sort_by_time)
        for r, file in enumerate(sort_files):
            if 'log' in file:
                continue

            try:
                dat = pd.read_csv(f"{_dir}/{file}", sep='\t')
                # dat['length'] += 1
                # _length = dat['length'].sum()
                # lengths[r, c] = _length
                nums[r, c] = len(dat)
            except Exception as e:
                print(f"{e} for {file}")
    return nums


def fitGrowthRateData(mu_0=4):
    """

    :param mu_0:
    :return:
    """
    mu_th = mu_0 / np.log(3)

    dat, avgs, all_data = getGrowthRates("/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput/SimOutput"
                                         "/GrowthRateTest/run0/")
    dat /= 0.5e-4
    avgs /= 0.5e-4
    mean = np.mean(avgs)
    std = np.std(avgs) / np.sqrt(len(avgs))
    mean_mu = 0.5 * (dat.max() + dat.min())
    print(f"Mean growth rate from data {mean} pm {std}  theoretical {mu_th}")

    bins = np.linspace(0.5 * mu_0, 1.5 * mu_0, 21, endpoint=True)
    # _, bins = np.histogram(dat, bins=300, density=True)

    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(1))
    avg_mu_sim = ax.axvline(mean, linestyle=':', color='m', alpha=0.8, label=r"$\langle \mu \rangle^{\mathrm{sim}}$")
    ax.axvspan(mean - std, mean + std, alpha=0.2, color='gray')
    avg_mu_th = ax.axvline(mu_th, linestyle='-.', color='k', alpha=0.8, zorder=5,
                           label=r"$\langle \mu \rangle^{\mathrm{th}}$")
    # mus = np.linspace(dat.min(), dat.max())
    mus = 0.5 * (bins[1:] + bins[:-1])

    # for avg in avgs:
    #     plt.axvline(avg, linestyle=':', color='r')
    # plt.show(block=False)

    ns = []
    for datum in all_data:
        _n, _mu = np.histogram(datum, bins=bins, density=True)
        ns.append(_n)
    ns = np.array(ns)

    n_mean = ns.mean(axis=0)

    stds = ns.std(axis=0)

    mu_h = 1.5 * mu_0
    mu_l = 0.5 * mu_0
    # p = lambda x, m, c: (m * x + c) / (0.5 * m * (mu_h ** 2 - mu_l ** 2) + c * (mu_h - mu_l))
    p = lambda x, m, c: (m * x + c)
    popt, pcov = scopt.curve_fit(p, mus, n_mean)
    print(popt, np.sqrt(np.diag(pcov)))

    linear_mean = lambda m, c: ((m / 3) * (mu_h ** 3 - mu_l ** 3) + 0.5 * c * (mu_h ** 2 - mu_l ** 2)) / (
            0.5 * m * (mu_h ** 2 - mu_l ** 2) + c * (mu_h - mu_l))
    avg_mu_lin = plt.axvline(linear_mean(popt[0], popt[1]), dashes=[3, 1, 1, 1, 1, 3, 1], color='c', alpha=0.7,
                             label=r"$\langle \mu \rangle^{\mathrm{lin}}$")
    print(f"linear fit mean : {linear_mean(popt[0], popt[1])}")

    # popt, pcov = scopt.curve_fit(lambda x, a: a / x, mus, n_mean)
    # print(popt, np.sqrt(np.diag(pcov)), popt[0]/np.log(3))
    # plt.plot(mus, popt[0]/mus)
    p_sim, = ax.plot(mus, n_mean, '-o', mfc='w', markevery=0.1, zorder=3, alpha=0.8,
                     label=r"$\mathcal{P}^{\mathrm{sim}}$")
    p_th, = ax.plot(mus, 1.0 / (mus * np.log(3)), zorder=6, color='k', linestyle='--',
                    label=r"$\mathcal{P}^{\mathrm{th}}$")
    p_lin, = plt.plot(mus, p(mus, *popt), '-x', label=r'$\mathcal{P}^{\mathrm{lin}}$', markevery=0.1)

    print(p(mus[0], *popt), p(mus[-1], *popt))

    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.fill_between(mus, n_mean - stds, n_mean + stds, alpha=0.2, color='gray', zorder=4)
    ax.set_xlim([mus[0], mus[-1]])
    ax.set_xlabel(r"$\mu$ (\si{\micro\metre\per\hour})")
    ax.set_ylabel(r"$\mathcal{P}\qty(\mu)$")

    plt.legend(handles=[avg_mu_sim, avg_mu_th, avg_mu_lin, p_sim, p_th, p_lin],
               labels=[r"$\langle \mu \rangle^{\mathrm{sim}}$",
                       r"$\langle \mu \rangle^{\mathrm{th}}$",
                       r"$\langle \mu \rangle^{\mathrm{lin}}$",
                       r"$\mathcal{P}^{\mathrm{sim}}$",
                       r"$\mathcal{P}^{\mathrm{th}}$",
                       r'$\mathcal{P}^{\mathrm{lin}}$'],
               loc='upper right',
               bbox_to_anchor=(1, 1), borderpad=0,
               columnspacing=1.0, handletextpad=0.4,
               borderaxespad=0.4, ncol=2)
    fig.savefig(f"{save_dir}/growth_rate.pdf",
                transparent=False,
                bbox_inches='tight'
                )
    plt.show()


def get_size_fn_t(_dirs):
    """

    :param _dirs:
    :return:
    """
    try:
        nums = np.loadtxt(f"{save_dir}/nums.txt")
    except FileNotFoundError:
        nums = getData(_dirs)
        np.savetxt(f"{save_dir}/nums.txt", nums)
    except Exception as e:
        print(e)
        quit()

    avg_n = np.array(list(row.mean() for row in nums if np.all(row > 0)))
    std_n = np.array(list(row.std() for row in nums if np.all(row > 0)))

    mu_0 = 4
    mu_th = mu_0 / np.log(3)

    # dat = getGrowthRates("/home/rory/PhD_sync/CUDA/cell_linked_list/output/")

    def t_12(mu):
        return 6.0 / (2.0 * mu)

    t = np.arange(0, avg_n.size) * 0.1

    def fit_mu(x, mu):
        lcl_t12 = t_12(mu)
        return 2 ** (x / lcl_t12)

    n2s = std_n > 0
    n0 = avg_n[np.argmax(n2s)]
    t0 = t[np.argmax(n2s)]
    popt, pcov = scopt.curve_fit(fit_mu, t[n2s] - t0, avg_n[n2s] / n0)
    print(f"mu={popt[0]} mu_th={mu_th} {100 * abs(popt[0] - mu_th) / mu_th}")

    t = t[n2s] - t0
    avg_n = avg_n[n2s] / n0
    std_n = std_n[n2s] / np.sqrt(n0)
    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(1))
    ax.plot(t, avg_n, zorder=1, marker='o', alpha=0.7, mec='k', mfc='w', markevery=0.2, label="sim")
    # ax.plot(t, fit_mu(t, *popt), label="fit")
    ax.fill_between(t, avg_n - std_n, avg_n + std_n, color='gray', alpha=0.2)
    ax.plot(t, 2 ** (t / t_12(mu_th)), 'k--',
            label=r"$\left\langle N  \right\rangle\qty(t;\langle \mu \rangle^{\mathrm{th}})$")
    ax.plot(t, 2 ** (t / t_12(mu_0)), ':',
            label=r"$\left\langle N \right\rangle\qty(t;\mu_0)$")
    ax.plot(t, 2 ** (t / t_12(17 * mu_0 / 18)), 'm-.', alpha=0.8,
            label=r"$\left\langle N \right\rangle\qty(t;\langle \mu \rangle^{\mathrm{lin}})$")

    # axi = inset_axes(ax, width="40%", height="40%", loc='lower right', borderpad=0.5)
    axi = zoomed_inset_axes(ax, zoom=2.5, loc='lower right', borderpad=0.2)
    # axi.set_xscale('log')
    # axi.set_yscale('log')
    axi.plot(t, avg_n, zorder=1, marker='o', alpha=0.7, mec='k', mfc='w', markevery=0.2)
    # ax.plot(t, fit_mu(t, *popt), label="fit")
    axi.fill_between(t, avg_n - std_n, avg_n + std_n, color='gray', alpha=0.2)
    axi.plot(t, 2 ** (t / t_12(mu_th)), 'k--')
    axi.plot(t, 2 ** (t / t_12(mu_0)), ':')
    axi.plot(t, 2 ** (t / t_12(17 * mu_0 / 18)), 'm-.', alpha=0.8)
    axi.set_xlim(7.5, np.max(t))
    axi.set_ylim(520, np.max(2 ** (t / t_12(mu_0))))
    axi.tick_params(labelleft=False, labelbottom=False)
    axi.set(xticks=[], yticks=[])
    axi.get_xaxis().set_ticks([])
    axi.get_yaxis().set_ticks([])
    mark_inset(ax, axi, loc1=2, loc2=1, fc="none", ec="0.5", alpha=0.5)

    n30 = np.argmax(avg_n > 30)
    t30 = t[n30]
    print((2 ** (t30 / t_12(mu_th)) - avg_n[n30]) / avg_n[n30])

    # ax.plot(t, fit_mu(t, *popt), ':')
    ax.set_ylabel(r"$\langle N \rangle$")
    ax.set_xlabel(r"$t$ (\si{\hour})")
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.legend(borderpad=0, columnspacing=1.0, handletextpad=0.4, borderaxespad=0.4)
    fig.savefig(f"{save_dir}/sizes.pdf",
                transparent=False,
                bbox_inches='tight'
                )
    plt.show()


if __name__ == "__main__":
    file_dir = "/home/rory/PhD_sync/BiofilmPhageDES/GeneratedOutput"
    file_dir += "/SimOutput/ChainingBucklingTransition/*/*"
    os.makedirs(save_dir, exist_ok=True)
    # fitGrowthRateData()
    get_size_fn_t(glob(file_dir))
