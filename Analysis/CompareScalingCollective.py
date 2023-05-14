# Standard modules
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint
import scipy
import scipy.optimize as scopt

# Third party
from shapely.geometry import Polygon
from tqdm import tqdm

#Custom modules
import utilities as ut
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
from DistributionFunctions import computeColonyContour,getColonyDensity
ut.setMPL()

def getDF(bdir):
    data=[]
    for dd in glob(f"{bdir}/*"):
        try:
            energies=np.loadtxt(f"{dd}/energies.dat",delimiter='\t',skiprows=1)
            hps=np.loadtxt(f"{dd}/LcDist.log",delimiter='\t',skiprows=1)
            # print(energies.shape)
            # print(np.std(energies,axis=0))
            # print(scipy.stats.sem(energies,axis=0))
            try:
                # row=[*hps,*energies.mean(axis=0),*scipy.stats.sem(energies,axis=0)]
                row=[*hps,*energies.mean(axis=0),*energies.std(axis=0)]

            except Exception as e:
                print(e)
                row=[*hps,*[ee for ee in energies],*[0 for ee in energies]]
            data.append(row)
            en_names=list(pd.read_csv(f"{dd}/energies.dat",sep='\t',nrows=1).columns)
            log_names=list(pd.read_csv(f"{dd}/LcDist.log",sep='\t',nrows=1).columns)
        except Exception as e:
            print(e)

    df=pd.DataFrame(data,columns=log_names+en_names+[f"{nn}_std" for nn in en_names])
    return df

def getPackingFraction(cells):
    gon=computeColonyContour(cells,add_links=True,ret_gon=True)
    # print(f"{gon.area=}")
    for geom in gon.interiors:
        if geom.length>=0.1*gon.exterior.length:
            inner=geom
    try:
        outer=Polygon(np.array([(z[0],z[1]) for z in zip(*gon.exterior.xy)]))
        inner=Polygon(np.array([(z[0],z[1]) for z in zip(*inner.xy)]))
        annulus=outer.difference(inner)
    except Exception as e:
        # print(e)
        annulus=Polygon(np.array([(z[0],z[1]) for z in zip(*gon.exterior.xy)]))

    print(f"{annulus.area=}")
    density=sum([cell.getArea() for cell in cells])/annulus.area
    return density

def getR2(y,y_fit):
    # residual sum of squares
    ss_res = np.sum((y - y_fit) ** 2)

    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)

    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    return r2

def getSmoothFit(xs,ys):
    n_interior_knots=5
    qs = np.linspace(0, 1, n_interior_knots+2)[1:-1]
    knots = np.quantile(xs, qs)
    # xxs=np.linspace(xs[imin],xs[imax],3*len(xs))
    tck = scipy.interpolate.splrep(xs, ys, t=knots, k=2)
    yys = scipy.interpolate.splev(xs, tck)
    return yys

def getFitToLinearPart(xs,ys,yerrs):
    # imin=np.argmin(ys)
    # plt.close()
    # plt.axvline(xs[imin],linestyle='--')

    # plt.plot(xs,ys,':')
    # n=5
    # smxs=np.convolve(xs,np.full(n,1/n),mode='valid')
    # smys=np.convolve(ys,np.full(n,1/n),mode='valid')
    # plt.plot(smxs,smys,'--')
    # pp=np.argmax(smys[imin:])
    # plt.axvline(smxs[pp+imin],ls='--')

    yys=getSmoothFit(xs,ys)
    imin=np.argmin(yys)
    imax=len(xs)
    data=[]
    if imin+5>=len(xs):
        print("Error, not enough points for fit")
        quit()
    for ii in range(imin+5,len(xs)):
        xxs=xs[imin:ii+1]
        yys=ys[imin:ii+1]
        yyerrs=yerrs[imin:ii+1]
        n=ii+1-imin+1
        m,c,err=getLinearFit(xxs,yys,yyerrs)
        r2=1-getR2(yys,c+m*xxs)
        data.append([r2,m,c,err,ii])

    data=np.array(data)
    print(data.shape)
    ix=np.argmin(data[:,0])
    print(data)
    # plt.close()
    # plt.plot(data[:,-1],data[:,0])
    # plt.show()
    print(f"im={data[ix,-1]} {imin=} {imax=}")



    # tck=scipy.interpolate.splrep(xs[imin:imax+1],ys[imin:imax+1],k=5,w=1/yerrs[imin:imax+1])



    # # plt.plot(xs,yys)
    # # plt.show()
    #
    # yprimes=scipy.interpolate.splev(xs,tck,der=1)
    # ii=np.argmax(yprimes)
    # xrs=xs[ii:]
    # f=scipy.interpolate.interp1d(xrs,yprimes[ii:],kind='linear')
    # sol=scipy.optimize.root_scalar(
    #     f,bracket=[xrs[0],xrs[np.argmin(yprimes[ii:])]],method='brentq'
    #     )
    # imax=np.argmin(np.abs(xs-sol.root))
    # plt.close()
    return data[ix,1],data[ix,2],data[ix,3],int(data[ix,4])

    # plt.axvline(sol.root)
    # plt.show()
    #
    # plt.plot(xrs,yprimes[ii:])
    # plt.plot(xrs,f(xrs))
    # plt.axvline(sol.root)
    # plt.ylim([-0.1,5])
    # plt.show()
    #
    # quit()

    # yys=scipy.interpolate.splev(xxs, tck, der=1)
    # plt.show()
    # quit()
    # data=[]
    # for ii in range(imin+2,len(xs)):
    #     m,c,err=getLinearFit(xs[imin:ii+1],ys[imin:ii+1],yerrs[imin:ii+1])
    #     data.append([ii+1,m,c,err])
    #     xx=xs[imin:ii+1]
    #     plt.plot(xx,c+m*xx,label=rf'$m={m:3.3f} \pm {err:3.3f}$')
    # data=np.array(data)
    # ix=np.argmin(data[:,-1])
    # print(data)
    # print(data[ix])
    # plt.legend(ncol=10)
    # plt.show()
    #
    # plt.errorbar(data[:,0],data[:,1],data[:,-1])
    # plt.show()
    #
    # plt.plot(xs,ys)
    # plt.axvline(xs[imin],linestyle='--')
    #
    # plt.show()
    # quit()

def getLinearFit(xs,ys,yerrs):
    f=lambda x,m,c: m*x+c
    popt,pcov=scopt.curve_fit(
        f,xs,ys,sigma=yerrs,
        p0=[0,0])
    m,c=popt
    errors=np.sqrt(np.diag(pcov))
    return m,c,errors[0]

def getSmoothedGradient(rs,es):
    kf=np.array([-1])
    sm_window=20
    while np.any(kf<=0):
        sm=np.full(sm_window,1/sm_window)
        smoothed_energy=np.convolve(es,sm,mode='valid')
        smoothed_rs=np.convolve(rs,sm,mode='valid')
        # ax.plot(1/smoothed_rs**2,smoothed_energy,'--',label=rf"$p_\ell={lp}^s$")
        kf=np.gradient(smoothed_energy,1/smoothed_rs**2,edge_order=2)
        sm_window*=2
    return kf

def getDependenceLinkingProb(bdir):
    df=getDF(bdir)
    df.sort_values(by=['annulus_radius','LinkingProb'],inplace=True)
    lps=df['LinkingProb'].unique()
    kfs=[]
    # kis=[]
    fit_params=[]
    # f=lambda x,x0,nu,k: k*(x+x0)**nu-k*x0**nu
    f=lambda x,m,c: m*x+c
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=1))
    for ii,lp in enumerate(lps):
        # if lp!=0:
        #     continue
        lpdf=df[df['LinkingProb']==lp]
        if len(lpdf)<=2:
            continue
        print(f"{len(lpdf)=}")
        rs=lpdf['annulus_radius'].to_numpy()
        ws=lpdf['annulus_width'].to_numpy()
        # energies_relaxed=lpdf['normalised_relaxed_delta_E'].to_numpy()
        # stds_relaxed=lpdf['normalised_relaxed_delta_E_std'].to_numpy()

        # energies_initial=lpdf['normalised_initial_delta_E'].to_numpy()
        # stds_initial=lpdf['normalised_initial_delta_E_std'].to_numpy()
        ix=np.argsort(rs)[::-1]
        bend_energy_annulus=lpdf['annulus_channel_initial_bend'].to_numpy()[ix]
        bend_energy_annulus_err=lpdf['annulus_channel_initial_bend_std'].to_numpy()[ix]

        # ix=energies!=np.nan
        # ax.plot(1/rs[ix]**2,energies[ix],label=rf"$p_\ell={lp}$")
        # energies_relaxed=energies_relaxed[ix]
        rs=rs[ix]
        xs=1/rs**2
        # ax.errorbar(1/rs**2,energies_relaxed,stds_relaxed,label=rf"$p_\ell={lp}$")
        # ax.errorbar(
        #     xs,energies_initial,stds_initial,
        #     marker='o',
        #     mfc='None',mec='k',capsize=5,label=rf"$p_\ell={lp}$"
        #     )
        ys=bend_energy_annulus/(2*np.pi*ws*rs)
        if ii%2==0:
            ax.plot(xs,ys,label=rf"$p_\ell={lp}$")
        m,c,m_err=getLinearFit(xs,ys,bend_energy_annulus_err)
        # ax.plot(xs,m*xs+c,'k')
        fit_params.append([lp,m,m_err])
        # kf=getSmoothedGradient(rs,bend_energy_annulus)
        # kfs.append([lp,kf.mean()])

    ax.set_xlabel("$R^{-2}$")
    # ax.set_xlabel("$R$")
    ax.set_ylabel(r"$\varepsilon_\beta$")
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_ylim([0,df['normalised_initial_delta_E'].max()])
    ax.legend(ncol=2)
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/bend_energy_dependence.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

    kfs=np.array(kfs)
    # print(kfs)
    # kis=np.array(kis)
    fit_params=np.array(fit_params)
    pprint(fit_params)
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    # ax.plot(
    #     kfs[:,0],kfs[:,1],
    #     label=r'$\left\langle \pdv{\tilde{\varepsilon}}{\tilde{\mathcal{K}^2}} \right\rangle$'
    #     )
    # ax.plot(kis[:,0],kis[:,1],'o--',label='initial')
    ax.errorbar(
        fit_params[:,0],fit_params[:,1],yerr=fit_params[:,2],
        linestyle=':',marker='o',mec='k',mfc='None',capsize=5,
        label=r'$K_\beta(p_\ell)$'
        # label=r'$\varepsilon(0)+\mathcal{K}^2 \eval{\pdv{\varepsilon}{\mathcal{K}^2}}_{\mathcal{K}=0}$'
        )
    m,c,m_err=getLinearFit(fit_params[:,0],fit_params[:,1],fit_params[:,2])
    ax.plot(fit_params[:,0],m*fit_params[:,0]+c,'k--',
            label=rf"$K_\beta={m:1.2f}p_\ell+{c:1.2f}$")
    print(f"fit to K_beta {m=}")
    ax.set_xlabel(r"$p_\ell$")
    ax.set_ylabel(r"$K_\beta$")
    ax.legend()
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_beta_vs_p_ell.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=True)
    quit()

def getDependenceWidth(bdir,energy_type='simple'):
    df=getDF(bdir)
    df.sort_values(by=['annulus_width','annulus_radius'],inplace=True)
    mod=f'_scal_mathcalK_{energy_type}'
    widths=df['annulus_width'].unique()
    lps=df['LinkingProb'].unique()
    lps.sort()
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=1))
    # use_w = [ ww for ww in widths if len(df[df['annulus_width']==ww])>=801 ]
    # for ww in use_w[::2]:
    for ww in widths:

        # if ww>8:
        #     continue
        wdf=df[ (df['annulus_width']==ww) & (df['LinkingProb']==lps[0]) ]
        rs=wdf['annulus_radius'].to_numpy()
        energies_initial=wdf['normalised_initial_delta_E'].to_numpy()
        ix=np.argsort(rs)[::-1]
        energies_initial=energies_initial[ix]
        rs=rs[ix]
        xs=1/rs
        # ix2=xs<0.05
        ax.plot(
            xs**2,energies_initial,
            marker='o',mfc='None',mec='k',
            label=rf"$w={ww:1.1f}$"
            )
        if energy_type!='simple':
            energies_relaxed=wdf['normalised_relaxed_delta_E'].to_numpy()
            energies_relaxed=energies_relaxed[ix]
            print(f"{energies_relaxed.mean()=}")
            ax.plot(
                xs**2,energies_relaxed,'--',
                marker='o',mfc='None',mec='k',
                label=rf"$w={ww:1.1f}$ relaxed"
                )
            ax.axvline(1e-4*3*4/(3.75**2),linestyle=':')

    ax.set_xlabel("$R^{-2}$")
    ax.set_ylabel(r"$\varepsilon_H$")
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend(ncol=2)
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/width_sweep_{energy_type}.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=1))
    # use_w = [ ww for ww in widths if len(df[df['annulus_width']==ww])>=801 ]
    # for ww in use_w[::2]:
        # if ww>8:
        #     continue
    for ww in widths:
        wdf=df[df['annulus_width']==ww]
        rs=wdf['annulus_radius'].to_numpy()
        energies_initial=wdf['normalised_initial_delta_E'].to_numpy()
        ix=np.argsort(rs)[::-1]
        energies_initial=energies_initial[ix]
        rs=rs[ix]
        xs=1/rs
        # ix2=xs<0.05
        ax.plot(
            xs,energies_initial,
            label=rf"$w={ww:1.1f}$"
            )
        if energy_type!='simple':
            energies_relaxed=wdf['normalised_relaxed_delta_E'].to_numpy()
            energies_relaxed=energies_relaxed[ix]
            print(f"{energies_relaxed.mean()=}")
            ax.plot(
                xs,energies_relaxed,'--',
                label=rf"$w={ww:1.1f}$ relaxed"
                )
    ax.set_xlabel("$R^{-1}$")
    ax.set_ylabel(r"$\varepsilon_H$")
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend(ncol=2)
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/width_sweep{mod}.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

    fit_params=[]
    kfs=[]
    for ww in widths:
        wdf=df[df['annulus_width']==ww]
        print(f"{len(wdf)=}")
        if len(wdf)<2:
            continue

        rs=wdf['annulus_radius'].to_numpy()
        energies_initial=wdf['normalised_initial_delta_E'].to_numpy()
        stds_initial=wdf['normalised_initial_delta_E_std'].to_numpy()

        ix=np.argsort(rs)[::-1]
        energies_initial=energies_initial[ix]
        rs=rs[ix]
        xs=1/rs
        # ax.errorbar(
        #     xs,energies_initial,stds_initial,
        #     marker='o',
        #     mfc='None',mec='k',capsize=5,label=rf"$w={ww:1.2f}$"
        #     )
        m,c,m_err=getLinearFit(xs,energies_initial,stds_initial)
        fit_params.append([ww,m,m_err])
        # kf=getSmoothedGradient(rs,energies_initial)
        # kfs.append([ww,kf.mean()])

    fit_params=np.array(fit_params)
    kfs=np.array(kfs)
    pprint(fit_params)
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    # ax.plot(
    #     kfs[:,0],kfs[:,1],
    #     label=r'$\left\langle \pdv{\tilde{\varepsilon}}{\tilde{\mathcal{K}^2}} \right\rangle$'
    #     )
    ax.errorbar(
        10.5/fit_params[:,0],fit_params[:,1],yerr=fit_params[:,2],
        linestyle=':',marker='o',mec='k',mfc='None',capsize=5,
        label=r'$\varepsilon(0)+\mathcal{K} \eval{\pdv{\varepsilon}{\mathcal{K}}}_{\mathcal{K}=0}$'
        )
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$K_H$")
    ax.legend()
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_w{mod}.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

def plotEnergies(bdir):
    df=getDF(bdir)
    df.sort_values(by=['annulus_width','annulus_radius'],inplace=True)
    mod='_energies'
    widths=df['annulus_width'].unique()
    fit_params=[]
    kfs=[]
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(cols=1))
    w0=widths[0]
    for ww in [w0]:

        wdf=df[df['annulus_width']==ww]
        print(wdf)
        if len(wdf)<=2:
            continue

        rs=wdf['annulus_radius'].to_numpy()
        ix=np.argsort(rs)[::-1]
        rs=rs[ix]
        hertzian_energy_annulus=wdf['annulus_channel_initial_hertzian'].to_numpy()[ix]
        hertzian_energy_annulus_err=wdf['annulus_channel_initial_hertzian_std'].to_numpy()[ix]
        hertzian_energy_straight=wdf['straight_channel_initial_hertzian'].to_numpy()[ix]
        hertzian_energy_straight_err=wdf['straight_channel_initial_hertzian_std'].to_numpy()[ix]

        xs=1/rs**2

        ax.plot(
            xs,hertzian_energy_annulus/(2*np.pi*ww*rs),'-',
            marker='o',
            mfc='None',mec='k',label=rf"$w={ww}$ initial annulus"
            )
        ax.plot(
            xs,hertzian_energy_straight/(2*np.pi*ww*rs),'-',
            marker='d',
            mfc='None',mec='k',label=rf"$w={ww}$ initial straight"
            )

        hertzian_energy_annulus=wdf['annulus_channel_relaxed_hertzian'].to_numpy()[ix]
        hertzian_energy_annulus_err=wdf['annulus_channel_relaxed_hertzian_std'].to_numpy()[ix]
        hertzian_energy_straight=wdf['straight_channel_relaxed_hertzian'].to_numpy()[ix]
        hertzian_energy_straight_err=wdf['straight_channel_relaxed_hertzian_std'].to_numpy()[ix]
        ax.plot(
            xs,hertzian_energy_annulus/(2*np.pi*ww*rs),'--',
            marker='o',
            mfc='None',mec='k',label=rf"$w={ww}$ relax annulus"
            )
        ax.plot(
            xs,hertzian_energy_straight/(2*np.pi*ww*rs),'--',
            marker='d',
            mfc='None',mec='k',label=rf"$w={ww}$ relax straight"
            )
        ax.axvline(1e-4*3*4/(3.75**2),linestyle=':')
        # m,c,m_err=getLinearFit(xs,energies_initial,stds_initial)
        # fit_params.append([ww,m,m_err])
        # kf=getSmoothedGradient(rs,energies_initial)
        # kfs.append([ww,kf.mean()])
    ax.set_xlabel("$1/R^2$")
    ax.set_ylabel("$E/2\pi w R$")
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.legend(ncol=2)
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/width_sweep{mod}.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

def getPackingFractionFit(bdir):
    data=[]
    channels=['annulus','straight']
    for dd in tqdm(glob(f"{bdir}/*")):
        try:
            # energies=np.loadtxt(f"{dd}/energies.dat",delimiter='\t',skiprows=1)
            hps=np.loadtxt(f"{dd}/LcDist.log",delimiter='\t',skiprows=1)
            # row=[*hps,*energies.mean(axis=0),*energies.std(axis=0)]
            pfs=[]
            for channel in channels:
                RodShapedBacterium.sus_vis_radius_factor=1
                cell_file=f"{dd}/biofilm_{channel}_channel_initial_00000.dat"
                cells=pd.read_csv(cell_file,delimiter='\t')
                cell_area=(2*cells['length']*cells['radius']+np.pi*cells['radius']**2).sum()
                # packing_fraction=getPackingFraction(cells)
                packing_fraction=cell_area/(2*np.pi*hps[-2]*hps[-1])
                pfs.append(packing_fraction)
            data.append([*hps,*pfs])
        except Exception as e:
            print(e)
            quit()
    log_names=list(pd.read_csv(f"{dd}/LcDist.log",sep='\t',nrows=1).columns)
    channel_names=[f"{channel}_pf" for channel in channels]
    df=pd.DataFrame(data,columns=log_names+channel_names)
    df.sort_values(by=['annulus_width','annulus_radius'],inplace=True)
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    widths=df['annulus_width'].unique()
    for ww in widths:
        wdf=df[df['annulus_width']==ww]
        rs=wdf['annulus_radius'].to_numpy()
        pfs=wdf['annulus_pf'].to_numpy()
        ax.plot(rs,pfs,label=rf"$w={ww:1.2f}$")
    ax.set_xlabel("$R$")
    ax.set_ylabel(r"$\rho$")
    ax.legend()
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/packing_fraction_vs_width_vs_R.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)

    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    data=[]
    for ww in widths:
        wdf=df[df['annulus_width']==ww]
        pf=wdf['annulus_pf'].mean()
        data.append([ww,pf])
    data=np.array(data)
    xs=data[:,0]
    ys=data[:,1]
    f=lambda x,a: a/x
    popt,pcov=scopt.curve_fit(
        f,xs,ys,
        p0=[1]
        )
    m=popt
    print(f"{m=}")
    ax.plot(xs,ys,label=r"$\rho$")
    ax.plot(xs,f(xs,*popt),'--',label=rf"$\rho={m[0]:2.2f}/w$")
    ax.set_xlabel("$w$")
    ax.set_ylabel(r"$\rho$")
    ax.legend()
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/packing_fraction_vs_width.pdf",
        bbox_inches='tight',transparent=True
        )
    plt.show(block=False)
    return m

def plotProfileAndGetFitParams(fdf,qty,lbl,ax,markervery=1):
    if len(fdf['annulus_width'].unique())!=1 or len(fdf['LinkingProb'].unique())!=1:
        print(f"{len(fdf['annulus_width'].unique())=}")
        print(f"{len(fdf['LinkingProb'].unique())=}")
        print("Error! There should only be one width and linking prob")
        quit()
    ww=fdf['annulus_width'].unique()[0]
    lp=fdf['LinkingProb'].unique()[0]
    rs=fdf['annulus_radius'].to_numpy()
    ix=np.argsort(rs)[::-1]
    energies=fdf[f'{qty}'].to_numpy()[ix]
    stds=fdf[f'{qty}_std'].to_numpy()[ix]

    # error_label=rf"""$w={ww:1.1f}$
    # $p_\ell={lp:1.1f}$"""
    error_label=rf"""$\rho={10.48/ww:1.2f}$"""

    rs=rs[ix]
    xs=1/rs**2
    imin=np.argmin(energies)
    # m,c,m_err=getLinearFit(xs[imin:],energies[imin:],stds[imin:])
    m,c,m_err,imax=getFitToLinearPart(xs,energies,stds)
    ebar=ax.errorbar(
            xs,energies,stds,
            marker='o',capsize=5,
            mfc='None',mec='k',
            markevery=markervery,errorevery=markervery,
            label=error_label,
            zorder=1,
            alpha=0.7
            )
    colour=ebar[0].get_color()
    ax.plot(xs[imin:imax],xs[imin:imax]*m+c,color='k',
            alpha=1,zorder=3,lw=1.5)
    # ax.axvline(xs[imax],linestyle='-.',color=colour,zorder=4,lw=1.3,alpha=0.6)
    fit_params=[ww,lp,m,m_err]
    return fit_params

def getRelaxedData(fdf,ww,lp,add_relax=True):

    print(fdf[['LinkingProb','annulus_width','normalised_initial_delta_E','normalised_relaxed_delta_E']])

    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())

    fit_params=plotProfileAndGetFitParams(
        fdf,'normalised_initial_delta_E','in',ax
        )
    if add_relax:
        fit_params_relaxed=plotProfileAndGetFitParams(
            fdf,'normalised_relaxed_delta_E','re',ax
            )
    else:
        fit_params_relaxed=None

    xc=1e-4*3*4/(3.75**2)
    ax.axvline(xc,linestyle=':',zorder=4,lw=1.3,
               label=r"$\frac{1}{R^{2}}=\frac{12\phi_c^2}{\langle\ell \rangle^2}$")
    ax.set_xlabel("$R^{-2}$")
    ax.set_ylabel(r"$\varepsilon_H$")
    ax.legend(ncol=1,loc="upper left")
    fig.savefig(
        (f"../GeneratedOutput/AnalysisResults/{data_dir}"+
        f"/width_sweep_relaxed_fit_partial_r2_{lp=}_{ww=}.pdf"),
        bbox_inches='tight',
        transparent=True
        )
    # plt.show()
    plt.close()
    return fit_params,fit_params_relaxed

def getRelaxedK(fit_params,fit_params_relaxed):
    if np.any(fit_params_relaxed==None):
        relax=False
    else:
        relax=True
    wdths=np.unique(fit_params[:,0])
    lps=np.unique(fit_params[:,1])
    print(wdths)
    print(lps)
    K_dependence=[]
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    for lp in lps:
        w_fp=fit_params[fit_params[:,1]==lp]

        xtmp=10.5/w_fp[:,0]
        fit=lambda x,k,a: k*(x-a)
        popt,pcov=scopt.curve_fit(
            fit,xtmp,w_fp[:,2],sigma=w_fp[:,3],
            p0=[1,1])
        m,rho_c=popt
        K_dependence.append([lp,m,np.sqrt(np.diag(pcov))[0]])
        ax.plot(xtmp,m*(xtmp-rho_c),'k--')
        ax.errorbar(
            10.5/w_fp[:,0],w_fp[:,2],yerr=w_fp[:,3],
            linestyle='-',marker='o',mec='k',mfc='None',capsize=5,
            label=rf'$p_\ell={lp}$',alpha=0.8
            )
        if relax:
            w_fpr=fit_params_relaxed[fit_params_relaxed[:,1]==lp]
            ax.errorbar(
                10.5/w_fpr[:,0],w_fpr[:,2],yerr=w_fpr[:,3],
                linestyle='--',marker='d',mec='k',mfc='None',capsize=5,
                label=rf'relaxed {lp=}'
                )
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$K$")
    # ax.set_yscale('log')
    ax.legend(ncol=1)
    # ax.set_ylim([None,11.5])

    if relax==False:
        fig.savefig(
            f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_w_scal_mathcalK_initial.pdf",
            bbox_inches='tight',transparent=True
            )
    else:
        fig.savefig(
            f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_w_scal_mathcalK_relaxed.pdf",
            bbox_inches='tight',transparent=True
            )
    plt.show(block=False)


    np.savetxt(f"../GeneratedOutput/AnalysisResults/{data_dir}/K0_dependence.txt",K_dependence)
    K_dependence=np.array(K_dependence)
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    ax.errorbar(K_dependence[:,0],K_dependence[:,1],K_dependence[:,2],
                marker='o',mec='k',mfc='None',capsize=5,alpha=0.8,
                label=r'$K_{F0}(p_\ell)$')
    m,c,err=getLinearFit(K_dependence[:,0],K_dependence[:,1],K_dependence[:,2])
    ax.plot(K_dependence[:,0],m*K_dependence[:,0]+c,'k--',
            label=rf'${m:1.2f}p_\ell+{c:1.2f}$',
            zorder=3)
    ax.set_xlabel('$p_\ell$')
    ax.set_ylabel('$K$')
    ax.legend()
    fig.savefig(
        f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_p_ell.pdf",
        bbox_inches='tight',transparent=True,format='pdf'
        )
    plt.show()
    print(f"{m=} {c=}")

    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize())
    for ww in wdths:
        w_fp=fit_params[fit_params[:,0]==ww]
        ax.errorbar(
            w_fp[:,1],w_fp[:,2],yerr=w_fp[:,3],
            linestyle='-',marker='o',mec='k',mfc='None',capsize=5,
            label=rf'$\rho={10.48/ww:1.2f}$',alpha=0.8
            )
        if relax:
            w_fpr=fit_params_relaxed[fit_params_relaxed[:,0]==ww]
            ax.errorbar(
                w_fpr[:,1],w_fpr[:,2],yerr=w_fpr[:,3],
                linestyle='--',marker='d',mec='k',mfc='None',capsize=5,
                label=rf'relaxed {ww=}'
                )
    ax.set_xlabel(r"$p_\ell$")
    # ax.set_yscale('log')
    ax.set_ylabel(r"$K$")
    ax.legend(ncol=2)
    # ax.set_ylim([None,11.5])

    if relax:
        fig.savefig(
            f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_p_ell_scal_mathcalK_relaxed.pdf",
            bbox_inches='tight',transparent=True
            )
    else:
        fig.savefig(
            f"../GeneratedOutput/AnalysisResults/{data_dir}/kappa_vs_p_ell_scal_mathcalK_initial.pdf",
            bbox_inches='tight',transparent=True
            )
    plt.show(block=False)

def getDependenciesRelaxed(bdir,add_relax=True):

    df=getDF(bdir)
    df.sort_values(by=['annulus_width','annulus_radius'],inplace=True)
    mod=f'_scal_mathcalK_relaxed'
    widths=df['annulus_width'].unique()
    lps=df['LinkingProb'].unique()
    lps.sort()
    print(widths)
    print(lps)

    fit_params=[]
    fit_params_relaxed=[]
    for ww in widths:
        for lp in lps:
            fdf=df[(df['LinkingProb']==lp)&(df['annulus_width']==ww)]
            print(fdf[['LinkingProb','annulus_width',
                       'normalised_initial_delta_E',
                       'normalised_relaxed_delta_E']])
            if len(fdf)==0:
                continue
            tmp_in,tmp_re=getRelaxedData(
                fdf,ww,lp,add_relax=add_relax
                )
            fit_params.append(tmp_in)
            fit_params_relaxed.append(tmp_re)
    fit_params=np.array(fit_params)
    fit_params_relaxed=np.array(fit_params_relaxed)
    getRelaxedK(fit_params,fit_params_relaxed)

    lp=lps[0]
    fig,ax=plt.subplots(1,1,figsize=ut.getFigsize(2))
    for ii,ww in enumerate(widths):
        if ww>8:
            continue
        fdf=df[(df['LinkingProb']==lp)&(df['annulus_width']==ww)]
        plotProfileAndGetFitParams(
            fdf,'normalised_initial_delta_E','in',ax,markervery=(ii%3,3)
            )
    xc=1e-4*3*4/(3.75**2)
    ax.axvline(xc,linestyle='-.',color='m',alpha=0.6,zorder=4,lw=1.5)
    ax.set_xlabel("$R^{-2}$")
    ax.set_ylabel(r"$\varepsilon_H$")
    ax.legend(ncol=4)
    fig.savefig(
        (f"../GeneratedOutput/AnalysisResults/{data_dir}"+
        f"/width_sweep_initial_fit_partial_r2.pdf"),
        bbox_inches='tight',
        transparent=True
        )
    plt.show(block=True)

if __name__ == "__main__":
    print("----------------------------------------------")

    data_dir="SCOWF_Perturbation"       # Standard directory
    data_dir="SCOWF_Perturbation_K0"    # Zoom in near 0

    bdir=f"/media/rory/Elements/ScalingOrientation/{data_dir}"
    # data_dir="ScalingOrientationWidthRelaxPackingFraction"
    # bdir=f"../GeneratedOutput/SimOutput/{data_dir}"
    # plotEnergies(bdir)
    # getDependenceWidth(bdir,energy_type='parallel')
    # getPackingFractionFit(bdir)
    # getDependenceLinkingProb(bdir)
    # getDependenciesRelaxed(bdir,add_relax=False)
    # quit()
    #
    # print("----------------------------------------------")
    # data_dir="ScalingOrientationWidthRelaxPackingFraction"
    # bdir=f"../GeneratedOutput/SimOutput/{data_dir}"
    # getDependenceWidth(bdir,energy_type='parallel')
    # quit()

    data_dir="ScalingOrientationWidthNoRelax"
    bdir=f"/media/rory/Elements/ScalingOrientation/{data_dir}"
    getDependenceLinkingProb(bdir)
    # getPackingFractionFit(bdir)

    data_dir="ScalingOrientationNoRelax"
    bdir=f"/media/rory/Elements/ScalingOrientation/{data_dir}"
    # getDependenceLinkingProb(bdir)

    data_dir="ScalingOrientationWidthNoRelax"
    bdir=f"/media/rory/Elements/ScalingOrientation/{data_dir}"
    # getDependenceWidth(bdir)
    # plotEnergies(bdir)
    quit()
