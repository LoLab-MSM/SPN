import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
#import wg
#import en
from full_model import model as model1
from full_Hill_ODEs import model as model2
from itertools import product
import sys
import re

nu_WGwg = [5]  # [5.2690034]
nu_CIwg = [3]  # [3.0448556]
nu_CNwg = [2]  # [2.2106833]
nu_WGen = [8] #[8.195488]
nu_CNen = [6] #[5.929026]
nu_CIptc= [3] # should it be 4?? #[3.5400295]
nu_CNptc= [9] #[9.436314]
nu_ENci= [3] #[3.1201758]

# since nu_cnptc > nu_cnen > nu_cnwg
# and (for now?) nu_ciwg == nu_ciptc
for nu_wgwg, nu_ci, nu_wgen, nu_cn in product(nu_WGwg, nu_CIwg, nu_WGen, nu_CNptc):
    model2.parameters['nu_WGwg'].value = model1.parameters['nu_WGwg'].value
    model2.parameters['nu_CIwg'].value = model1.parameters['nu_CIwg'].value
    model2.parameters['nu_CNwg'].value = model1.parameters['nu_CNwg'].value
    model2.parameters['nu_WGen'].value = model1.parameters['nu_WGen'].value
    model2.parameters['nu_CNen'].value = model1.parameters['nu_CNen'].value
    model2.parameters['nu_CIptc'].value = model1.parameters['nu_CIptc'].value
    model2.parameters['nu_CNptc'].value = model1.parameters['nu_CNptc'].value
    model2.parameters['nu_ENci'].value = model1.parameters['nu_ENci'].value

    H_wg = model1.parameters['H_wg']
    H_IWG = model1.parameters['H_IWG']
    H_EWG = model1.parameters['H_EWG']
    H_ci = model1.parameters['H_ci']
    H_CI = model1.parameters['H_CI']
    H_CN = model1.parameters['H_CN']
    H_en = model1.parameters['H_en']
    H_EN = model1.parameters['H_EN']
    H_ptc = model1.parameters['H_ptc']

    tspan = np.linspace(0, 1000, 501)
    sim = ScipyOdeSimulator(model1, tspan, compiler='cython', verbose=False)
    sim2 = ScipyOdeSimulator(model2, tspan, compiler='cython', verbose=False, integrator_options={'atol': 1e-6, 'rtol': 1e-6})

    #print ('MODEL 1: MA')
    #specdict = {}
    #for n,i in enumerate(model1.species):
    #    specdict[n] = i
    #print (specdict)
    # repldict = {}
    # for n,i in enumerate(model1.odes):
    #     repldict[n] = str(i)
    # #print (repldict)
    # # do everything backwards to make sure you don't replace something like the following:
    # #   '__s16*...' re.sub('__s'+str(j),specdict[j],repldict[key]) where j = 1, so it returns
    # #   'pCN(cn_l=None,cn_r=None,gX2=None)6' -> it replaced '__s1' because it found that first
    # for key, repl in sorted(list(repldict.items()), reverse=True):
    #     print ("Checking {0} in {1}".format(key, repl))
    #     for j, spec in sorted(list(specdict.items()), reverse=True):
    #         repldict[key] = re.sub('__s' + str(j), str(spec), repldict[key])
    # for i in repldict:
    #     print (repldict[i])
    #     print ('')
    #for i in model1.reactions:
    #    print (i)
    #sys.exit(0)

    n_samples = 10
    # CI_init = np.linspace(0.001, 1., n_samples)
    # CN_init = np.linspace(0.1, 1, n_samples)
    # IWG_init = np.linspace(0.001, 1, n_samples)
    # EWG_init = np.linspace(0.001, 1., n_samples)
    CI_init = np.linspace(0.001, 1, n_samples) # CI more? probably
#    CI_init = np.linspace(0.001, 0.005, n_samples) # CI more? probably
    CN_init = np.linspace(0.001, .01, n_samples)
    IWG_init = np.linspace(0.001, .1, n_samples) #IWG less?? or broader range?
#    IWG_init = np.linspace(0.01, .02, n_samples) #IWG less?? or broader range?
    EWG_init = np.linspace(0.1, 1, n_samples)
    #EWG_init = np.linspace(0.001, 1, n_samples)
    #
    mWG_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    mWG_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    mWG_final_CI_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mWG_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    mWG_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mEN_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    mEN_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    mEN_final_CI_EWG = np.zeros((len(CI_init), len(EWG_init)))
    mEN_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    mEN_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mPTC_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    mPTC_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    mPTC_final_CI_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mPTC_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    mPTC_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mCI_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    mCI_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))

    wg_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    wg_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    wg_final_CI_EWG = np.zeros((len(CI_init), len(EWG_init)))
    wg_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    wg_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    en_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    en_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    en_final_CI_EWG = np.zeros((len(CI_init), len(EWG_init)))
    en_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    en_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    ptc_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    ptc_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    ptc_final_CI_EWG = np.zeros((len(CI_init), len(EWG_init)))
    ptc_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    ptc_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    ci_final_CN_EWG = np.zeros((len(CN_init), len(EWG_init)))
    ci_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))

    for i, w in enumerate(IWG_init):
        for j, c in enumerate(CI_init):
            for k, cl in enumerate(CN_init):
                for l, e in enumerate(EWG_init):
                    print('(%s, %g, %g, %g): %d.%d.%d.%d' % (nu_wgwg, nu_ci, nu_wgen, nu_cn,i,j,k,l))
                    # Simulate full model
                    x = sim.run(initials={model1.initial_conditions[2][0]: w / (H_IWG.value * H_wg.value), \
                                          model1.initial_conditions[0][0]: c / (H_CI.value * H_ci.value), \
                                          model1.initial_conditions[1][0]: cl / (H_CN.value * H_ci.value), \
                                          #model1.initial_conditions[23][0]: e / (H_en.value * H_EN.value)})#,
                                          model1.initial_conditions[3][0]: e / (H_EWG.value * H_wg.value)})
                    #
                    mWG = x.observables['mWng_obs']
                    scaled_mWG = x.observables['mWng_obs'] * H_wg.value
                    mEN = x.observables['mEN_obs']
                    scaled_mEN = x.observables['mEN_obs'] * H_en.value
                    mPTC = x.observables['mPTC_obs']
                    scaled_mPTC = x.observables['mPTC_obs'] * H_ptc.value
                    mCI = x.observables['mCI_obs']
                    scaled_mCI = x.observables['mCI_obs'] * H_ci.value
                    plt.plot(tspan,scaled_mCI,color='red')
                    #
                    pEWG = x.observables['pEWG_free']
                    scaled_pEWG = pEWG * (H_EWG.value * H_wg.value)  # H's are already divisors (see en.py)
                    pCN = x.observables['pCN_free']
                    scaled_pCN = pCN * (H_CN.value * H_ci.value)  # also already divisors
                    pIWG = x.observables['pIWG_free']
                    scaled_pIWG = pIWG * (H_IWG.value * H_wg.value)  # H's are already divisors (see en.py)
                    pCI = x.observables['pCI_free']
                    scaled_pCI = pCI * (H_CI.value * H_ci.value)  # also already divisors
                    pEN = x.observables['pEN_free']
                    scaled_pEN = pEN * (H_EN.value * H_en.value)
                    plt.plot(tspan,scaled_pEN,color='green')
                    pB = x.observables['B_free']
                    scaled_pB = pB * .4
                    #plt.plot(tspan,scaled_pB,color='blue')
                    #exactly which of these to use only matters for plotting
                    mWG_final_CI_CN[j][k] = scaled_mWG[-1]
                    mWG_final_CI_IWG[j][i] = scaled_mWG[-1]
                    mWG_final_CN_IWG[k][i] = scaled_mWG[-1]
                    mEN_final_CI_CN[j][k] = scaled_mEN[-1]
                    mEN_final_CI_IWG[j][i] = scaled_mEN[-1]
                    mEN_final_CN_IWG[k][i] = scaled_mEN[-1]
                    mPTC_final_CI_CN[j][k] = scaled_mPTC[-1]
                    mPTC_final_CI_IWG[j][i] = scaled_mPTC[-1]
                    mPTC_final_CN_IWG[k][i] = scaled_mPTC[-1]
                    mCI_final_CN_IWG[k][i] = scaled_mCI[-1]
                    mCI_final_CN_EWG[k][l] = scaled_mCI[-1]
                    # Simulate reduced model
                    if not np.isnan(scaled_pIWG[-1]):
                        param_values = np.array([p.value for p in model2.parameters])
                        param_values[25] = scaled_pB[-1]
                        x = sim2.run(initials={model2.initial_conditions[5][0]: scaled_pIWG[-1],
                                               model2.initial_conditions[6][0]: scaled_pCI[-1],
                                               model2.initial_conditions[1][0]: scaled_pEWG[-1],
                                               model2.initial_conditions[2][0]: scaled_pCN[-1],#})
                                               model2.initial_conditions[3][0]: scaled_pEN[-1]},
                                     param_values=param_values)
                        wg_final_CI_CN[j][k] = x.observables['wg_obs'][-1]
                        wg_final_CI_IWG[j][i] = x.observables['wg_obs'][-1]
                        wg_final_CN_IWG[k][i] = x.observables['wg_obs'][-1]
                        en_final_CI_CN[j][k] = x.observables['en_obs'][-1]
                        en_final_CI_IWG[j][i] = x.observables['en_obs'][-1]
                        en_final_CN_IWG[k][i] = x.observables['en_obs'][-1]
                        ptc_final_CI_CN[j][k] = x.observables['ptc_obs'][-1]
                        ptc_final_CI_IWG[j][i] = x.observables['ptc_obs'][-1]
                        ptc_final_CN_IWG[k][i] = x.observables['ptc_obs'][-1]
                        ci_final_CN_IWG[k][i] = x.observables['ci_obs'][-1]
                        ci_final_CN_EWG[k][l] = x.observables['ci_obs'][-1]
                        #print('low EN '+str(scaled_pCN[-1])+' ('+str(x.observables['EN_obs'][-1])+')')
                        plt.plot(tspan, x.observables['EN_obs'], color='pink', linestyle='--')
                        #plt.plot(tspan, x.observables['EN_obs'], color='blue', linestyle='--')
                        #plt.hlines(x.param_values[0][25],tspan[0],tspan[-1],color='turquoise',linestyle='--')
                    else:
                        wg_final_CI_CN[j][k] = float('nan')
                        wg_final_CI_IWG[j][i] = float('nan')
                        wg_final_CN_IWG[k][i] = float('nan')
                        en_final_CI_CN[j][k] = float('nan')
                        en_final_CI_IWG[j][i] = float('nan')
                        en_final_CN_IWG[k][i] = float('nan')
                        ptc_final_CI_CN[j][k] = float('nan')
                        ptc_final_CI_IWG[j][i] = float('nan')
                        ptc_final_CN_IWG[k][i] = float('nan')
                        ci_final_CN_IWG[k][i] = float('nan')
                        ci_final_CN_EWG[k][l] = float('nan')
    plt.show()
#    sys.exit(0)
    #
    # subplot arrangements
    #
    left = 0.075 # the value on the x axis where the left margin should be
    right = 0.95
    bottom = 0.09 # the value on the y axis where the bottom margin should be
    top = 0.95
    wspace = 0.4 # the proportion of the average of the left value and the right value -> for example, 0.2*((0.075+0.95)/2)
    hspace = 0.25

    # wg
    fig, ax = plt.subplots(3, 2,figsize=(12,12))
    fig.suptitle('wingless',fontsize=24)
    # CI, CN
    cm = ax[0,0].pcolormesh(CI_init, CN_init, mWG_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[0,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,0].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,0])
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[0,1].pcolormesh(CI_init, CN_init, wg_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[0,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,1].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,1])
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
    # CI, IWG
    cm = ax[1,0].pcolormesh(CI_init, IWG_init, mWG_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[1,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,0])
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[1,1].pcolormesh(CI_init, IWG_init, wg_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[1,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,1])
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
    #CN, IWG
    cm = ax[2,0].pcolormesh(CN_init, IWG_init, mWG_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[2,0].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,0].set_xscale('log')
    ax[2,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,0])
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[2,1].pcolormesh(CN_init, IWG_init, wg_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[2,1].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,1].set_xscale('log')
    ax[2,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,1])
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
    # en
    fig, ax = plt.subplots(3, 2,figsize=(12,12))
    fig.suptitle('engrailed',fontsize=24)
    # CI, CN
    cm = ax[0,0].pcolormesh(CI_init, CN_init, mEN_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[0,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,0].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,0])
    cb.ax.set_ylabel(r'[mEN]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[0,1].pcolormesh(CI_init, CN_init, en_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[0,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,1].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,1])
    cb.ax.set_ylabel(r'[en]$_\infty$ (scaled conc)', labelpad=20)
    # CI, IWG
    cm = ax[1,0].pcolormesh(CI_init, IWG_init, mEN_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[1,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,0])
    cb.ax.set_ylabel(r'[mEN]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[1,1].pcolormesh(CI_init, IWG_init, en_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[1,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,1])
    cb.ax.set_ylabel(r'[en]$_\infty$ (scaled conc)', labelpad=20)
    # CN, IWG
    cm = ax[2,0].pcolormesh(CN_init, IWG_init, mEN_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[2,0].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,0].set_xscale('log')
    ax[2,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,0])
    cb.ax.set_ylabel(r'[mEN]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[2,1].pcolormesh(CN_init, IWG_init, en_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[2,1].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,1].set_xscale('log')
    ax[2,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,1])
    cb.ax.set_ylabel(r'[en]$_\infty$ (scaled conc)', labelpad=20)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
    # ptc
    fig, ax = plt.subplots(3, 2,figsize=(12,12))
    fig.suptitle('patched',fontsize=24)
    # CI, CN
    cm = ax[0,0].pcolormesh(CI_init, CN_init, mPTC_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[0,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,0].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,0])
    cb.ax.set_ylabel(r'[mPTC]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[0,1].pcolormesh(CI_init, CN_init, ptc_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[0,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[0,1].set_ylabel(r'[CN]$_0$ (scaled conc)')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,1])
    cb.ax.set_ylabel(r'[ptc]$_\infty$ (scaled conc)', labelpad=20)
    # CI, IWG
    cm = ax[1,0].pcolormesh(CI_init, IWG_init, mPTC_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[1,0].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,0])
    cb.ax.set_ylabel(r'[mPTC]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[1,1].pcolormesh(CI_init, IWG_init, ptc_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[1,1].set_xlabel(r'[CI]$_0$ (scaled conc)')
    ax[1,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,1])
    cb.ax.set_ylabel(r'[ptc]$_\infty$ (scaled conc)', labelpad=20)
    # CN, IWG
    cm = ax[2,0].pcolormesh(CN_init, IWG_init, mPTC_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[2,0].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,0].set_xscale('log')
    ax[2,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,0])
    cb.ax.set_ylabel(r'[mPTC]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[2,1].pcolormesh(CN_init, IWG_init, ptc_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[2,1].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[2,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[2,1].set_xscale('log')
    ax[2,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[2,1])
    cb.ax.set_ylabel(r'[ptc]$_\infty$ (scaled conc)', labelpad=20)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
    # ci
    fig, ax = plt.subplots(2, 2,figsize=(12,12))
    fig.suptitle('cubitus interruptus',fontsize=24)
    # CN, EWG
    cm = ax[0,0].pcolormesh(CN_init, EWG_init, mCI_final_CN_EWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[0,0].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[0,0].set_ylabel(r'[EWG]$_0$ (scaled conc)')
    ax[0,0].set_xscale('log')
    ax[0,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,0])
    cb.ax.set_ylabel(r'[mCI]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[0,1].pcolormesh(CN_init, EWG_init, ci_final_CN_EWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[0,1].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[0,1].set_ylabel(r'[EWG]$_0$ (scaled conc)')
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[0,1])
    cb.ax.set_ylabel(r'[ci]$_\infty$ (scaled conc)', labelpad=20)
    # CN, IWG
    cm = ax[1,0].pcolormesh(CN_init, IWG_init, mCI_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    ax[1,0].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[1,0].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,0])
    cb.ax.set_ylabel(r'[mCI]$_\infty$ (scaled conc)', labelpad=20)
    #
    cm = ax[1,1].pcolormesh(CN_init, IWG_init, ci_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    #ax[1,1].set_xlabel(r'[CN]$_0$ (scaled conc)')
    ax[1,1].set_ylabel(r'[IWG]$_0$ (scaled conc)')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    cb = fig.colorbar(cm,ax=ax[1,1])
    cb.ax.set_ylabel(r'[ci]$_\infty$ (scaled conc)', labelpad=20)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
plt.show()







