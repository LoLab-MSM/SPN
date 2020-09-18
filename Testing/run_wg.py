import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
import wg
from wg_Hill import model as model2
from itertools import product
import sys
import re

nu_WGwg = [5]  # [5.2690034]
nu_CIwg = [3]  # [3.0448556]
nu_CNwg = [2]  # [2.2106833]

for nu_wg, nu_ci, nu_cn in product(nu_WGwg, nu_CIwg, nu_CNwg):
    #print (nu_wg, nu_ci, nu_cn)

    wg.create_model(nu_wg, nu_ci, nu_cn)
    model1 = wg.model

    model2.parameters['nu_WGwg'].value = model1.parameters['nu_WGwg'].value
    model2.parameters['nu_CIwg'].value = model1.parameters['nu_CIwg'].value
    model2.parameters['nu_CNwg'].value = model1.parameters['nu_CNwg'].value

    H_wg = model1.parameters['H_wg']
    H_IWG = model1.parameters['H_IWG']
    H_EWG = model1.parameters['H_EWG']
    H_ci = model1.parameters['H_ci']
    H_CI = model1.parameters['H_CI']
    H_CN = model1.parameters['H_CN']

    tspan = np.linspace(0, 1000, 501)
    sim = ScipyOdeSimulator(model1, tspan, verbose=False)
    sim2 = ScipyOdeSimulator(model2, tspan, verbose=False, integrator_options={'atol': 1e-6, 'rtol': 1e-6})

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

    n_samples = 30
    CI_init = np.linspace(0.01, 1., n_samples)
    CN_init = np.linspace(0.1, .2, n_samples)
    IWG_init = np.linspace(0.01, 1, n_samples)

    mWG_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    mWG_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    mWG_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))
    pIWG_final_CI = np.zeros((len(IWG_init), len(CI_init)))
    pIWG_final_CN = np.zeros((len(IWG_init), len(CN_init)))
    pCN_final_CI = np.zeros((len(CN_init), len(CI_init)))
    pCN_final_IWG = np.zeros((len(CN_init), len(IWG_init)))
    pCI_final_CN = np.zeros((len(CI_init), len(CN_init)))
    pCI_final_IWG = np.zeros((len(CI_init), len(IWG_init)))

    wg_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    wg_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    wg_final_CN_IWG = np.zeros((len(CN_init), len(IWG_init)))

    for i, w in enumerate(IWG_init):
        for j, c in enumerate(CI_init):
            for k, cl in enumerate(CN_init):
                print('(%s, %g, %g): %d.%d.%d' % (nu_wg, nu_cn, nu_ci, i, j, k))
                # Simulate full model
                x = sim.run(initials={model1.initial_conditions[2][0]: w / (H_IWG.value * H_wg.value), \
                                      model1.initial_conditions[0][0]: c / (H_CI.value * H_ci.value), \
                                      model1.initial_conditions[1][0]: cl / (H_CN.value * H_ci.value)})

                mWG = x.observables['mWng_obs']
                scaled_mWG = x.observables['mWng_obs'] * H_wg.value
                #plt.plot(tspan,scaled_mWG,color='red')
                pIWG = x.observables['pIWG_free']
                scaled_pIWG = pIWG * (H_IWG.value * H_wg.value)  # H's are already divisors (see en.py)
                pCN = x.observables['pCN_free']
                scaled_pCN = pCN * (H_CN.value * H_ci.value)  # also already divisors
                pCI = x.observables['pCI_free']
                scaled_pCI = pCI * (H_CI.value * H_ci.value)  # also already divisors

                # exactly which one of these to use only matters for plotting
                mWG_final_CI_CN[j][k] = scaled_mWG[-1]
                mWG_final_CI_IWG[j][i] = scaled_mWG[-1]
                mWG_final_CN_IWG[k][i] = scaled_mWG[-1]
                pIWG_final_CI[i][j] = scaled_pIWG[-1]
                pIWG_final_CI[i][k] = scaled_pIWG[-1]
                pCN_final_CI[k][j] = scaled_pCN[-1]
                pCN_final_IWG[k][i] = scaled_pCN[-1]
                pCI_final_CN[j][k] = scaled_pCI[-1]
                pCI_final_IWG[j][i] = scaled_pCI[-1]

                # Simulate reduced model
                if not np.isnan(scaled_pIWG[-1]):
                    x = sim2.run(initials={model2.initial_conditions[1][0]: scaled_pIWG[-1],
                                           model2.initial_conditions[2][0]: scaled_pCI[-1],
                                           model2.initial_conditions[3][0]: scaled_pCN[-1]})
                    wg_final_CI_CN[j][k] = x.observables['wg_obs'][-1]
                    wg_final_CI_IWG[j][i] = x.observables['wg_obs'][-1]
                    wg_final_CN_IWG[k][i] = x.observables['wg_obs'][-1]
                    #plt.plot(tspan,x.observables['wg_obs'],color='green',linestyle='--')
                    #plt.show()
                else:
                    wg_final_CI_CN[j][k] = float('nan')
                    wg_final_CI_IWG[j][i] = float('nan')
                    wg_final_CN_IWG[k][i] = float('nan')
    #plt.show()

    plt.figure()
    # CI, CN
    cm = plt.pcolormesh(CI_init, CN_init, mWG_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CI]$_0$ (scaled conc)')
    plt.ylabel(r'[CN]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    plt.figure()
    cm = plt.pcolormesh(CI_init, CN_init, wg_final_CI_CN, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CI]$_0$ (scaled conc)')
    plt.ylabel(r'[CN]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
    # CI, IWG
    plt.figure()
    cm = plt.pcolormesh(CI_init, IWG_init, mWG_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CI]$_0$ (scaled conc)')
    plt.ylabel(r'[IWG]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    plt.figure()
    cm = plt.pcolormesh(CI_init, IWG_init, wg_final_CI_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CI]$_0$ (scaled conc)')
    plt.ylabel(r'[IWG]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
    #CN, IWG
    plt.figure()
    cm = plt.pcolormesh(CN_init, IWG_init, mWG_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CN]$_0$ (scaled conc)')
    plt.ylabel(r'[IWG]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[mWG]$_\infty$ (scaled conc)', labelpad=20)
    #
    plt.figure()
    cm = plt.pcolormesh(CN_init, IWG_init, wg_final_CN_IWG, vmin=0, vmax=1, cmap='bwr')
    cm.get_cmap().set_bad('k')
    plt.xlabel(r'[CN]$_0$ (scaled conc)')
    plt.ylabel(r'[IWG]$_0$ (scaled conc)')
    plt.xscale('log')
    plt.yscale('log')
    cb = plt.colorbar()
    cb.ax.set_ylabel(r'[wg]$_\infty$ (scaled conc)', labelpad=20)
plt.show()







