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
    specdict = {}
    for n,i in enumerate(model1.species):
        specdict[n] = i
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

    # NEED TO FIGURE OUT THE RIGHT BOUNDS ETC SO CAN GET GOOD NUMBERS
    n_samples = 5 #10 #100
    CI_init = np.linspace(0.01, 1., n_samples)
    CN_init = np.linspace(0.001, .01, n_samples)
    IWG_init = np.linspace(0.01, 1, n_samples)

    mWG_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    mWG_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))
    pIWG_final_CI = np.zeros((len(IWG_init), len(CI_init)))
    pIWG_final_CN = np.zeros((len(IWG_init), len(CN_init)))
    pCN_final_CI = np.zeros((len(CN_init), len(CI_init)))
    pCN_final_IWG = np.zeros((len(CN_init), len(IWG_init)))
    pCI_final_CN = np.zeros((len(CI_init), len(CN_init)))
    pCI_final_IWG = np.zeros((len(CI_init), len(IWG_init)))

    wg_final_CI_CN = np.zeros((len(CI_init), len(CN_init)))
    wg_final_CI_IWG = np.zeros((len(CI_init), len(IWG_init)))

    for i, w in enumerate(IWG_init):
        for j, c in enumerate(CI_init):
            for k, cl in enumerate(CN_init):
                print('(%s, %g, %g): %d.%d.%d' % (nu_wg, nu_cn, nu_ci, i, j, k))
                #c = 0
                #w = 0
                # print (e)
                # print (c)
                # print (str(float((e/(H_EWG.value*H_wg.value))/(c/(H_CN.value*H_ci.value)))))
                # this ratio (directly above) is about double the e to c ratio
                # Simulate full model
                x = sim.run(initials={model1.initial_conditions[2][0]: w / (H_IWG.value * H_wg.value), \
                                      model1.initial_conditions[0][0]: c / (H_CI.value * H_ci.value), \
                                      model1.initial_conditions[1][0]: cl / (H_CN.value * H_ci.value)})

                mWG = x.observables['mWng_obs']
                scaled_mWG = x.observables['mWng_obs'] * H_wg.value
                plt.plot(tspan,scaled_mWG,color='red')
                #plt.plot(tspan[-3:], x.observables['gWng_boundby_A'][-3:]+x.observables['gWng_boundby_both'][-3:]+x.observables['gWng_boundby_Y'][-3:], color='red')
                #plt.plot(tspan[-3:], x.observables['gWng_bound_total'][-3:], color='green', linestyle='--')
                #plt.plot(tspan[-3:], x.observables['gWng_boundby_A'][-3:]+x.observables['gWng_boundby_Y'][-3:], color='blue', linestyle=':')
                #plt.show()
                pIWG = x.observables['pIWG_free']
                scaled_pIWG = pIWG * (H_IWG.value * H_wg.value)  # H's are already divisors (see en.py)
                pCN = x.observables['pCN_free']
                scaled_pCN = pCN * (H_CN.value * H_ci.value)  # also already divisors
                pCI = x.observables['pCI_free']
                scaled_pCI = pCI * (H_CI.value * H_ci.value)  # also already divisors

                # exactly which one of these to use only matters for plotting
                mWG_final_CI_CN[j][k] = scaled_mWG[-1]
                mWG_final_CI_IWG[j][i] = scaled_mWG[-1]
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
                    plt.plot(tspan,x.observables['wg_obs'],color='green',linestyle='--')
#                    plt.title(str(w)+','+str(c)+','+str(cl))
                    #plt.show()
                else:
                    wg_final_CI_CN[j][k] = float('nan')
                    wg_final_CI_IWG[j][i] = float('nan')
#                break
#            break
#        break
    plt.show()

    # FIX THE BELOW LATER
    # even though it was never super important to plot them one-dimensionally...

    sys.exit(0)
    plt.figure()
    if len(CN_init) == 1:
        plt.plot(pEWG_final, mEN_final, color='0.5', lw=2)
        #         plt.plot(EWG_init, mEN_final, color='0.5', lw=2)
        idx = np.arange(0, len(EWG_init), int(len(EWG_init) / 100))
        plt.plot(pEWG_final[idx], en_final[idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg, nu_cn))
        #         plt.plot(EWG_init[idx], en_final[idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg,nu_cn))
        plt.xlabel(r'[pEWG]$_\infty$ (scaled conc.)')
        plt.ylabel(r'[mEN]$_\infty$ (scaled conc.)')
        plt.xscale('log')
        plt.legend(loc=0)
    elif len(EWG_init) == 1:
        plt.plot(pCN_final[0], mEN_final[0], color='0.5', lw=2)
        #         plt.plot(CN_init, mEN_final[0], color='0.5', lw=2)
        idx = np.arange(0, len(CN_init), int(len(CN_init) / 100))
        plt.plot(pCN_final[0][idx], en_final[0][idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg, nu_cn))
        #         plt.plot(CN_init[idx], en_final[0][idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg,nu_cn))
        plt.xlabel(r'[pCN]$_\infty$ (scaled conc.)')
        plt.ylabel(r'[mEN]$_\infty$ (scaled conc.)')
        plt.xscale('log')
        plt.legend(loc=0)
    else:
        cm = plt.pcolormesh(CN_init, EWG_init, mEN_final, vmin=0, vmax=1, cmap='bwr')
        cm.get_cmap().set_bad('k')
        plt.xlabel(r'[CN]$_0$ (scaled conc)')
        plt.ylabel(r'[EWG]$_0$ (scaled conc)')
        plt.xscale('log')
        plt.yscale('log')
        cb = plt.colorbar()
        cb.ax.set_ylabel(r'[mEN]$_\infty$ (scaled conc)', labelpad=20)
        #
        plt.figure()
        cm = plt.pcolormesh(CN_init, EWG_init, en_final, vmin=0, vmax=1, cmap='bwr')
        cm.get_cmap().set_bad('k')
        plt.xlabel(r'[CN]$_0$ (scaled conc)')
        plt.ylabel(r'[EWG]$_0$ (scaled conc)')
        plt.xscale('log')
        plt.yscale('log')
        cb = plt.colorbar()
        cb.ax.set_ylabel(r'[en]$_\infty$ (scaled conc)', labelpad=20)
plt.show()







