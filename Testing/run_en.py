import numpy as np
import matplotlib.pyplot as plt
from pysb.simulator import ScipyOdeSimulator
import en
from en_Hill import model as model2
from itertools import product
import sys

nu_WGen = [8] #[8.195488]
nu_CNen = [6] #[5.929026]

for nu_wg, nu_cn in product(nu_WGen, nu_CNen):
#     print nu_wg, nu_cn
    
    en.create_model(nu_wg, nu_cn)
    model1 = en.model


    model2.parameters['nu_WGen'].value = model1.parameters['nu_WGen'].value
    model2.parameters['nu_CNen'].value = model1.parameters['nu_CNen'].value
    
    H_en = model1.parameters['H_en']
    ktr_en = model1.parameters['ktr_en']
    H_wg = model1.parameters['H_wg']
    H_WG = model1.parameters['H_WG']
    H_EWG = model1.parameters['H_EWG']
    H_ci = model1.parameters['H_ci']
    H_CI = model1.parameters['H_CI']
    H_CN = model1.parameters['H_CN']

    tspan = np.linspace(0, 1000, 501)
    sim = ScipyOdeSimulator(model1, tspan, verbose=False)
    sim2 = ScipyOdeSimulator(model2, tspan, verbose=False, integrator_options = {'atol' : 1e-6, 'rtol' : 1e-6})

    n_samples = 10 #100
    #EWG_init = np.linspace(1./n_samples, 1, n_samples) #1
    #EWG_init = np.linspace(0,1,11) #[1]
    EWG_init = np.linspace(0.01, 1., n_samples)
    #CN_init = [6] #np.linspace(0,6,101)
    #CN_init = np.linspace(6./n_samples, 6, n_samples) #6
    #CN_init = np.linspace(0,1,11) #np.linspace(6./n_samples, 6, n_samples) #6
    CN_init = np.linspace(0.001, .01, n_samples)
    
    mEN_final = np.zeros((len(EWG_init), len(CN_init)))
    pEWG_final = np.zeros((len(EWG_init), len(CN_init)))
    pCN_final = np.zeros((len(EWG_init), len(CN_init)))
    
    en_final = np.zeros((len(EWG_init), len(CN_init)))

    for i,e in enumerate(EWG_init):
        for j,c in enumerate(CN_init):
            print('(%s, %g): %d.%d' % (nu_wg, nu_cn, i, j))
            # print (e)
            # print (c)
            # print (str(float((e/(H_EWG.value*H_wg.value))/(c/(H_CN.value*H_ci.value)))))
            # this ratio (directly above) is about double the e to c ratio
            # Simulate full model
            x = sim.run(initials={model1.initial_conditions[0][0] : e / (H_EWG.value*H_wg.value), \
                                  model1.initial_conditions[1][0] : c / (H_CN.value*H_ci.value)})
            mEN = x.observables['mEN_obs']
            scaled_mEN = x.observables['mEN_obs'] * H_en.value / ktr_en.value
            pEWG = x.observables['pEWG_free']
            scaled_pEWG = pEWG * (H_EWG.value*H_wg.value) # H's are already divisors (see en.py)
#             scaled_pEWG = (x.observables['pEWG_free'] + x.observables['pEWG_pCN'])*H_en.value
            pCN = x.observables['pCN_free']
            scaled_pCN = pCN * (H_CN.value*H_ci.value) # also already divisors
            mEN_final[i][j] = scaled_mEN[-1]
            pEWG_final[i][j] = scaled_pEWG[-1]
            pCN_final[i][j] = scaled_pCN[-1]
            # Simulate reduced model
            if not np.isnan(scaled_pEWG[-1]):
                x = sim2.run(initials={model2.initial_conditions[1][0] : scaled_pEWG[-1],
                                       model2.initial_conditions[2][0] : scaled_pCN[-1]})
                en_final[i][j] = x.observables['en_obs'][-1]
            else:
                en_final[i][j] = float('nan')

    plt.figure()
    if len(CN_init) == 1:
        plt.plot(pEWG_final, mEN_final, color='0.5', lw=2)
#         plt.plot(EWG_init, mEN_final, color='0.5', lw=2)
        idx = np.arange(0,len(EWG_init),int(len(EWG_init)/100))
        plt.plot(pEWG_final[idx], en_final[idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg,nu_cn))
#         plt.plot(EWG_init[idx], en_final[idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg,nu_cn))
        plt.xlabel(r'[pEWG]$_\infty$ (scaled conc.)')
        plt.ylabel(r'[mEN]$_\infty$ (scaled conc.)')
        plt.xscale('log')
        plt.legend(loc=0)
    elif len(EWG_init) == 1:
        plt.plot(pCN_final[0], mEN_final[0], color='0.5', lw=2)
#         plt.plot(CN_init, mEN_final[0], color='0.5', lw=2)
        idx = np.arange(0,len(CN_init),int(len(CN_init)/100))
        plt.plot(pCN_final[0][idx], en_final[0][idx], 'o', label='nu_WGen = %g\nnu_CNen = %g' % (nu_wg,nu_cn))
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







