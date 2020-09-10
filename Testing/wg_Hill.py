from pysb import *

Model()

Parameter('H_wg', 1.0 / 65.74263)
Parameter('K_WGwg', 0.010839762)
Parameter('nu_WGwg', 5) #5.2690034)
Parameter('alpha_WGwg', 1.2893039)
Parameter('K_CIwg', 0.1358894)
Parameter('nu_CIwg', 3) #3.0448556)
Parameter('alpha_CIwg', 5.354691)
Parameter('K_CNwg', 0.08410451)
Parameter('nu_CNwg', 2) #2.2106833)

Monomer('wg') # wingless mRNA
Monomer('IWG') # intracellular wingless protein
Monomer('CI') # cubitus interruptus protein
Monomer('CN') # cleaved cubitus interruptus protein

Initial(wg(), Parameter('wg_0'))
Initial(IWG(), Parameter('IWG_0'))
Initial(CI(), Parameter('CI_0'))
Initial(CN(), Parameter('CN_0'))

Observable('wg_obs', wg())
Observable('CN_obs', CN())
Observable('CI_obs', CI())
Observable('IWG_obs', IWG())

Expression('Hill_CNwg', (K_CNwg ** nu_CNwg) / (K_CNwg ** nu_CNwg + CN_obs ** nu_CNwg))

Expression('Hill_CIwg',
           alpha_CIwg * (CI_obs * Hill_CNwg) ** nu_CIwg / (K_CIwg ** nu_CIwg + (CI_obs * Hill_CNwg) ** nu_CIwg))

Expression('Hill_IWGwg',
           alpha_WGwg * (IWG_obs ** nu_WGwg) / (K_WGwg ** nu_WGwg + IWG_obs ** nu_WGwg))

Expression('Hill_CI_IWG_combined',
           H_wg * (Hill_CIwg + Hill_IWGwg) / (1 + Hill_CIwg + Hill_IWGwg))

Rule('synthesize_wg_1', None >> wg(), Hill_CI_IWG_combined)
Rule('degrade_wg_1', wg() >> None, Parameter('deg',H_wg.value))

#print (model.rules)