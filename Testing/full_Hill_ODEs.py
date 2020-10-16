from pysb import *

Model()

# Engrailed
Parameter('H_en', 1.0 / 56.195484)
Parameter('K_WGen', 0.016870342)
Parameter('nu_WGen', 8)  # 2) #8.195488)
Parameter('K_CNen', 0.0015239986)
Parameter('nu_CNen', 6)  # 2) # 5.929026)
Parameter('H_EN', 1.0/26.655169)
# Wingless
Parameter('H_wg', 1.0 / 65.74263)
Parameter('K_WGwg', 0.010839762)
Parameter('nu_WGwg', 5) #5.2690034)
Parameter('alpha_WGwg', 1.2893039)
Parameter('K_CIwg', 0.1358894)
Parameter('nu_CIwg', 3) #3.0448556)
Parameter('alpha_CIwg', 5.354691)
Parameter('K_CNwg', 0.08410451)
Parameter('nu_CNwg', 2) #2.2106833)
# Patched
Parameter('K_CIptc', 0.0012549695)
Parameter('nu_CIptc', 3.) # should it be 4??
Parameter('H_ptc', 1.0/62.299126)
Parameter('K_CNptc', 0.15335482)
Parameter('nu_CNptc', 9)
# Cubitus Interruptus
Parameter('K_Bci', 0.016921775)
Parameter('nu_Bci', 1.0)
Parameter('H_ci', 1.0/88.05598)
Parameter('K_ENci', 0.05230557)
#Parameter('nu_ENci', 3.1201758)
Parameter('nu_ENci', 3)
Parameter('B', 0.4) # determines ci expression

# en
Monomer('en')  # engrailed mRNA
Monomer('EWG')  # extracellular wingless protein
Monomer('CN')  # cleaved cubitus interuptus protein
#
Monomer('EN')  # engrailed protein
# wg
Monomer('wg') # wingless mRNA
Monomer('IWG') # intracellular wingless protein
Monomer('CI') # cubitus interruptus protein
# already here from Engrailed
#Monomer('CN') # cleaved cubitus interruptus protein
# ptc
Monomer('ptc') # patched mRNA
# already here from engrailed, wingless
#Monomer('CN') # cleaved cubitus interruptus protein
#Monomer('CI') # cubitus interruptus protein
# ci
Monomer('ci')

# en
Initial(en(), Parameter('en_0'))
Initial(EWG(), Parameter('EWG_0'))
Initial(CN(), Parameter('CN_0'))
#
Initial(EN(), Parameter('EN_0'))
# wg
Initial(wg(), Parameter('wg_0'))
Initial(IWG(), Parameter('IWG_0'))
Initial(CI(), Parameter('CI_0'))
#Initial(CN(), Parameter('CN_0'))
# ptc
Initial(ptc(), Parameter('ptc_0'))
#Initial(CI(), Parameter('CI_0'))
#Initial(CN(), Parameter('CN_0'))
# ci
Initial(ci(), Parameter('ci_0'))
# en
Observable('en_obs', en())
Observable('EWG_nT', EWG())
Observable('CN_obs', CN())
#
Observable('EN_obs', EN())
# wg
Observable('wg_obs', wg())
#Observable('CN_obs', CN())
Observable('CI_obs', CI())
Observable('IWG_obs', IWG())
# ptc
Observable('ptc_obs', ptc())
#Observable('CN_obs', CN())
#Observable('CI_obs', CI())
#ci
Observable('ci_obs', ci())

# engrailed
Expression('Hill_CNen', (K_CNen ** nu_CNen) / (K_CNen ** nu_CNen + CN_obs ** nu_CNen))
Expression('Hill_EWG_nT',
           H_en * (EWG_nT * Hill_CNen) ** nu_WGen / (K_WGen ** nu_WGen + (EWG_nT * Hill_CNen) ** nu_WGen))

Rule('synthesize_en_1', None >> en(), Hill_EWG_nT)
Rule('degrade_en_1', en() >> None, Parameter('en_deg', H_en.value))

Expression('EN_translation',
           H_EN * en_obs)
Rule('synthesize_EN', en() >> en() + EN(), EN_translation)
Rule('degrade_EN', EN() >> None, Parameter('EN_deg', H_EN.value))

# wingless
Expression('Hill_CNwg', (K_CNwg ** nu_CNwg) / (K_CNwg ** nu_CNwg + CN_obs ** nu_CNwg))

Expression('Hill_CIwg',
           alpha_CIwg * (CI_obs * Hill_CNwg) ** nu_CIwg / (K_CIwg ** nu_CIwg + (CI_obs * Hill_CNwg) ** nu_CIwg))

Expression('Hill_IWGwg',
           alpha_WGwg * (IWG_obs ** nu_WGwg) / (K_WGwg ** nu_WGwg + IWG_obs ** nu_WGwg))

Expression('Hill_CI_IWG_combined',
           H_wg * (Hill_CIwg + Hill_IWGwg) / (1 + Hill_CIwg + Hill_IWGwg))

Rule('synthesize_wg_1', None >> wg(), Hill_CI_IWG_combined)
Rule('degrade_wg_1', wg() >> None, Parameter('wg_deg',H_wg.value))

# # patched
# Expression('Hill_CNptc', (K_CNptc ** nu_CNptc) / (K_CNptc ** nu_CNptc + CN_obs ** nu_CNptc))
# Expression('Hill_CIptc',
#            H_ptc * (CI_obs * Hill_CNptc) ** nu_CIptc / (K_CIptc ** nu_CIptc + (CI_obs * Hill_CNptc) ** nu_CIptc))
#
# Rule('synthesize_ptc_1', None >> ptc(), Hill_CIptc)
# Rule('degrade_ptc_1', ptc() >> None, Parameter('ptc_deg', H_ptc.value))
#
# #cubitus interruptus
# Expression('Hill_ENci', (K_ENci ** nu_ENci) / (K_ENci ** nu_ENci + EN_obs ** nu_ENci))
# #Parameter('Hill_ENci',1)
# Expression('Hill_Bci',
#           H_ci * (B * Hill_ENci) ** nu_Bci / (K_Bci ** nu_Bci + (B * Hill_ENci) ** nu_Bci))
#
# Rule('synthesize_ci_1', None >> ci(), Hill_Bci)
# Rule('degrade_ci_1', ci() >> None, Parameter('ci_deg', H_ci.value))

#print (model.rules)