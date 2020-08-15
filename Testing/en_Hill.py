from pysb import *

Model()

Parameter('H_en', 1.0/56.195484)
Parameter('K_WGen', 0.016870342)
Parameter('nu_WGen', 8)#2) #8.195488)
Parameter('K_CNen', 0.0015239986)
Parameter('nu_CNen', 6)#2) # 5.929026)

Monomer('en') # engrailed mRNA
Monomer('EWG') # extracellular wingless protein
Monomer('CN') # cleaved cubitus interuptus protein

Initial(en(), Parameter('en_0'))
Initial(EWG(), Parameter('EWG_0'))
Initial(CN(), Parameter('CN_0'))

Observable('en_obs', en())
Observable('EWG_nT', EWG())
Observable('CN_obs', CN())
    
Expression('Hill_CNen', (K_CNen**nu_CNen) / (K_CNen**nu_CNen + CN_obs**nu_CNen))
Expression('Hill_EWG_nT', H_en*(EWG_nT*Hill_CNen)**nu_WGen / (K_WGen**nu_WGen + (EWG_nT*Hill_CNen)**nu_WGen))

Rule('synthesize_en_1', None >> en(), Hill_EWG_nT)
Rule('degrade_en_1', en() >> None, Parameter('deg',H_en.value))


