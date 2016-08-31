
from pysb import *
from SPN_macros import transcription
from SPN_macros import transcription_switch
from SPN_macros import translation
from SPN_macros import translation_persistant
from SPN_macros import binary_transformation
from SPN_macros import dimerization
from SPN_macros import proteolysis
from pysb.macros import bind

Model()

# Monomers

Monomer('SLP_1', ['en_g', 'wg_g'])
Monomer('SLP_2', ['en_g', 'wg_g'])
Monomer('SLP_3', ['en_g', 'wg_g'])
Monomer('SLP_4', ['en_g', 'wg_g'])

Monomer('wg_g_1', ['CIA','CIR','SLP'])
Monomer('wg_g_2', ['CIA','CIR','SLP'])
Monomer('wg_g_3', ['CIA','CIR','SLP'])
Monomer('wg_g_4', ['CIA','CIR','SLP'])

Monomer('wg_1')
Monomer('wg_2')
Monomer('wg_3')
Monomer('wg_4')

Monomer('WG_1', ['en_g'])
Monomer('WG_2', ['en_g'])
Monomer('WG_3', ['en_g'])
Monomer('WG_4', ['en_g'])

Monomer('en_g_1', ['WG','SLP'])
Monomer('en_g_2', ['WG','SLP'])
Monomer('en_g_3', ['WG','SLP'])
Monomer('en_g_4', ['WG','SLP'])

Monomer('en_1')
Monomer('en_2')
Monomer('en_3')
Monomer('en_4')

Monomer('EN_1', ['hh_g','ptc_g','ci_g'])
Monomer('EN_2', ['hh_g','ptc_g','ci_g'])
Monomer('EN_3', ['hh_g','ptc_g','ci_g'])
Monomer('EN_4', ['hh_g','ptc_g','ci_g'])

Monomer('hh_g_1', ['EN','CIR'])
Monomer('hh_g_2', ['EN','CIR'])
Monomer('hh_g_3', ['EN','CIR'])
Monomer('hh_g_4', ['EN','CIR'])

Monomer('hh_1')
Monomer('hh_2')
Monomer('hh_3')
Monomer('hh_4')

Monomer('HH_1', ['PTC'])
Monomer('HH_2', ['PTC'])
Monomer('HH_3', ['PTC'])
Monomer('HH_4', ['PTC'])

Monomer('ptc_g_1', ['CIA','CIR','EN'])
Monomer('ptc_g_2', ['CIA','CIR','EN'])
Monomer('ptc_g_3', ['CIA','CIR','EN'])
Monomer('ptc_g_4', ['CIA','CIR','EN'])

Monomer('ptc_1')
Monomer('ptc_2')
Monomer('ptc_3')
Monomer('ptc_4')

Monomer('PTC_1', ['HH', 'SMO'])
Monomer('PTC_2', ['HH', 'SMO'])
Monomer('PTC_3', ['HH', 'SMO'])
Monomer('PTC_4', ['HH', 'SMO'])

Monomer('PH_1', ['SMO'])
Monomer('PH_2', ['SMO'])
Monomer('PH_3', ['SMO'])
Monomer('PH_4', ['SMO'])

Monomer('SMO_1', ['PTC','PH','CI'])
Monomer('SMO_2', ['PTC','PH','CI'])
Monomer('SMO_3', ['PTC','PH','CI'])
Monomer('SMO_4', ['PTC','PH','CI'])

Monomer('ci_g_1', ['EN'])
Monomer('ci_g_2', ['EN'])
Monomer('ci_g_3', ['EN'])
Monomer('ci_g_4', ['EN'])

Monomer('ci_1')
Monomer('ci_2')
Monomer('ci_3')
Monomer('ci_4')

Monomer('CI_1',['SMO'])
Monomer('CI_2',['SMO'])
Monomer('CI_3',['SMO'])
Monomer('CI_4',['SMO'])

Monomer('CIA_1', ['ptc_g', 'wg_g'])
Monomer('CIA_2', ['ptc_g', 'wg_g'])
Monomer('CIA_3', ['ptc_g', 'wg_g'])
Monomer('CIA_4', ['ptc_g', 'wg_g'])

Monomer('CIR_1', ['hh_g','ptc_g', 'wg_g'])
Monomer('CIR_2', ['hh_g','ptc_g', 'wg_g'])
Monomer('CIR_3', ['hh_g','ptc_g', 'wg_g'])
Monomer('CIR_4', ['hh_g','ptc_g', 'wg_g'])

# Initializations

Initial(SLP_1(en_g=None, wg_g=None), Parameter('SLP_1_0', 0))
Initial(SLP_2(en_g=None, wg_g=None), Parameter('SLP_2_0', 0))
Initial(SLP_3(en_g=None, wg_g=None), Parameter('SLP_3_0', 10**5))
Initial(SLP_4(en_g=None, wg_g=None), Parameter('SLP_4_0', 10**5))

Initial(wg_g_1(CIA=None, CIR=None, SLP=None), Parameter('wg_g_1_0', 1))
Initial(wg_g_2(CIA=None, CIR=None, SLP=None), Parameter('wg_g_2_0', 1))
Initial(wg_g_3(CIA=None, CIR=None, SLP=None), Parameter('wg_g_3_0', 1))
Initial(wg_g_4(CIA=None, CIR=None, SLP=None), Parameter('wg_g_4_0', 1))

Initial(wg_1(), Parameter('wg_1_0', 0))
Initial(wg_2(), Parameter('wg_2_0', 0))
Initial(wg_3(), Parameter('wg_3_0', 0))
Initial(wg_4(), Parameter('wg_4_0', 10**5))

Initial(WG_1(en_g=None), Parameter('WG_1_0', 0))
Initial(WG_2(en_g=None), Parameter('WG_2_0', 0))
Initial(WG_3(en_g=None), Parameter('WG_3_0', 0))
Initial(WG_4(en_g=None), Parameter('WG_4_0', 0))

Initial(en_g_1(WG=None, SLP=None), Parameter('en_g_1_0', 1))
Initial(en_g_2(WG=None, SLP=None), Parameter('en_g_2_0', 1))
Initial(en_g_3(WG=None, SLP=None), Parameter('en_g_3_0', 1))
Initial(en_g_4(WG=None, SLP=None), Parameter('en_g_4_0', 1))

Initial(en_1(), Parameter('en_1_0', 10**5))
Initial(en_2(), Parameter('en_2_0', 0))
Initial(en_3(), Parameter('en_3_0', 0))
Initial(en_4(), Parameter('en_4_0', 0))

Initial(EN_1(hh_g=None, ptc_g=None, ci_g=None), Parameter('EN_1_0', 0))
Initial(EN_2(hh_g=None, ptc_g=None, ci_g=None), Parameter('EN_2_0', 0))
Initial(EN_3(hh_g=None, ptc_g=None, ci_g=None), Parameter('EN_3_0', 0))
Initial(EN_4(hh_g=None, ptc_g=None, ci_g=None), Parameter('EN_4_0', 0))

Initial(hh_g_1(EN=None, CIR=None), Parameter('hh_g_1_0', 1))
Initial(hh_g_2(EN=None, CIR=None), Parameter('hh_g_2_0', 1))
Initial(hh_g_3(EN=None, CIR=None), Parameter('hh_g_3_0', 1))
Initial(hh_g_4(EN=None, CIR=None), Parameter('hh_g_4_0', 1))

Initial(hh_1(), Parameter('hh_1_0', 10**5))
Initial(hh_2(), Parameter('hh_2_0', 0))
Initial(hh_3(), Parameter('hh_3_0', 0))
Initial(hh_4(), Parameter('hh_4_0', 0))

Initial(HH_1(PTC=None), Parameter('HH_1_0', 0))
Initial(HH_2(PTC=None), Parameter('HH_2_0', 0))
Initial(HH_3(PTC=None), Parameter('HH_3_0', 0))
Initial(HH_4(PTC=None), Parameter('HH_4_0', 0))

Initial(ptc_g_1(CIA=None, CIR=None, EN=None), Parameter('ptc_g_1_0', 1))
Initial(ptc_g_2(CIA=None, CIR=None, EN=None), Parameter('ptc_g_2_0', 1))
Initial(ptc_g_3(CIA=None, CIR=None, EN=None), Parameter('ptc_g_3_0', 1))
Initial(ptc_g_4(CIA=None, CIR=None, EN=None), Parameter('ptc_g_4_0', 1))

Initial(ptc_1(), Parameter('ptc_1_0', 0))
Initial(ptc_2(), Parameter('ptc_2_0', 10**5))
Initial(ptc_3(), Parameter('ptc_3_0', 10**5))
Initial(ptc_4(), Parameter('ptc_4_0', 10**5))

Initial(PTC_1(HH=None, SMO=None), Parameter('PTC_1_0', 0))
Initial(PTC_2(HH=None, SMO=None), Parameter('PTC_2_0', 0))
Initial(PTC_3(HH=None, SMO=None), Parameter('PTC_3_0', 0))
Initial(PTC_4(HH=None, SMO=None), Parameter('PTC_4_0', 0))

Initial(PH_1(SMO=None), Parameter('PH_1_0', 0))
Initial(PH_2(SMO=None), Parameter('PH_2_0', 0))
Initial(PH_3(SMO=None), Parameter('PH_3_0', 0))
Initial(PH_4(SMO=None), Parameter('PH_4_0', 0))

# In the Boolean model SMO is constitutively expressed. Thus,
# there it has no synthesis or degradation equations, only 
# deactivation/activation via binding/unbinding to PTC 

Initial(SMO_1(PTC=None, PH=None, CI=None), Parameter('SMO_1_0', 10**5))
Initial(SMO_2(PTC=None, PH=None, CI=None), Parameter('SMO_2_0', 10**5))
Initial(SMO_3(PTC=None, PH=None, CI=None), Parameter('SMO_3_0', 10**5))
Initial(SMO_4(PTC=None, PH=None, CI=None), Parameter('SMO_4_0', 10**5))

Initial(ci_g_1(EN=None), Parameter('ci_g_1_0', 1))
Initial(ci_g_2(EN=None), Parameter('ci_g_2_0', 1))
Initial(ci_g_3(EN=None), Parameter('ci_g_3_0', 1))
Initial(ci_g_4(EN=None), Parameter('ci_g_4_0', 1))

Initial(ci_1(), Parameter('ci_1_0', 0))
Initial(ci_2(), Parameter('ci_2_0', 10**5))
Initial(ci_3(), Parameter('ci_3_0', 10**5))
Initial(ci_4(), Parameter('ci_4_0', 10**5))

Initial(CI_1(SMO=None), Parameter('CI_1_0', 0))
Initial(CI_2(SMO=None), Parameter('CI_2_0', 0))
Initial(CI_3(SMO=None), Parameter('CI_3_0', 0))
Initial(CI_4(SMO=None), Parameter('CI_4_0', 0))

Initial(CIA_1(ptc_g=None, wg_g=None), Parameter('CIA_1_0', 0))
Initial(CIA_2(ptc_g=None, wg_g=None), Parameter('CIA_2_0', 0))
Initial(CIA_3(ptc_g=None, wg_g=None), Parameter('CIA_3_0', 0))
Initial(CIA_4(ptc_g=None, wg_g=None), Parameter('CIA_4_0', 0))

Initial(CIR_1(hh_g=None, ptc_g=None, wg_g=None), Parameter('CIR_1_0', 0))
Initial(CIR_2(hh_g=None, ptc_g=None, wg_g=None), Parameter('CIR_2_0', 0))
Initial(CIR_3(hh_g=None, ptc_g=None, wg_g=None), Parameter('CIR_3_0', 0))
Initial(CIR_4(hh_g=None, ptc_g=None, wg_g=None), Parameter('CIR_4_0', 0))

# Observables

Observable('SLP_1_obs', SLP_1())
Observable('SLP_2_obs', SLP_2())
Observable('SLP_3_obs', SLP_3())
Observable('SLP_4_obs', SLP_4())

Observable('wg_1_obs', wg_1())
Observable('wg_2_obs', wg_2())
Observable('wg_3_obs', wg_3())
Observable('wg_4_obs', wg_4())

Observable('WG_1_obs', WG_1())
Observable('WG_2_obs', WG_2())
Observable('WG_3_obs', WG_3())
Observable('WG_4_obs', WG_4())

Observable('en_1_obs', en_1())
Observable('en_2_obs', en_2())
Observable('en_3_obs', en_3())
Observable('en_4_obs', en_4())

Observable('EN_1_obs', EN_1())
Observable('EN_2_obs', EN_2())
Observable('EN_3_obs', EN_3())
Observable('EN_4_obs', EN_4())

Observable('hh_1_obs', hh_1())
Observable('hh_2_obs', hh_2())
Observable('hh_3_obs', hh_3())
Observable('hh_4_obs', hh_4())

Observable('HH_1_obs', HH_1())
Observable('HH_2_obs', HH_2())
Observable('HH_3_obs', HH_3())
Observable('HH_4_obs', HH_4())

Observable('ptc_1_obs',ptc_1())
Observable('ptc_2_obs',ptc_2())
Observable('ptc_3_obs',ptc_3())
Observable('ptc_4_obs',ptc_4())

Observable('PTC_1_obs',PTC_1())
Observable('PTC_2_obs',PTC_2())
Observable('PTC_3_obs',PTC_3())
Observable('PTC_4_obs',PTC_4())

Observable('PH_1_obs', PH_1())
Observable('PH_2_obs', PH_2())
Observable('PH_3_obs', PH_3())
Observable('PH_4_obs', PH_4())

Observable('SMO_1_obs', SMO_1())
Observable('SMO_2_obs', SMO_2())
Observable('SMO_3_obs', SMO_3())
Observable('SMO_4_obs', SMO_4())

Observable('ci_1_obs', ci_1())
Observable('ci_2_obs', ci_2())
Observable('ci_3_obs', ci_3())
Observable('ci_4_obs', ci_4())

Observable('CI_1_obs', CI_1())
Observable('CI_2_obs', CI_2())
Observable('CI_3_obs', CI_3())
Observable('CI_4_obs', CI_4())

Observable('CIA_1_obs', CIA_1())
Observable('CIA_2_obs', CIA_2())
Observable('CIA_3_obs', CIA_3())
Observable('CIA_4_obs', CIA_4())

Observable('CIR_1_obs', CIR_1())
Observable('CIR_2_obs', CIR_2())
Observable('CIR_3_obs', CIR_3())
Observable('CIR_4_obs', CIR_4())

# Parameters

# wg TF binding
Parameter('CIA_wg_kf', 1)
Parameter('CIA_wg_kr', 1)
Parameter('CIR_wg_kf', 1)
Parameter('CIR_wg_kr', 1)
Parameter('SLP_wg_kf', 1)
Parameter('SLP_wg_kr', 1)

# wg transcription
Parameter('wg_ts', 1)
# Parameter('wg_ts_d', 1) # used with 'transcription' macro
Parameter('wg_l', 1)
Parameter('wg_deg', 1)

# WG translation
Parameter('WG_tl', 1)
Parameter('WG_deg', 1)

# en TF binding
Parameter('WG_en_kf', 1)
Parameter('WG_en_kr', 1)
Parameter('SLP_en_kf', 1)
Parameter('SLP_en_kr', 1)

# en transcription
Parameter('en_ts', 1)
Parameter('en_l', 1)
Parameter('en_deg', 1)

# EN translation
Parameter('EN_tl', 1)
Parameter('EN_deg', 1)

# hh TF binding
Parameter('EN_hh_kf', 1)
Parameter('EN_hh_kr', 1)
Parameter('CIR_hh_kf', 1)
Parameter('CIR_hh_kr', 1)

# hh transcription
Parameter('hh_ts', 1)
Parameter('hh_l', 1)
Parameter('hh_deg', 1)

# HH translation
Parameter('HH_tl', 1)
Parameter('HH_deg', 1)

# ptc TF binding
Parameter('CIA_ptc_kf', 1)
Parameter('CIA_ptc_kr', 1)
Parameter('CIR_ptc_kf', 1)
Parameter('CIR_ptc_kr', 1)
Parameter('EN_ptc_kf', 1)
Parameter('EN_ptc_kr', 1)

# ptc transcription
Parameter('ptc_ts', 1)
Parameter('ptc_l', 1)
Parameter('ptc_deg', 1)

# PTC translation
Parameter('PTC_tl', 1)
Parameter('PTC_deg', 1)

# PTC-HH binding and proteolysis
Parameter('PTC_HH_kf', 1)
Parameter('PTC_HH_kr', 1)
Parameter('PTC_HH_proteolysis', 1)

# PTC HH dimerization
Parameter('PH_kf', 1)
Parameter('PH_kr', 1)
Parameter('PH_deg', 1)

# ci TF binding
Parameter('EN_ci_kf', 1)
Parameter('EN_ci_kr', 1)

# ci transcription
Parameter('ci_ts', 1)
Parameter('ci_l', 1)
Parameter('ci_deg', 1)

# CI translation
Parameter('CI_tl', 1)
Parameter('CI_deg', 1)

# SMO-PTC binding
Parameter('SMO_PTC_kf', 1)
Parameter('SMO_PTC_kr', 1)

# SMO-PH binding
Parameter('SMO_PH_kf', 1)
Parameter('SMO_PH_kr', 1)

# SMO-PTC-HH binding
Parameter('SMO_PTC_HH_kf', 1)
Parameter('SMO_PTC_HH_kr', 1)

# CI cleavage
Parameter('SMO_CI_kf', 1)
Parameter('SMO_CI_kr', 1)
Parameter('CIA_kc1', 1)
Parameter('CIA_kd1', 1)
Parameter('CIR_kc2', 1)
Parameter('CIR_kd2', 1)

# Rules

# wg
# ==
# Boolean function: wg = (CIA and SLP and not CIR) or (wg and (CIA or SLP) and not CIR)
# Boolean rules:
#     wg(1) + CIR(1) >> wg(0) + CIR(1)
#     wg(1) + CIR(0) + SLP(0) + CIA(0) >> wg(0) + CIR(0) + SLP(0) + CIA(0)
#     wg(0) + CIR(0) + SLP(1) + CIA(1) >> wg(1) + CIR(0) + SLP(1) + CIA(1)

# Wingless transcription is inhibited by CIR (dominantly) and activated by CIA and SLP.
# All are transcription factors. Note the second part of the Boolean function. It says that the 
# absense of CIR and presence of one of CIA or SLP is enough to preserve, but not change, the state
# of wg. This can be represented as transcription with distinct transcription rates for each
# combination of trascription factors. There are 2^n combinations for n TFs.  

# transcription(gene, g_sites, TF_list, TF_sites, klist, product, tlist, deg_r):

# transcription(wg_g_1, ['CIA','CIR','SLP'], [CIA_1, CIR_1, SLP_1], ['CIA','CIR','SLP'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_1, [wg_l, wg_ts, wg_l, wg_l, wg_ts, wg_ts_d, wg_l, wg_l], wg_deg)
# transcription(wg_g_2, ['CIA','CIR','SLP'], [CIA_2, CIR_2, SLP_2], ['CIA','CIR','SLP'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_2, [wg_l, wg_ts, wg_l, wg_l, wg_ts, wg_ts_d, wg_l, wg_l], wg_deg)
# transcription(wg_g_3, ['CIA','CIR','SLP'], [CIA_3, CIR_3, SLP_3], ['CIA','CIR','SLP'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_3, [wg_l, wg_ts, wg_l, wg_l, wg_ts, wg_ts_d, wg_l, wg_l], wg_deg)
# transcription(wg_g_4, ['CIA','CIR','SLP'], [CIA_4, CIR_4, SLP_4], ['CIA','CIR','SLP'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_4, [wg_l, wg_ts, wg_l, wg_l, wg_ts, wg_ts_d, wg_l, wg_l], wg_deg)

# However, the continuation of state is not explicit in the Boolean rules representation generated by 
# our algorithm. In this case we can use a simple 'switch-like' macro in which there is a transcription
# rate for the 'correct' combination of TFs and a leakage rate for all other combinations. The 'correct'
# combination is provided as a binary string.

# transcription_switch(gene, g_sites, TF_list, TF_sites, klist, product, tlist, deg_r, switch_state):

transcription_switch(wg_g_1, ['CIA','CIR','SLP'], [CIA_1, CIR_1, SLP_1], ['wg_g','wg_g','wg_g'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_1, [wg_ts, wg_l], wg_deg, '101')
transcription_switch(wg_g_2, ['CIA','CIR','SLP'], [CIA_2, CIR_2, SLP_2], ['wg_g','wg_g','wg_g'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_2, [wg_ts, wg_l], wg_deg, '101')
transcription_switch(wg_g_3, ['CIA','CIR','SLP'], [CIA_3, CIR_3, SLP_3], ['wg_g','wg_g','wg_g'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_3, [wg_ts, wg_l], wg_deg, '101')
transcription_switch(wg_g_4, ['CIA','CIR','SLP'], [CIA_4, CIR_4, SLP_4], ['wg_g','wg_g','wg_g'], [CIA_wg_kf, CIA_wg_kr, CIR_wg_kf, CIR_wg_kr, SLP_wg_kf, SLP_wg_kr], wg_4, [wg_ts, wg_l], wg_deg, '101')

# WG
# ==
# Boolean function: WG* = wg
# Boolean rules:
#     WG(0) + wg(1) >> WG(1) + wg(1)
#     WG(1) + wg(0) >> WG(0) + wg(0)

# This can be interpreted as simple translation of the mRNA wg into the protein WG.

# translation(mRNA, protein, klist):

translation(wg_1, WG_1, [WG_tl, WG_deg])
translation(wg_2, WG_2, [WG_tl, WG_deg])
translation(wg_3, WG_3, [WG_tl, WG_deg])
translation(wg_4, WG_4, [WG_tl, WG_deg])

# en
# ==
# Boolean function: en* = WG+- and not SLP
# Boolean rules:
#     en(1) + SLP(1) >> en(0) + SLP(1)
#     en(0) + SLP(0) + WG+(1) >> en(1) + SLP(0) + WG+(1)
#     en(0) + SLP(0) + WG+(0) + WG-(1) >> en(1) + SLP(0) + WG+(0) + WG-(1)
#     en(1) + SLP(0) + WG+(0) + WG-(0) >> en(0) + SLP(0) + WG+(0) + WG-(0)

# This is interpreted as transcription. In this case the wingless protein WG from either adjacent cell
# can bind the gene en_g. In reality this is mediated through the Frizzled protein, but this is not 
# represented in the model. The macro accounts for the fact that WG from either adjacent cell binds 
# to the same transcription factor binding site on the gene.

transcription_switch(en_g_1, ['WG','SLP'], [[WG_4, WG_2], SLP_1], ['en_g','en_g'], [WG_en_kf, WG_en_kr, SLP_en_kf, SLP_en_kr], en_1, [en_ts, en_l], en_deg, '10')
transcription_switch(en_g_2, ['WG','SLP'], [[WG_1, WG_3], SLP_2], ['en_g','en_g'], [WG_en_kf, WG_en_kr, SLP_en_kf, SLP_en_kr], en_2, [en_ts, en_l], en_deg, '10')
transcription_switch(en_g_3, ['WG','SLP'], [[WG_2, WG_4], SLP_3], ['en_g','en_g'], [WG_en_kf, WG_en_kr, SLP_en_kf, SLP_en_kr], en_3, [en_ts, en_l], en_deg, '10')
transcription_switch(en_g_4, ['WG','SLP'], [[WG_3, WG_1], SLP_4], ['en_g','en_g'], [WG_en_kf, WG_en_kr, SLP_en_kf, SLP_en_kr], en_4, [en_ts, en_l], en_deg, '10')

# EN
# ==
# Boolean function: EN* = en
# Boolean rules:
#     EN(0) + en(1) >> EN(1) + en(1)
#     EN(1) + en(0) >> EN(0) + en(0)

translation(en_1, EN_1, [EN_tl, EN_deg])
translation(en_2, EN_2, [EN_tl, EN_deg])
translation(en_3, EN_3, [EN_tl, EN_deg])
translation(en_4, EN_4, [EN_tl, EN_deg])

# hh
# ==
# Boolean function: hh* = EN and not CIR
# Boolean rules:
#     hh(1) + CIR(1) >> hh(0) + CIR(1)
#     hh(0) + CIR(0) + EN(1) >> hh(1) + CIR(0) + EN(1)
#     hh(1) + CIR(0) + EN(0) >> hh(0) + CIR(0) + EN(0)

transcription_switch(hh_g_1, ['EN', 'CIR'], [EN_1, CIR_1], ['hh_g', 'hh_g'], [EN_hh_kf, EN_hh_kr, CIR_hh_kf, CIR_hh_kr], hh_1, [hh_ts, hh_l], hh_deg, '01')
transcription_switch(hh_g_2, ['EN', 'CIR'], [EN_2, CIR_2], ['hh_g', 'hh_g'], [EN_hh_kf, EN_hh_kr, CIR_hh_kf, CIR_hh_kr], hh_2, [hh_ts, hh_l], hh_deg, '01')
transcription_switch(hh_g_3, ['EN', 'CIR'], [EN_3, CIR_3], ['hh_g', 'hh_g'], [EN_hh_kf, EN_hh_kr, CIR_hh_kf, CIR_hh_kr], hh_3, [hh_ts, hh_l], hh_deg, '01')
transcription_switch(hh_g_4, ['EN', 'CIR'], [EN_4, CIR_4], ['hh_g', 'hh_g'], [EN_hh_kf, EN_hh_kr, CIR_hh_kf, CIR_hh_kr], hh_4, [hh_ts, hh_l], hh_deg, '01')

# HH
# ==
# Boolean function: HH* = hh
# Boolean rules:
#     HH(0) + hh(1) >> HH(1) + hh(1)
#     HH(1) + hh(0) >> HH(0) + hh(0)

translation(hh_1, HH_1, [HH_tl, HH_deg])
translation(hh_2, HH_2, [HH_tl, HH_deg])
translation(hh_3, HH_3, [HH_tl, HH_deg])
translation(hh_4, HH_4, [HH_tl, HH_deg])

# ptc
# ===
# Boolean function: ptc* = CIA and not CIR and not EN
# Boolean rules:
#     ptc(1) + CIA(1) + CIR(1) >> ptc(0) + CIA(1) + CIR(1)
#     ptc(1) + CIA(1) + CIR(0) + EN(1) >> ptc(0) + CIA(1) + CIR(0) + EN(1)
#     ptc(0) + CIA(1) + CIR(0) + EN(0) >> ptc(1) + CIA(1) + CIR(0) + EN(0)
#     ptc(1) + CIA(0) >> ptc(0) + CIA(0)

transcription_switch(ptc_g_1, ['CIA','CIR','EN'], [CIA_1, CIR_1, EN_1], ['ptc_g','ptc_g','ptc_g'], [CIA_ptc_kf, CIA_ptc_kr, CIR_ptc_kf, CIR_ptc_kr, EN_ptc_kf, EN_ptc_kr], ptc_1, [ptc_ts, ptc_l], ptc_deg, '100')
transcription_switch(ptc_g_2, ['CIA','CIR','EN'], [CIA_2, CIR_2, EN_2], ['ptc_g','ptc_g','ptc_g'], [CIA_ptc_kf, CIA_ptc_kr, CIR_ptc_kf, CIR_ptc_kr, EN_ptc_kf, EN_ptc_kr], ptc_2, [ptc_ts, ptc_l], ptc_deg, '100')
transcription_switch(ptc_g_3, ['CIA','CIR','EN'], [CIA_3, CIR_3, EN_3], ['ptc_g','ptc_g','ptc_g'], [CIA_ptc_kf, CIA_ptc_kr, CIR_ptc_kf, CIR_ptc_kr, EN_ptc_kf, EN_ptc_kr], ptc_3, [ptc_ts, ptc_l], ptc_deg, '100')
transcription_switch(ptc_g_4, ['CIA','CIR','EN'], [CIA_4, CIR_4, EN_4], ['ptc_g','ptc_g','ptc_g'], [CIA_ptc_kf, CIA_ptc_kr, CIR_ptc_kf, CIR_ptc_kr, EN_ptc_kf, EN_ptc_kr], ptc_4, [ptc_ts, ptc_l], ptc_deg, '100')

# PTC
# ===
# Boolean function: PTC = ptc or (PTC and not HH+-)
# Boolean rules:
#     PTC(0) + ptc(1) >> PTC(1) + ptc(0)
#     PTC(1) + ptc(0) + HH+(1) >> PTC(0) + ptc(0) + HH+(1)
#     PTC(1) + ptc(0) + HH+(0) + HH-(1) >> PTC(1) + ptc(0) + HH+(0) + HH-(1)

# Here they expect PTC to be translated from ptc and to persist until bound to HH. PTC and HH bind
# to form the dimer PH. This is not apparent given this Boolean function and these Boolean rules
# for PTC. Thus it appears as though HH is degrading PTC giving us the following hypothesized
# mass action rules.

#    ptc >> ptc + PTC
#    PTC + HH <> PTC.HH
#    PTC.HH >> HH

# The translation of ptc is subsumed into the interpretation of PH below and the
# degradation of PTC by HH directly conflicts with it. We thus disregard this
# interpretation of PTC and instead consider these rules only when interpreting
# the rules for PH. 

# translation_persistant(mRNA, protein, klist):

# translation_persistant(ptc_1, PTC_1, [PTC_tl, PTC_deg])
# translation_persistant(ptc_2, PTC_2, [PTC_tl, PTC_deg])
# translation_persistant(ptc_3, PTC_3, [PTC_tl, PTC_deg])
# translation_persistant(ptc_4, PTC_4, [PTC_tl, PTC_deg])

# proteolysis(enzyme, e_site, target, t_site, k_list):

# proteolysis(HH_4, 'PTC', PTC_1, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_2, 'PTC', PTC_1, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_1, 'PTC', PTC_2, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_3, 'PTC', PTC_2, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_2, 'PTC', PTC_3, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_4, 'PTC', PTC_3, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_3, 'PTC', PTC_4, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])
# proteolysis(HH_1, 'PTC', PTC_4, 'HH', [PTC_HH_kf, PTC_HH_kr, PTC_HH_proteolysis])

# PH
# ==
# Boolean function: PH = PTC and HH+- (instantaneous)
#                      = (ptc or (PTC and not HH+-)) and hh+- (substituted)
# Boolean rules: (instantaneous)
#     PH(0) + PTC(1) + HH+(1) >> PH(1) + PTC(1) + HH+(1)
#     PH(0) + PTC(1) + HH+(0) + HH-(1) >> PH(1) + PTC(1) + HH+(0) + HH-(1)
#     PH(1) + PTC(1) + HH+(0) + HH-(0) >> PH(0) + PTC(1) + HH+(0) + HH-(0)
#     PH(1) + PTC(0) >> PH(0) + PTC(0)
# Boolean rules: (instantaneous-condensed): (hh+- --> hh, HH+- --> HH)
#     PH(0) + PTC(1) + HH(1) >> PH(1) + PTC(1) + HH(1)
#     PH(1) + PTC(1) + HH(0) >> PH(0) + PTC(1) + HH(0)
#     PH(1) + PTC(0) >> PH(0) + PTC(0)
# Boolean rules: (substituted)
#     PH(0) + ptc(1) + hh+(1) >> PH(1) + ptc(1) + hh+(1)
#     PH(0) + ptc(1) + hh+(0) + hh-(1) >> PH(1) + ptc(1) + hh+(0) + hh-(1)
#     PH(1) + ptc(1) + hh+(0) + hh-(0) >> PH(0) + ptc(1) + hh+(0) + hh-(0)
#     PH(1) + ptc(0) + HH+(1) >> PH(0) + ptc(0) + HH+(1)
#     PH(1) + ptc(0) + HH+(0) + HH-(1) >> PH(0) + ptc(0) + HH+(0) + HH-(1)
#     PH(0) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(1) >> PH(1) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(1)
#     PH(0) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(0) + hh-(1) >> PH(1) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(0) + hh-(1)
#     PH(1) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(0) + hh-(0) >> PH(0) + ptc(0) + HH+(0) + HH-(0) + PTC(1) + hh+(0) + hh-(0)
#     PH(1) + ptc(0) + HH+(0) + HH-(0) + PTC(0) >> PH(0) + ptc(0) + HH+(0) + HH-(0) + PTC(0)
# Boolean rules (substituted-condensed): (hh+- --> hh, HH+- --> HH)
#     PH(0) + ptc(1) + hh(1) >> PH(1) + ptc(1) + hh(1)
#     PH(1) + ptc(1) + hh(0) >> PH(0) + ptc(1) + hh(0)
#     PH(1) + ptc(0) + HH(1) >> PH(0) + ptc(0) + HH(1)
#     PH(0) + ptc(0) + HH(0) + PTC(1) + hh(1) >> PH(1) + ptc(0) + HH(0) + PTC(1) + hh(1)
#     PH(1) + ptc(0) + HH(0) + PTC(1) + hh(0) >> PH(0) + ptc(0) + HH(0) + PTC(1) + hh(0)
#     PH(1) + ptc(0) + HH(0) + PTC(0) >> PH(0) + ptc(0) + HH(0) + PTC(0)

# The Boolean functions for PH is taken to be instantaneous by Albert et al. This is 
# because they expect reactions like binding to be fast compared to transcription
# and translation. Because Booleannet does not allow for instantaneous updates we must
# substitute the update functions for each of the variables.

# Example:
#     A* = B         (normal)
#     C* = A* = B    (instantaneous and substituted)

# We can therefore interpret the function in either its instantaneous or substituted form. 
# The instantaneous update function (PH = PTC and HH+-) can be easily interpreted as a 
# dimerization between PTC and HH to produce PH. The mass action rules are

#     PTC + HH  <> PH
#     PH >> None

# If we consider the Boolean rules for the substituted function 
# (PH = (ptc or (PTC and not HH+-)) and hh+-) and assume we are dealing with a two 
# step process we can say that PH is produced when

#     1. both ptc and hh are actively translated
#     2. hh is actively translated and PTC is present but HH is not

# The active ptc and hh translation at the previous timestep would seem to indicate 
# that PTC and HH are interacting to form PH. Translated hh can also 
# interact with existing PTC but only when HH is not present. Now consider the rule
# for PTC: PTC(1) + ptc(0) + HH(1) >> PTC(0) + ptc(0) + HH(1). One interpretation 
# of this is that existing HH has bound all currently available PTC leaving
# none to bind to the newly formed PTC. Altogether then we have the transcription
# of ptc and hh, and the dimerization of PTC and HH to form PH as seen above. 
# The translation of hh is already accounted for. Thus we have

translation(ptc_1, PTC_1, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4, [PTC_tl, PTC_deg])

# dimerization(s1, site_1, s2, site_2, product, k_list):

dimerization(PTC_1, 'HH', [HH_4, HH_2], 'PTC', PH_1, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2, 'HH', [HH_1, HH_3], 'PTC', PH_2, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3, 'HH', [HH_2, HH_4], 'PTC', PH_3, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4, 'HH', [HH_3, HH_1], 'PTC', PH_4, [PH_kf, PH_kr, PH_deg])

# SMO
# ===
# Boolean function: SMO = not PTC or HH+- (instantaneous)
#                         not (ptc or (PTC and not HH+-)) or hh+- (substituted)
# Boolean rules (instantaneous
#     SMO(0) + HH+(1) >> SMO(1) + HH+(1)
#     SMO(0) + HH+(0) + HH-(1) >> SMO(1) + HH+(0) + HH-(1)
#     SMO(1) + HH+(0) + HH-(0) + PTC(1) >> SMO(0) + HH+(0) + HH-(0) + PTC(1)
#     SMO(0) + HH+(0) + HH-(0) + PTC(0) >> SMO(1) + HH+(0) + HH-(0) + PTC(0)
# Boolean rules (instantaneous-condensed): (hh+- --> hh, HH+- --> HH)
#     SMO(0) + HH(1) >> SMO(1) + HH(1)
#     SMO(1) + HH(0) + PTC(1) >> SMO(0) + HH(0) + PTC(1)
#     SMO(0) + HH(0) + PTC(0) >> SMO(1) + HH(0) + PTC(0)
# Boolean rules (substituted)
#     SMO(0) + hh+(1) >> SMO(1) + hh+(1)
#     SMO(0) + hh+(0) + hh-(1) >> SMO(1) + hh+(0) + hh-(1)
#     SMO(1) + hh+(0) + hh-(0) + ptc(1) >> SMO(0) + hh+(0) + hh-(0) + ptc(1)
#     SMO(0) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(1) >> SMO(1) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(1)
#     SMO(0) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(0) + HH+(1) >> SMO(1) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(0) + HH+(1)
#     SMO(1) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(0) + HH+(0) >> SMO(0) + hh+(0) + hh-(0) + ptc(0) + PTC(1) + HH-(0) + HH+(0)
#     SMO(0) + hh+(0) + hh-(0) + ptc(0) + PTC(0) >> SMO(1) + hh+(0) + hh-(0) + ptc(0) + PTC(0)
# Boolean rules (substituted-condensed): (hh+- --> hh, HH+- --> HH)
#     SMO(0) + hh(1) >> SMO(1) + hh(1)
#     SMO(1) + hh(0) + ptc(1) >> SMO(0) + hh(0) + ptc(1)
#     SMO(0) + hh(0) + ptc(0) + PTC(1) + HH(1) >> SMO(1) + hh(0) + ptc(0) + PTC(1) + HH(1)
#     SMO(1) + hh(0) + ptc(0) + PTC(1) + HH(0) >> SMO(0) + hh(0) + ptc(0) + PTC(1) + HH(0)
#     SMO(0) + hh(0) + ptc(0) + PTC(0) >> SMO(1) + hh(0) + ptc(0) + PTC(0)

# SMO is another one where the reaction is taken to be instantaneous. From
# the rules above we have that PTC and HH dimerize (fast) to form PH. Given 
# this, we can interpret the (instantaneous) Boolean rules for SMO as inactive
# only in the presence of free PTC. Here, activity status is equivalent
# to its bound state and will only be apparent when considering those species 
# that are targets of SMO. Thus we need only consider the various bound states
# of SMO.

#     PTC + SMO <> PTC.SMO
#     PTC.SMO + HH <> PH.SMO
#     PH + SMO <> PH.SMO

# For the substituted Boolean function (and rules) we again assume a two step process. 
# We have already determined that PTC and HH bind to form the dimer PH. Fron the 
# Boolean rules for SMO is see that SMO is inactivated when

#     1. ptc is actively translated but hh is not. 
#     2. PTC but not HH currently present and neither is currently transcribed

# SMO is activated whenever

#     3. hh is currently translated
#     4. PTC and HH are both present but neither is translated
#     5. PTC is both absent and untranslated

# Given this and our interpretation for PH above we conclude that SMO is inactivated
# in the presence of free PTC. As above this will only be apparent when considering 
# targets for SMO. For now we need only consider all the binding states. Of course,
# we also don't need to account for the translation of ptc and hh again. Thus we have

bind(SMO_1, 'PTC', PTC_1, 'SMO', [SMO_PTC_kf, SMO_PTC_kr])
bind(SMO_2, 'PTC', PTC_2, 'SMO', [SMO_PTC_kf, SMO_PTC_kr])
bind(SMO_3, 'PTC', PTC_3, 'SMO', [SMO_PTC_kf, SMO_PTC_kr])
bind(SMO_4, 'PTC', PTC_4, 'SMO', [SMO_PTC_kf, SMO_PTC_kr])

bind(SMO_1, 'PH', PH_1, 'SMO', [SMO_PH_kf, SMO_PH_kr])
bind(SMO_2, 'PH', PH_2, 'SMO', [SMO_PH_kf, SMO_PH_kr])
bind(SMO_3, 'PH', PH_3, 'SMO', [SMO_PH_kf, SMO_PH_kr])
bind(SMO_4, 'PH', PH_4, 'SMO', [SMO_PH_kf, SMO_PH_kr])

# In order to be consistent with the Boolean model I have renamed the bound entity PTC % HH as PH throughout.

Rule('bound_dimerization_SMO_1_PTC_1_HH_4', SMO_1(PTC=1) % PTC_1(SMO=1, HH=None) + HH_4(PTC=None) <> SMO_1(PTC=1) % PH_1(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_1_PTC_1_HH_2', SMO_1(PTC=1) % PTC_1(SMO=1, HH=None) + HH_2(PTC=None) <> SMO_1(PTC=1) % PH_1(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_2_PTC_2_HH_1', SMO_2(PTC=1) % PTC_2(SMO=1, HH=None) + HH_1(PTC=None) <> SMO_2(PTC=1) % PH_2(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_2_PTC_2_HH_3', SMO_2(PTC=1) % PTC_2(SMO=1, HH=None) + HH_3(PTC=None) <> SMO_2(PTC=1) % PH_2(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_3_PTC_3_HH_2', SMO_3(PTC=1) % PTC_3(SMO=1, HH=None) + HH_2(PTC=None) <> SMO_3(PTC=1) % PH_3(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_3_PTC_3_HH_4', SMO_3(PTC=1) % PTC_3(SMO=1, HH=None) + HH_4(PTC=None) <> SMO_3(PTC=1) % PH_3(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_4_PTC_4_HH_3', SMO_4(PTC=1) % PTC_4(SMO=1, HH=None) + HH_3(PTC=None) <> SMO_4(PTC=1) % PH_4(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)
Rule('bound_dimerization_SMO_4_PTC_4_HH_1', SMO_4(PTC=1) % PTC_4(SMO=1, HH=None) + HH_1(PTC=None) <> SMO_4(PTC=1) % PH_4(SMO=1), SMO_PTC_HH_kf, SMO_PTC_HH_kr)

# ci
# ==
# Boolean function: ci* = not EN
# Boolean rules:
#     ci(1) + EN(1) >> ci(0) + EN(1)
#     ci(0) + EN(0) >> ci(1) + EN(0)

transcription_switch(ci_g_1, ['EN'], [EN_1], ['ci_g'], [EN_ci_kf, EN_ci_kr], ci_1, [ci_ts, ci_l], ci_deg, '0')
transcription_switch(ci_g_2, ['EN'], [EN_2], ['ci_g'], [EN_ci_kf, EN_ci_kr], ci_2, [ci_ts, ci_l], ci_deg, '0')
transcription_switch(ci_g_3, ['EN'], [EN_3], ['ci_g'], [EN_ci_kf, EN_ci_kr], ci_3, [ci_ts, ci_l], ci_deg, '0')
transcription_switch(ci_g_4, ['EN'], [EN_4], ['ci_g'], [EN_ci_kf, EN_ci_kr], ci_4, [ci_ts, ci_l], ci_deg, '0')

# CI
# ==
# Boolean function: CI* = ci
# Boolean rules:
#     CI(0) + ci(1) >> CI(1) + ci(1)
#     CI(1) + ci(0) >> CI(0) + ci(0)

translation(ci_1, CI_1, [CI_tl, CI_deg])
translation(ci_2, CI_2, [CI_tl, CI_deg])
translation(ci_3, CI_3, [CI_tl, CI_deg])
translation(ci_4, CI_4, [CI_tl, CI_deg])

# CIA
# ===
# Boolean function: CIA = CI and (SMO or hh+-)
# Boolean rules:
#     CIA(0) + CI(1) + hh+(1) >> CIA(1) + CI(1) + hh+(1)
#     CIA(0) + CI(1) + hh+(0) + SMO(1) >> CIA(1) + CI(1) + hh+(0) + SMO(1)
#     CIA(0) + CI(1) + hh+(0) + SMO(0) + hh-(1) >> CIA(1) + CI(1) + hh+(0) + SMO(0) + hh-(1)
#     CIA(1) + CI(1) + hh+(0) + SMO(0) + hh-(0) >> CIA(0) + CI(1) + hh+(0) + SMO(0) + hh-(0)
#     CIA(1) + CI(0) >> CIA(0) + CI(0)
# Boolean rules: condensed (hh+- --> hh)
#     CIA(0) + CI(1) + hh(1) >> CIA(1) + CI(1) + hh(1)
#     CIA(0) + CI(1) + hh(0) + SMO(1) >> CIA(1) + CI(1) + hh(0) + SMO(1)
#     CIA(1) + CI(1) + hh(0) + SMO(0) >> CIA(0) + CI(1) + hh(0) + SMO(0)
#     CIA(1) + CI(0) >> CIA(0) + CI(0)

# CIA activation clearly requires CI so, for now, let's interpret this as some 
# sort of transformation from CI to CIA. In the presence of CI and absence of 
# hh we can interpret active SMO as an inducer of this transformation. In interpreting
# the rules for SMO and PTC above we determined that actively translating hh 
# results in HH that binds PTC and eliminates PTCs ability to inhibit SMO. Thus the
# presence of hh is interpreted here as activating currently present SMO. hh 
# translation, HH binding to PTC, and the interaction of the resulting complex
# with SMO have been covered above. Here we need only consider the binding of 
# SMO to CI and the resulting transformation. 

#     CI + SMO > CIA + SMO

# CIR
# ===
# Boolean function: CIR = CI and not SMO and not hh+-)
# Boolean rules:
#     CIR(1) + SMO(1) >> CIR(0) + SMO(1)
#     CIR(1) + SMO(0) + CI(1) + hh+(1) >> CIR(0) + SMO(0) + CI(1) + hh+(1)
#     CIR(1) + SMO(0) + CI(1) + hh+(0) + hh-(1) >> CIR(0) + SMO(0) + CI(1) + hh+(0) + hh-(1)
#     CIR(0) + SMO(0) + CI(1) + hh+(0) + hh-(0) >> CIR(1) + SMO(0) + CI(1) + hh+(0) + hh-(0)
#     CIR(1) + SMO(0) + CI(0) >> CIR(0) + SMO(0) + CI(0)
# Boolean rules: condensed (hh+- --> hh)
#     CIR(1) + SMO(1) >> CIR(0) + SMO(1)
#     CIR(1) + SMO(0) + CI(1) + hh(1) >> CIR(0) + SMO(0) + CI(1) + hh(1)
#     CIR(0) + SMO(0) + CI(1) + hh(0) >> CIR(1) + SMO(0) + CI(1) + hh(0)
#     CIR(1) + SMO(0) + CI(0) >> CIR(0) + SMO(0) + CI(0)

# We recognize CIR as the exact opposite of CIA. It too requires CI but is formed
# only when neither SMO or hh are active. Apparently there is some other player(s)
# that carries this out, but are not present in this model. Thus, assuming this
# convertion is uni-directional
# we have 

#     CI > CIR

# Given that CIA and CIR are the two possible fates for CI (simple degradation is covered above)
# the following macro for a binary transformation is used to cover both. Note that all other 
# binding sites are assumed to have a value of None, which may not be general enough for other
# cases.

# binary_transformation(enzyme, e_site, substrate, s_site, product_1, product_2, k_list):

binary_transformation(SMO_1, 'CI', CI_1, 'SMO', CIA_1, CIR_1, [SMO_CI_kf, SMO_CI_kr, CIA_kc1, CIR_kc2, CIA_kd1, CIR_kd2])
binary_transformation(SMO_2, 'CI', CI_2, 'SMO', CIA_2, CIR_2, [SMO_CI_kf, SMO_CI_kr, CIA_kc1, CIR_kc2, CIA_kd1, CIR_kd2])
binary_transformation(SMO_3, 'CI', CI_3, 'SMO', CIA_3, CIR_3, [SMO_CI_kf, SMO_CI_kr, CIA_kc1, CIR_kc2, CIA_kd1, CIR_kd2])
binary_transformation(SMO_4, 'CI', CI_4, 'SMO', CIA_4, CIR_4, [SMO_CI_kf, SMO_CI_kr, CIA_kc1, CIR_kc2, CIA_kd1, CIR_kd2])

# print model.rules

