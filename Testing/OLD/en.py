from pysb import *
from pysb.core import as_complex_pattern, MonomerPattern, ComplexPattern,\
    ReactionPattern, RuleExpression

def create_model(nu_wg, nu_cn):

    Model()

    ########### FROM REDUCED MODEL ##########
    Parameter('H_en', 1.0/56.195484)
    Parameter('K_WGen', 0.016870342)
    Parameter('nu_WGen', nu_wg) #8.195488)
    Parameter('K_CNen', 0.0015239986)
    Parameter('nu_CNen', nu_cn) #5.929026)
    #########################################
    
    ########### MECHANISTIC MODEL ###########
    
    order_WGen = int(round(nu_WGen.value))
#     print order_WGen
#     print

    order_CNen = int(round(nu_CNen.value))
#     print order_CNen
#     print
    
    mon_gEN = Monomer('gEN', ['ewg%d' % i for i in range(order_WGen)])
    Monomer('mEN')
    mon_pEWG = Monomer('pEWG', ['ewg_l', 'ewg_r', 'gen'] + ['cn%d' % i for i in range(order_CNen)])
    Monomer('pCN', ['cn_l', 'cn_r', 'ewg'])
    
#     print model.monomers
#     print
    
    Initial(MonomerPattern(mon_gEN, {s : None for s in mon_gEN.sites}, None), Parameter('gEN_0', 1))
#     Initial(pEWG(ewg_l=None, ewg_r=None, gen=None), Parameter('pEWG_0'))
    Initial(MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None), Parameter('pEWG_0'))
    Initial(pCN(cn_l=None, cn_r=None, ewg=None), Parameter('pCN_0'))
    
#     for ic in model.initial_conditions:
#         print ic
#     print
    
    Observable('pEWG_tot', pEWG())
    Observable('pEWG_free', pEWG(ewg_l=None, ewg_r=None, gen=None, cn0=None))
#     Observable('pEWG_free', MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None))    
    Observable('pEWG_pCN', pEWG(cn0=ANY))
    Observable('pCN_tot', pCN())
    Observable('pCN_free', pCN(cn_l=None, cn_r=None, ewg=None))
    Observable('mEN_obs', mEN())    
    
#     print model.observables
#     print
    
    Parameter('kf_ewg_olig', 1e8)
    Parameter('kr_ewg_olig', 1e9)
    Parameter('kf_g_ewg', 1e8*10**(order_WGen-1))
    # kr = kf * (Kappa/H)^n / K1^(n-1)
    Parameter('kr_g_ewg', kf_g_ewg.value * \
              (K_WGen.value / H_en.value)**order_WGen * \
              (kf_ewg_olig.value / kr_ewg_olig.value)**(order_WGen-1)) 
    Parameter('ktr_en', 1)
    
#     print kr_g_ewg.value
#     K1 = kr_ewg_olig.value / kf_ewg_olig.value
#     K2 = kr_g_ewg.value / kf_g_ewg.value
#     print (K1**(order_WGen-1)*K2)**(1./order_WGen) * H_en.value
#     print K_WGen.value
#     print
    
    # pEWG oligomerization rules (x + y <> z)
    x = pEWG(ewg_l=None, ewg_r=None, cn0=None)
    for i in range(order_WGen-1):
        reactants = []
        products = [pEWG(ewg_l=None, ewg_r=i, cn0=None)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pEWG(ewg_l=left, ewg_r=right, cn0=None))
            # product pattern
            left = i if j==0 else j-1
            products.append(pEWG(ewg_l=left, ewg_r=right, cn0=None))
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('EWG_oligomerize_1_%d_to_%d' % (i+1, i+2), x + y | z, kf_ewg_olig, kr_ewg_olig)
    
    # Promoter binding rule
    olig_unbound = z.copy()
    bound_complex = z.copy()
    g_unbound = MonomerPattern(mon_gEN, {s : None for s in mon_gEN.sites}, None)
    g_bound = MonomerPattern(mon_gEN, {s : None for s in mon_gEN.sites}, None)
    j = order_WGen-1
    for mp1,mp2,site in zip(olig_unbound.monomer_patterns, 
                            bound_complex.monomer_patterns,
                            g_bound.monomer.sites):
        mp1.site_conditions['gen'] = None
        mp2.site_conditions['gen'] = j
        g_bound.site_conditions[site] = j
        j += 1
    bound_complex.monomer_patterns.append(g_bound)
    Rule('EWG_binds_promoter', olig_unbound + g_unbound | bound_complex, kf_g_ewg, kr_g_ewg)
    
    # Transcription rule
    mp = MonomerPattern(mon_gEN, {s : ANY for s in mon_gEN.sites}, None)
    Rule('EN_transcription', mp >> mp + mEN(), ktr_en)
    
    # Degradation rule
    Rule('mRNA_degradation', mEN() >> None, H_en)
    
#     print model.rules
    
    Parameter('kf_cn_olig', 1e8)
    Parameter('kr_cn_olig', 1e9)
    Parameter('kf_cn_ewg', 1e8*10**(order_CNen-1))
    # kr = kf * (Kappa/H)^n / K1^(n-1)
    Parameter('kr_cn_ewg', kf_cn_ewg.value * \
              (K_CNen.value / H_en.value)**order_CNen * \
              (kf_cn_olig.value / kr_cn_olig.value)**(order_CNen-1))
      
#     print kr_cn_ewg.value
#     K1 = kr_cn_olig.value / kf_cn_olig.value
#     K2 = kr_cn_ewg.value / kf_cn_ewg.value
#     print (K1**(order_CNen-1)*K2)**(1./order_CNen) * H_en.value
#     print K_CNen.value
#     print
    
    # pCN oligomerization rules (x + y <> z)
    x = pCN(cn_l=None, cn_r=None)
    for i in range(order_CNen-1):
        reactants = []
        products = [pCN(cn_l=None, cn_r=i)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pCN(cn_l=left, cn_r=right))
            # product pattern
            left = i if j==0 else j-1
            products.append(pCN(cn_l=left, cn_r=right))
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('CN_oligomerize_1_%d_to_%d' % (i+1, i+2), x + y | z, kf_cn_olig, kr_cn_olig)
    
    # pCN-pEWG binding rule
    olig_unbound = z.copy()
    bound_complex = z.copy()
    ewg_unbound = MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None)
    ewg_bound = MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None)
    j = order_CNen-1
    for mp1,mp2,site in zip(olig_unbound.monomer_patterns, 
                            bound_complex.monomer_patterns,
                            ewg_bound.monomer.sites[3:]):
        mp1.site_conditions['ewg'] = None
        mp2.site_conditions['ewg'] = j
        ewg_bound.site_conditions[site] = j
        j += 1
    bound_complex.monomer_patterns.append(ewg_bound)
    Rule('CN_binds_EWG', olig_unbound + ewg_unbound | bound_complex, kf_cn_ewg, kr_cn_ewg)
    
#     print(model.rules)
#     quit()
