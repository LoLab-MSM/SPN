from pysb import *
from pysb.core import as_complex_pattern, MonomerPattern, ComplexPattern, \
    ReactionPattern, RuleExpression, as_reaction_pattern
import sys


def create_model(nu_wg, nu_ci, nu_cn):
    Model()

    ########### FROM REDUCED MODEL ##########
    Parameter('H_wg', 1.0 / 65.74263)
    Parameter('K_WGwg', 0.010839762)
    Parameter('nu_WGwg', nu_wg) #5.2690034)
    Parameter('alpha_WGwg', 1.2893039)
    Parameter('K_CIwg', 0.1358894)
    Parameter('nu_CIwg', nu_ci) #3.0448556)
    Parameter('alpha_CIwg', 5.354691)
    Parameter('K_CNwg', 0.08410451)
    Parameter('nu_CNwg', nu_cn) #2.2106833)
    ### Extras needed for comparing reduced to mechanistic
    #Parameter('H_wg', 1.0 / 65.74263) #already above
    Parameter('H_IWG', 1.0 / 95.70869)
    Parameter('H_EWG', 1.0 / 47.52389)
    Parameter('H_ci', 1.0 / 88.05598)
    Parameter('H_CI', 1.0 / 20.967995)
    Parameter('H_CN', 1.0 / 17.765957)
    #########################################

    ########### MECHANISTIC MODEL ###########

    order_WGwg = int(round(nu_WGwg.value))  # 5
    #     print order_WGen
    #     print

    order_CIwg = int(round(nu_CIwg.value))  # 3
    #     print order_CNen
    #     print

    order_CNwg = int(round(nu_CNwg.value))  # 2
    #     print order_CNen
    #     print

    # instead of doing x for x in range(whatever), need to get the right number into the binding reaction
    # (so not the monomer but the rule)
    #Monomer('pEWG', ['ewg_l', 'ewg_r', 'gEN', 'X1'])  # + ['cn%d' % i for i in range(order_CNen)])
    Monomer('pCI', ['ci_l', 'ci_r', 'X2', 'gA'])  # binding sites relevant for en only
    Monomer('pCN', ['cn_l', 'cn_r', 'gX2'])  # binding sites relevant for wg only
    Monomer('gX2', ['CN_oligomer'])
    Monomer('mX2')
    Monomer('pX2', ['x2_l', 'x2_r', 'CI'])
    Monomer('CI_X2_oligomer', ['gA'])
    Monomer('gA', ['ci_X2'])
    Monomer('mA')
    Monomer('pA', ['gWng'])
    Monomer('pIWG',['iwg_l', 'iwg_r', 'enzyme'])
    Monomer('IWG_oligomer',['enzyme'])
    Monomer('enzyme',['IWG_oligomer'])
    Monomer('pY', ['gWng'])
    Monomer('gWng', ['pA', 'pY'])
    Monomer('wng')

    # print (model.monomers)

    Initial(pCI(ci_l=None, ci_r=None, X2=None, gA=None), Parameter('pCI_0'))
    Initial(pCN(cn_l=None, cn_r=None, gX2=None), Parameter('pCN_0'))
    Initial(pIWG(iwg_l=None, iwg_r=None, enzyme=None), Parameter('pIWG_0'))
    Initial(gX2(CN_oligomer=None), Parameter('gX2_0',1))
    Initial(mX2(), Parameter('mX2_0'))
    Initial(pX2(x2_l=None, x2_r=None, CI=None), Parameter('pX2_0'))
    Initial(CI_X2_oligomer(gA=None), Parameter('CI_X2_olig_0'))
    Initial(gA(ci_X2=None), Parameter('gA_0',1))
    Initial(mA(), Parameter('mA_0'))
    Initial(pA(gWng=None), Parameter('pA_0'))
    Initial(IWG_oligomer(enzyme=None), Parameter('IWG_olig_0'))
    Initial(enzyme(IWG_oligomer=None), Parameter('enzyme_0',1))
    Initial(pY(gWng=None), Parameter('pY_0'))
    Initial(gWng(pA=None, pY=None), Parameter('gWng_0',1))
    Initial(wng(), Parameter('wng_0'))

    # for ic in model.initial_conditions:
    #    print ic

    Observable('pCI_tot', pCI())
    Observable('pCI_free', pCI(ci_l=None, ci_r=None, X2=None, gA=None))
    Observable('pCI_X2_dimer', pCI(ci_l=None, ci_r=None, X2=0, gA=None) % pX2(x2_l=None, x2_r=None, CI=0))
    Observable('pCI_X2_free', CI_X2_oligomer(gA=None))
    Observable('pCI_X2_bound', CI_X2_oligomer(gA=ANY))
    Observable('pCN_tot', pCN())
    Observable('pCN_free', pCN(cn_l=None, cn_r=None, gX2=None))
    Observable('pX2_tot', pX2())
    Observable('pX2_free', pX2(x2_l=None, x2_r=None, CI=None))
    Observable('pA_tot', pA())
    Observable('pA_free', pA(gWng=None))
    Observable('pIWG_tot', pIWG())
    Observable('pIWG_free', pIWG(iwg_l=None, iwg_r=None, enzyme=None))
    Observable('pIWG_olig_free', IWG_oligomer(enzyme=None))
    Observable('pIWG_olig_bound', IWG_oligomer(enzyme=ANY))
    Observable('pY_tot', pY())
    Observable('pY_free', pY(gWng=None))
    Observable('gWng_boundby_A', gWng(pA=ANY, pY=None))
    Observable('gWng_boundby_Y', gWng(pA=None, pY=ANY))
    Observable('gWng_boundby_both', gWng(pA=ANY, pY=ANY))
    Observable('gWng_empty', gWng(pA=None, pY=None))
    Observable('gWng_any', gWng())
    Observable('mWng_obs', wng())

    #
    ## pCN oligomerization rules (x + x <> x:x)
    Parameter('kf_cn_olig', 1e8)
    Parameter('kr_cn_olig', 1e9)
    Parameter('kf_cn_x2', 1e8)# * 10 ** (order_CNen - 1))
    Parameter('kr_cn_x2', kf_cn_x2.value * \
              (K_CNwg.value/(H_CN.value*H_ci.value)) ** order_CNwg * \
              (kf_cn_olig.value / kr_cn_olig.value) ** (order_CNwg - 1))
    Parameter('ktr_x2', 1e8)
    Parameter('ktl_x2', 1e9)
    Parameter('kdeg_mx2', 1e9)
    Parameter('kdeg_px2', 1e8)

#    print (kf_cn_x2.value)
#    print (kr_cn_x2.value)
#    K1 = kr_cn_olig.value / kf_cn_olig.value
#    K2 = kr_cn_x2.value / kf_cn_x2.value
#    print (K1 ** (order_CNwg - 1) * K2) ** (1. / order_CNwg)
#    print (K_CNwg.value / (H_CN.value*H_ci.value))


    # pCN oligomerization rules (x + x <> x:x)
    x = pCN(cn_l=None, cn_r=None, gX2=None)
    for i in range(order_CNwg-1):
        reactants = []
        products = [pCN(cn_l=None, cn_r=i, gX2=None)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pCN(cn_l=left, cn_r=right, gX2=None))
            # product pattern
            left = i if j==0 else j-1
            products.append(pCN(cn_l=left, cn_r=right, gX2=None))
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('CN_oligomerize_%d_to_%d' % (i+1, i+2), x + y | z, kf_cn_olig, kr_cn_olig)

    # # Promoter binding rule
    olig_unbound = z.copy()
    olig_tobind = z.copy()
    g_unbound = gX2(CN_oligomer=None)
    g_bound = gX2(CN_oligomer=50)
    j = order_CNwg-1
    olig_bound = [x for x in olig_tobind.monomer_patterns]
    olig_bound[j].site_conditions['gX2']=50
    olig_bound.append(g_bound)
    bound_complex = ComplexPattern(olig_bound,None)
    Rule('CN_binds_X2_promoter', olig_unbound + g_unbound | bound_complex, kf_cn_x2, kr_cn_x2)
    Rule('X2_transcription', gX2(CN_oligomer=None) >> gX2(CN_oligomer=None) + mX2(), ktr_x2)

    # Degradation rule
    Rule('mX2_degradation', mX2() >> None, kdeg_mx2)

    # Translation rule
    Rule('X2_translation', mX2() >> mX2() + pX2(x2_l=None, x2_r=None, CI=None), ktl_x2)

    # Degradation rule
    Rule('pX2_degradation', pX2(x2_l=None, x2_r=None, CI=None) >> None, kdeg_px2)

    #
    ## CI:pX2 oligomerization rule
    Parameter('kf_ci_olig', 1e8)
    Parameter('kr_ci_olig', 1e9)
    Parameter('kf_g_ci', 1e8)
    # kr = kf * (Kappa/H)^n / K1^(sum([n+n-1]))
    Parameter('kr_g_ci', kf_g_ci.value * \
              (K_CIwg.value/(H_CI.value*H_ci.value))**order_CIwg * \
              (kf_ci_olig.value / kr_ci_olig.value)**(sum([order_CIwg + (order_CIwg - 1)])))
    Parameter('ktr_a', 1e8)
    Parameter('ktl_a', 1e9)
    Parameter('kdeg_ma', 1e9)
    Parameter('kdeg_pa', 1e8)

#    print (sum([order_CIwg + (order_CIwg - 1)]))
#    print (kf_g_ci.value)
#    print (kr_g_ci.value)
#    K1 = kr_ci_olig.value / kf_ci_olig.value
#    K2 = kr_g_ci.value / kf_g_ci.value
#    print (K1**((sum([order_CIwg + (order_CIwg - 1)])))*K2)**(1./order_CIwg)
#    print (K_CIwg.value / (H_CI.value*H_ci.value))

    Rule('CI_X2_heterodimerize',pCI(ci_l=None, ci_r=None, gA=None, X2=None) + pX2(x2_l=None, x2_r=None, CI=None)
         | pCI(ci_l=None, ci_r=None, gA=None, X2=0) % pX2(x2_l=None, x2_r=None, CI=0), kf_ci_olig,kr_ci_olig)

    # pCI:pX2 oligomerization rules (x + y <> x:y)
    x = pCI(ci_l=None, ci_r=None, gA=None, X2=50) % pX2(x2_l=None, x2_r=None, CI=50)
    for i in range(order_CIwg-1):
        reactants = []
        products = [pCI(ci_l=None, ci_r=i, gA=None, X2=50) % pX2(x2_l=None, x2_r=None, CI=50)]
        for j in range(i+1):
            x2ci_site = 50 + (j + 1)
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pCI(ci_l=left, ci_r=right, gA=None, X2=x2ci_site) % pX2(x2_l=None, x2_r=None, CI=x2ci_site))
            # product pattern
            left = i if j==0 else j-1
            products.append(pCI(ci_l=left, ci_r=right, gA=None, X2=x2ci_site) % pX2(x2_l=None, x2_r=None, CI=x2ci_site))
            reactant_monomers = [y for z in [k.monomer_patterns for k in reactants] for y in z]
        if i == (order_CIwg)-2:
            product_monomers = [CI_X2_oligomer(gA=None)]
        else:
            product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        y = ComplexPattern(reactant_monomers,None)
        z = ComplexPattern(product_monomers, None)
        Rule('CI_X2_oligomerize_%dmer_to_%dmer' % ((i+1)*2, (i+1)*2+2), x + y | z, kf_ci_olig,kr_ci_olig)

    # Promoter binding rule for Wng
    Rule('CI_X2_heterooligomer_binds_A_promoter', CI_X2_oligomer(gA=None) + gA(ci_X2=None)
         | CI_X2_oligomer(gA=0) % gA(ci_X2=0), kf_g_ci, kr_g_ci)

    ## New version, can use the bound complex outright
    # Transcription rule
    Rule('A_transcription', CI_X2_oligomer(gA=0) % gA(ci_X2=0)
         >> CI_X2_oligomer(gA=0) % gA(ci_X2=0) + mA(), ktr_a)

    # Degradation rule
    Rule('mA_degradation', mA() >> None, kdeg_ma)

    # Translation rule
    Rule('A_translation', mA() >> mA() + pA(gWng=None), ktl_a)

    # Degradation rule
    Rule('pA_degradation', pA(gWng=None) >> None, kdeg_pa)

    #
    ## IWG oligomerization rules (x + x <> x:x)
    ## Plus catalysis leading to Y
    Parameter('kf_iwg_olig', 1e8)
    Parameter('kr_iwg_olig', 1e9)
    Parameter('kf_enz', 1e8)
    Parameter('kcat', 1e7) # think about with regard to alpha_WGwg
    # K^(nu_wg) = k_m               * K_d_olig                 ^(nu_wg - 1)
    #           = (k_r + k_cat)/k_f * (kr_iwg_olig / kf_iwg_olig)^(nu_wg - 1)
    #           |
    #       k_r = k_f               * (kf_iwg_olig / kr_iwg_olig)^(nu_wg - 1) * K^(nu_wg) - k_cat
    # Remember K_Wgwg.value = kappa_WGwg, so the 'real' K_WGwg = kappa_WGwg / [IWG]_o
    # and [IWG]_o = k_tr_IWG * k_deg_IWG * k_tr_wg * k_deg_wg
    #             = 1        * H_IWG.value * 1      * H_wg.value
    Parameter('kr_enz', kf_enz.value * \
              (kf_iwg_olig.value / kr_iwg_olig.value) ** (order_WGwg - 1) * \
              (K_WGwg.value/(H_IWG.value*H_wg.value)) ** order_WGwg \
              - kcat.value)
    Parameter('kdeg_y', 1e5)

#    print (kf_enz.value)
#    print (kr_enz.value)
#    print (kcat.value)
#    K1 = kr_iwg_olig.value / kf_iwg_olig.value
#    K2 = (kr_enz.value + kcat.value) / kf_enz.value
#    print (K1 ** (order_WGwg - 1) * K2) ** (1. / order_WGwg)
#    print (K_WGwg.value / (H_IWG.value*H_wg.value))

    # Testing
    # # pIWG oligomerization rules (x + x <> x:x)
    # # plus catalysis leading to Y
    # x = pIWG(iwg_l=None, iwg_r=None, enzyme=None)
    # for i in range(order_WGwg-1):
    #     reactants = []
    #     products = [pIWG(iwg_l=None, iwg_r=i, enzyme=None)]
    #     for j in range(i+1):
    #         # reactant pattern
    #         left = None if j == 0 else j-1
    #         right = None if j == i else j
    #         reactants.append(pIWG(iwg_l=left, iwg_r=right, enzyme=None))
    #         # product pattern
    #         left = i if j==0 else j-1
    #         products.append(pIWG(iwg_l=left, iwg_r=right, enzyme=None))
    #     if i == (order_WGwg)-2:
    #         print (i)
    #         products = [IWG_oligomer(enzyme=None)]
    #     #else:
    #     #    product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
    #     y = ComplexPattern(reactants, None)
    #     z = ComplexPattern(products, None)
    #     Rule('IWG_oligomerize_%dmer_to_%dmer' % (i+1, i+2), x + y | z, kf_iwg_olig, kr_iwg_olig)
    #
    # Rule('IWG_olig_binds_enz', IWG_oligomer(enzyme=None) + enzyme(IWG_oligomer=None)
    #      | IWG_oligomer(enzyme=0) % enzyme(IWG_oligomer=0), kf_enz, kr_enz)
    #
    # Rule('IWG_olig_enz_generate_Y', IWG_oligomer(enzyme=0) % enzyme(IWG_oligomer=0)
    #      >> IWG_oligomer(enzyme=None) + enzyme(IWG_oligomer=None) + pY(gWng=None), kcat)

    # pCN oligomerization rules (x + x <> x:x)
    x = pIWG(iwg_l=None, iwg_r=None, enzyme=None)
    for i in range(order_WGwg-1):
        reactants = []
        products = [pIWG(iwg_l=None, iwg_r=i, enzyme=None)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pIWG(iwg_l=left, iwg_r=right, enzyme=None))
            # product pattern
            left = i if j==0 else j-1
            products.append(pIWG(iwg_l=left, iwg_r=right, enzyme=None))
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('IWG_oligomerize_%d_to_%d' % (i+1, i+2), x + y | z, kf_iwg_olig, kr_iwg_olig)

    # # Promoter binding rule
    olig_unbound = z.copy()
    olig_tobind = z.copy()
    g_unbound = enzyme(IWG_oligomer=None)
    g_bound = enzyme(IWG_oligomer=50)
    j = order_WGwg-1
    olig_bound = [x for x in olig_tobind.monomer_patterns]
    olig_bound[j].site_conditions['enzyme']=50
    olig_bound.append(g_bound)
    bound_complex = ComplexPattern(olig_bound,None)
    Rule('IWG_binds_enzyme', olig_unbound + g_unbound | bound_complex, kf_enz, kr_enz)
    # Rule('X2_transcription', gX2(CN_oligomer=None) >> gX2(CN_oligomer=None) + mX2(), ktr_x2)
    Rule('IWG_olig_enz_generate_Y', bound_complex
          >> olig_unbound + enzyme(IWG_oligomer=None) + pY(gWng=None), kcat)


    # Degradation rule
    Rule('pY_degradation', pY(gWng=None) >> None, kdeg_y)

    #
    ## Transcription and translation of wng, with A and/or Y binding at enhancers of wng
    Parameter('ktr_wg',1) # set the overall transcription rate

    Parameter('alpha_contribution_wg_Y',alpha_WGwg.value / \
              (alpha_WGwg.value + alpha_CIwg.value))

    Parameter('alpha_contribution_wg_A', alpha_CIwg.value / \
              (alpha_WGwg.value + alpha_CIwg.value))

    Parameter('kf_wg_a', 1e8)
    # alpha_CIwg.value = ('our'alpha / K_dwg_A) * (ktl_a / kdeg_pa) * (ktr_a / kdeg_ma)
    # 'our'alpha must be between 0 and 1
    # alpha_CIwg.value / (ktl_a * ktr_a / kdeg_pa * kdeg_ma)= ('our'alpha / (kr_wg_a / kf_wg_a))
    # (alpha_CIwg.value * kdeg_pa * kdeg_ma / ktl_a * ktr_a) = 'our'alpha*kf_wg_a / kr_wg_a
    # kr_wg_a          = 'our'alpha*kf_wg_a / (alpha_CIwg.value * kdeg_pa * kdeg_ma / ktl_a * ktr_a)
    # kr_wg_a          = 'our'alpha*kf_wg_a*ktl_a*ktr_a / (alpha_CIwg.value * kdeg_pa * kdeg_ma)
    Parameter('kr_wg_a', alpha_contribution_wg_A.value * \
              kf_wg_a.value * ktl_a.value * ktr_a.value / \
              (alpha_CIwg.value * kdeg_pa.value * kdeg_ma.value)
              )

    Parameter('kf_wg_y', 1e8)
    # alpha_WGwg.value = ('our'alpha / K_dwg_Y) * (kcat / kdeg_y)
    # 'our'alpha must be between 0 and 1
    # alpha_WGwg.value / (kcat / kdeg_y) = ('our'alpha / (kr_wg_y / kf_wg_y))
    # (alpha_WGwg.value * kdeg_y / kcat) = 'our'alpha*kf_wg_y / kr_wg_y
    # kr_wg_y          = 'our'alpha*kf_wg_y / (alpha_WGwg.value * kdeg_y / kcat)
    # kr_wg_y          = 'our'alpha*kf_wg_y*kcat / (alpha_WGwg.value * kdeg_y)
    Parameter('kr_wg_y', alpha_contribution_wg_Y.value * \
              kf_wg_y.value * kcat.value / \
              (alpha_WGwg.value * kdeg_y.value)
              )
 

    Parameter('ktr_wg_y', alpha_contribution_wg_Y.value * \
               ktr_wg.value)

    Parameter('ktr_wg_a', alpha_contribution_wg_A.value * \
              ktr_wg.value)


    print ((alpha_contribution_wg_A.value / (kr_wg_a.value / kf_wg_a.value)) *
           (ktl_a.value / kdeg_pa.value) * (ktr_a.value / kdeg_ma.value))
    print (alpha_CIwg.value)

    print ((alpha_contribution_wg_Y.value / (kr_wg_y.value/kf_wg_y.value)) *
           (kcat.value / kdeg_y.value))
    print (alpha_WGwg.value)


    Rule('A_binds_gWng', pA(gWng=None) + gWng(pA=None)
         | pA(gWng=0) % gWng(pA=0), kf_wg_a, kr_wg_a)

    Rule('gWng_transcription_by_A_promotion', pA(gWng=0) % gWng(pA=0)
         >> pA(gWng=0) % gWng(pA=0) + wng(), ktr_wg_a)

    Rule('Y_binds_gWng', pY(gWng=None) + gWng(pY=None)
         | pY(gWng=1) % gWng(pY=1), kf_wg_y, kr_wg_y)

    Rule('gWng_transcription_by_Y_promotion', pY(gWng=1) % gWng(pY=1)
         >> pY(gWng=1) % gWng(pY=1) + wng(), ktr_wg_y)

    Rule('wng_degradation', wng() >> None, H_wg)
