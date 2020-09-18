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
    Monomer('pCI', ['ci_l', 'ci_r', 'X2'])  # binding sites relevant for wg only
    Monomer('pCN', ['cn_l', 'cn_r', 'gX2'])  # binding sites relevant for wg only
    Monomer('gX2', ['CN_oligomer'])
    Monomer('mX2')
    Monomer('pX2', ['x2_l', 'x2_r', 'CI'])
    Monomer('CI_X2_oligomer', ['e1'])
    Monomer('pIWG',['iwg_l', 'iwg_r'])
    Monomer('IWG_oligomer',['e2'])
    Monomer('gWng', ['pB'])
    Monomer('e1', ['CI_X2_oligomer'])
    Monomer('e2', ['IWG_oligomer'])
    Monomer('mB')
    Monomer('pB',['gWng'])
    Monomer('wng')

    # print (model.monomers)

    Initial(pCI(ci_l=None, ci_r=None, X2=None), Parameter('pCI_0'))
    Initial(pCN(cn_l=None, cn_r=None, gX2=None), Parameter('pCN_0'))
    Initial(pIWG(iwg_l=None, iwg_r=None), Parameter('pIWG_0'))
    Initial(gX2(CN_oligomer=None), Parameter('gX2_0',1))
    Initial(mX2(), Parameter('mX2_0'))
    Initial(pX2(x2_l=None, x2_r=None, CI=None), Parameter('pX2_0'))
    Initial(CI_X2_oligomer(e1=None), Parameter('CI_X2_olig_0'))
    Initial(IWG_oligomer(e2=None), Parameter('IWG_olig_0'))
    Initial(gWng(pB=None), Parameter('gWng_0',1))
    Initial(wng(), Parameter('wng_0'))
    Initial(e1(CI_X2_oligomer=None), Parameter('e1_0',1))
    Initial(e2(IWG_oligomer=None), Parameter('e2_0',1))
    Initial(mB(), Parameter('mB_0'))
    Initial(pB(gWng=None), Parameter('B_0'))

    # for ic in model.initial_conditions:
    #    print ic

    Observable('pCI_tot', pCI())
    Observable('pCI_free', pCI(ci_l=None, ci_r=None, X2=None))
    Observable('pCI_X2_dimer', pCI(ci_l=None, ci_r=None, X2=0) % pX2(x2_l=None, x2_r=None, CI=0))
    Observable('pCI_X2_free', CI_X2_oligomer(e1=None))
    Observable('pCI_X2_bound', CI_X2_oligomer(e1=ANY))
    Observable('pCN_tot', pCN())
    Observable('pCN_free', pCN(cn_l=None, cn_r=None, gX2=None))
    Observable('pX2_tot', pX2())
    Observable('pX2_free', pX2(x2_l=None, x2_r=None, CI=None))
    Observable('pIWG_tot', pIWG())
    Observable('pIWG_free', pIWG(iwg_l=None, iwg_r=None))
    Observable('pIWG_olig_free', IWG_oligomer(e2=None))
    Observable('pIWG_olig_bound', IWG_oligomer(e2=ANY))
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
    # shouldnt need this anymore?
    # Parameter('kf_g_ci', 1e8)
    # # kr = kf * (Kappa/H)^n / K1^(sum([n+n-1]))
    # Parameter('kr_g_ci', kf_g_ci.value * \
    #           (K_CIwg.value/(H_CI.value*H_ci.value))**order_CIwg * \
    #           (kf_ci_olig.value / kr_ci_olig.value)**(sum([order_CIwg + (order_CIwg - 1)])))

    Rule('CI_X2_heterodimerize',pCI(ci_l=None, ci_r=None, X2=None) + pX2(x2_l=None, x2_r=None, CI=None)
         | pCI(ci_l=None, ci_r=None, X2=0) % pX2(x2_l=None, x2_r=None, CI=0), kf_ci_olig,kr_ci_olig)

    # pCI:pX2 oligomerization rules (x + y <> x:y)
    x = pCI(ci_l=None, ci_r=None, X2=50) % pX2(x2_l=None, x2_r=None, CI=50)
    for i in range(order_CIwg-1):
        reactants = []
        products = [pCI(ci_l=None, ci_r=i, X2=50) % pX2(x2_l=None, x2_r=None, CI=50)]
        for j in range(i+1):
            x2ci_site = 50 + (j + 1)
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pCI(ci_l=left, ci_r=right, X2=x2ci_site) % pX2(x2_l=None, x2_r=None, CI=x2ci_site))
            # product pattern
            left = i if j==0 else j-1
            products.append(pCI(ci_l=left, ci_r=right, X2=x2ci_site) % pX2(x2_l=None, x2_r=None, CI=x2ci_site))
            reactant_monomers = [y for z in [k.monomer_patterns for k in reactants] for y in z]
        if i == (order_CIwg)-2:
            product_monomers = [CI_X2_oligomer(e1=None)]
        else:
            product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        y = ComplexPattern(reactant_monomers,None)
        z = ComplexPattern(product_monomers, None)
        Rule('CI_X2_oligomerize_%dmer_to_%dmer' % ((i+1)*2, (i+1)*2+2), x + y | z, kf_ci_olig,kr_ci_olig)
    #
    # Now we have the species CI_X2_oligomer()
    #
    ## IWG oligomerization rules (x + x <> x:x)
    ## Plus catalysis leading to Y
    Parameter('kf_iwg_olig', 1e8)
    Parameter('kr_iwg_olig', 1e9)

    # pIWG oligomerization rules (x + x <> x:x)
    x = pIWG(iwg_l=None, iwg_r=None)
    for i in range(order_WGwg-1):
        reactants = []
        products = [pIWG(iwg_l=None, iwg_r=i)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pIWG(iwg_l=left, iwg_r=right))
            # product pattern
            left = i if j==0 else j-1
            products.append(pIWG(iwg_l=left, iwg_r=right))
        if i == (order_WGwg)-2:
            products = [IWG_oligomer(e2=None)]
        #else:
        #    product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('IWG_oligomerize_%dmer_to_%dmer' % (i+1, i+2), x + y | z, kf_iwg_olig, kr_iwg_olig)
    #
    # Now we have the species IWG_oligomer()
    #
    ## Transcription and translation of B, with either oligomer binding at enhancers of B
    Parameter('ktr_wg',1) # set the overall transcription rate

    Parameter('kf_enh_1', 1e8)
    Parameter('kr_enh_1', kf_enh_1.value * \
              (K_CIwg.value/(H_CI.value*H_ci.value))**order_CIwg * \
              (kf_ci_olig.value / kr_ci_olig.value)**(sum([order_CIwg + (order_CIwg - 1)])))
    Parameter('kf_enh_2', 1e8)
    Parameter('kr_enh_2', kf_enh_2.value * \
              (K_WGwg.value/(H_IWG.value*H_wg.value)) ** order_WGwg * \
              (kf_iwg_olig.value / kr_iwg_olig.value) ** (order_WGwg - 1))

    # print (sum([order_CIwg + (order_CIwg - 1)]))
    # print (kf_enh_1.value)
    # print (kr_enh_1.value)
    # K1 = kr_ci_olig.value / kf_ci_olig.value
    # K2 = kr_enh_1.value / kf_enh_1.value
    # print (K1**((sum([order_CIwg + (order_CIwg - 1)])))*K2)**(1./order_CIwg)
    # print (K_CIwg.value / (H_CI.value*H_ci.value))

    # print (order_WGwg,(order_WGwg - 1))
    # print (kf_enh_2.value)
    # print (kr_enh_2.value)
    # K1 = kr_iwg_olig.value / kf_iwg_olig.value
    # K2 = kr_enh_2.value / kf_enh_2.value
    # print (K1**(order_WGwg - 1)*K2)**(1./order_WGwg)
    # print (K_WGwg.value / (H_IWG.value*H_wg.value))

    # alpha_1, alpha_2
    Parameter('kdeg_mb', 1e9)
    Parameter('ktl_b',1e8)
    Parameter('kdeg_pb',1e9)
    Parameter('kf_b',1e9)
    Parameter('kr_b',1e8)

    # alpha_1_mass_axn is the rate at which the bound enhancer_1 results in mB
    # alpha_CIwg = (ktl_b * alpha_1_mass_axn) / ((kr_b / kf_b) * kdeg_pb * kdeg_mb)
    # (ktl_b*alpha_1_mass_axn) = alpha_CI_wg * ((kr_b / kf_b) * kdeg_pb * kdeg_mb)
    # alpha_1_mass_axn = { alpha_CI_wg * ((kr_b / kf_b) * kdeg_pb * kdeg_mb)} / ktl_B
    # alpha_1_mass_axn = { alpha_CI_wg * kr_b * kdeg_pb * kdeg_mb)} / (kf_b * ktl_B)
    Parameter('alpha_1_mass_axn',
              alpha_CIwg.value * (kr_b.value * \
                                 kdeg_mb.value*kdeg_pb.value) \
                               / (kf_b.value * ktl_b.value)
              )

    # alpha_2_mass_axn is the rate at which the bound enhancer_2 results in mB
    # alpha_WGwg = alpha_2_mass_axn / ((kr_b / kf_b) * kdeg_pb * kdeg_mb)
    # alpha_2_mass_axn = alpha_CI_wg * ((kr_b / kf_b) * kdeg_pb * kdeg_mb)
    Parameter('alpha_2_mass_axn',
                            alpha_WGwg.value * (kr_b.value * \
                            kdeg_mb.value * kdeg_pb.value) \
                                / (kf_b.value * ktl_b.value)

    )

    # print (alpha_1_mass_axn.value)
    # print (alpha_1_mass_axn.value * ktl_b.value / ((kr_b.value / kf_b.value) * kdeg_pb.value * kdeg_mb.value))
    # print (alpha_CIwg.value)
    #
    # print (alpha_2_mass_axn.value)
    # print (alpha_2_mass_axn.value * ktl_b.value / ((kr_b.value / kf_b.value) * kdeg_pb.value * kdeg_mb.value))
    # print (alpha_WGwg.value)

    Rule('CI_X2_olig_binds_e1', CI_X2_oligomer(e1=None) + e1(CI_X2_oligomer=None)
         | CI_X2_oligomer(e1=0) % e1(CI_X2_oligomer=0), kf_enh_1, kr_enh_1)

    Rule('IWG_olig_binds_e2', IWG_oligomer(e2=None) + e2(IWG_oligomer=None)
         | IWG_oligomer(e2=1) % e2(IWG_oligomer=1), kf_enh_2, kr_enh_2)

    Rule('mB_production_by_e1', CI_X2_oligomer(e1=0) % e1(CI_X2_oligomer=0)
         >> CI_X2_oligomer(e1=0) % e1(CI_X2_oligomer=0) + mB(), alpha_1_mass_axn)

    Rule('mB_production_by_e2', IWG_oligomer(e2=1) % e2(IWG_oligomer=1)
         >> IWG_oligomer(e2=1) % e2(IWG_oligomer=1) + mB(), alpha_2_mass_axn)

    Rule('mB_degrades', mB() >> None, kdeg_mb)

    Rule('pB_translation', mB() >> mB() + pB(gWng=None), ktl_b)

    Rule('pB_degrades', pB(gWng=None) >> None, kdeg_pb)

    Rule('B_binds_gWng', pB(gWng=None) + gWng(pB=None)
         | pB(gWng=0) % gWng(pB=0), kf_b, kr_b)

    Rule('wg_transcription', pB(gWng=0) % gWng(pB=0)
         >> pB(gWng=0) % gWng(pB=0) + wng(), ktr_wg)

    Rule('wng_degradation', wng() >> None, H_wg)
