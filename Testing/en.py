from pysb import *
from pysb.core import as_complex_pattern, MonomerPattern, ComplexPattern,\
    ReactionPattern, RuleExpression, as_reaction_pattern
import sys

def create_model(nu_wg, nu_cn):

    Model()

    ########### FROM REDUCED MODEL ##########
    Parameter('H_en', 1.0/56.195484)
    Parameter('K_WGen', 0.016870342)
    Parameter('nu_WGen', nu_wg) #8.195488)
    Parameter('K_CNen', 0.0015239986)
    Parameter('nu_CNen', nu_cn) #5.929026)
    ### Extras needed for comparing reduced to mechanistic
    Parameter('H_wg', 1.0 / 65.74263)
    Parameter('H_WG', 1.0 / 95.70869)     #Parameter('H_IWG', 1.0 / 95.70869), Parameter('H_EWG', 1.0 / 47.52389)
    Parameter('H_EWG', 1.0 / 47.52389)
    Parameter('H_ci', 1.0 / 88.05598)
    Parameter('H_CI', 1.0 / 20.967995)  #Parameter('H_CI', 1.0 / 20.967995), Parameter('H_CN', 1.0 / 17.765957)
    Parameter('H_CN', 1.0 / 17.765957)
    #########################################
    
    ########### MECHANISTIC MODEL ###########
    
    order_WGen = int(round(nu_WGen.value)) #8
#     print order_WGen
#     print

    order_CNen = int(round(nu_CNen.value)) #6
#     print order_CNen
#     print

    # instead of doing x for x in range(whatever), need to get the right number into the binding reaction
    # (so not the monomer but the rule)
    #mon_gEN = Monomer('gEN', ['ewg%d' % i for i in range(order_WGen)])
    Monomer('pCN', ['cn_l', 'cn_r', 'X1']) # binding sites relevant for en only
    Monomer('gX1',['CN_oligomer'])
    Monomer('mX1')
    Monomer('pX1', ['x1_l','x1_r','EWG'])
    #mon_pEWG = Monomer('pEWG', ['ewg_l', 'ewg_r', 'gen'] + ['cn%d' % i for i in range(order_CNen)])
    Monomer('pEWG', ['ewg_l', 'ewg_r', 'gEN', 'X1'])# + ['cn%d' % i for i in range(order_CNen)])
    Monomer('EWG_X1_oligomer',['gEN'])
    Monomer('gEN', ['ewg_X1'])
    Monomer('en')

    Monomer('gEWG')
    Monomer('mEWG')
    Monomer('gCN')
    Monomer('mCN')

    #print (model.monomers)
#     print

    #Initial(MonomerPattern(mon_gEN, {s : None for s in mon_gEN.sites}, None), Parameter('gEN_0', 1))
    Initial(pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None), Parameter('pEWG_0'))
#    Initial(MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None), Parameter('pEWG_0'))
    #Initial(pCN(cn_l=None, cn_r=None, ewg=None), Parameter('pCN_0'))
    Initial(pCN(cn_l=None, cn_r=None, X1=None), Parameter('pCN_0'))
    Initial(gX1(CN_oligomer=None), Parameter('gX1_0',1))
    Initial(mX1(), Parameter('mX1_0'))
    Initial(pX1(x1_l=None, x1_r=None, EWG=None), Parameter('pX1_0'))
    Initial(EWG_X1_oligomer(gEN=None),Parameter('EWG_X1_olig_0'))
    Initial(gEN(ewg_X1=None), Parameter('gEN_0',1))
    Initial(en(), Parameter('en_0'))

    Initial(gEWG,Parameter('gEWG_0',1))
    Initial(gCN,Parameter('gCN_0',1))
    #for ic in model.initial_conditions:
    #    print ic
#     print

    #Observable('mEWG',mEWG())
    #Observable('mCN',mCN())

    Observable('pEWG_tot', pEWG())
    Observable('pEWG_free', pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None))
#     Observable('pEWG_free', MonomerPattern(mon_pEWG, {s : None for s in mon_pEWG.sites}, None))    
    #Observable('pEWG_pCN', pEWG(cn0=ANY)) #doesnt exist
    Observable('pEWG_X1_dimer',pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=0) % pX1(x1_l=None, x1_r=None, EWG=0))
    Observable('pEWG_X1_free', EWG_X1_oligomer(gEN=None))
    Observable('pEWG_X1_bound', EWG_X1_oligomer(gEN=ANY))
    Observable('pCN_tot', pCN())
    Observable('pCN_free', pCN(cn_l=None, cn_r=None, X1=None))
    Observable('pX1_tot', pX1())
    Observable('pX1_free', pX1(x1_l=None, x1_r=None, EWG=None))
    Observable('mEN_obs', en())
    
#     print model.observables
#     print

    # pCN oligomerization rules (x + y <> z)
    Parameter('kf_cn_olig', 1e8)
    Parameter('kr_cn_olig', 1e9)
    Parameter('kf_cn_x1', 1e8)# * 10 ** (order_CNen - 1))
    Parameter('kr_cn_x1', kf_cn_x1.value * \
              (K_CNen.value/(H_CN.value*H_ci.value)) ** order_CNen * \
              (kf_cn_olig.value / kr_cn_olig.value) ** (order_CNen - 1))
    Parameter('ktr_x1', 1e8)
    Parameter('ktl_x1', 1e9)
    Parameter('kdeg_mx1', 1e9)
    Parameter('kdeg_px1', 1e8)

    print (kf_cn_x1.value)
    print (kr_cn_x1.value)
    K1 = kr_cn_olig.value / kf_cn_olig.value
    K2 = kr_cn_x1.value / kf_cn_x1.value
    print (K1 ** (order_CNen - 1) * K2) ** (1. / order_CNen)
    print (K_CNen.value / (H_CN.value*H_ci.value))


    # pCN oligomerization rules (x + y <> z)
    x = pCN(cn_l=None, cn_r=None, X1=None)
    for i in range(order_CNen-1):
        reactants = []
        products = [pCN(cn_l=None, cn_r=i, X1=None)]
        for j in range(i+1):
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pCN(cn_l=left, cn_r=right, X1=None))
            # product pattern
            left = i if j==0 else j-1
            products.append(pCN(cn_l=left, cn_r=right, X1=None))
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        Rule('CN_oligomerize_%d_to_%d' % (i+1, i+2), x + y | z, kf_cn_olig, kr_cn_olig)

# #
#     # # Promoter binding rule
#     # olig_unbound = z.copy()
#     # bound_complex = z.copy()
#     # g_unbound = MonomerPattern(mon_gEN, {s: None for s in mon_gEN.sites}, None)
#     # g_bound = MonomerPattern(mon_gEN, {s: None for s in mon_gEN.sites}, None)
#     # j = order_WGen - 1
#     # for mp1, mp2, site in zip(olig_unbound.monomer_patterns,
#     #                           bound_complex.monomer_patterns,
#     #                           g_bound.monomer.sites):
#     #     mp1.site_conditions['gen'] = None
#     #     mp2.site_conditions['gen'] = j
#     #     g_bound.site_conditions[site] = j
#     #     j += 1
#     # bound_complex.monomer_patterns.append(g_bound)
#     # Rule('EWG_binds_promoter', olig_unbound + g_unbound | bound_complex, kf_g_ewg, kr_g_ewg)
#
    # # Promoter binding rule
    olig_unbound = z.copy()
    olig_tobind = z.copy()
    g_unbound = gX1(CN_oligomer=None)
    g_bound = gX1(CN_oligomer=50)
    j = order_CNen-1
    olig_bound = [x for x in olig_tobind.monomer_patterns]
    olig_bound[j].site_conditions['X1']=50
    olig_bound.append(g_bound)
    bound_complex = ComplexPattern(olig_bound,None)
    Rule('CN_binds_X1_promoter', olig_unbound + g_unbound | bound_complex, kf_cn_x1, kr_cn_x1)
    Rule('X1_transcription', gX1(CN_oligomer=None) >> gX1(CN_oligomer=None) + mX1(), ktr_x1)

    # Degradation rule
    Rule('mX1_degradation', mX1() >> None, kdeg_mx1)

    # Translation rule
    Rule('X1_translation', mX1() >> mX1() + pX1(x1_l=None, x1_r=None, EWG=None), ktl_x1)

    # Degradation rule
    Rule('pX1_degradation', pX1(x1_l=None, x1_r=None, EWG=None) >> None, kdeg_px1)

   ## EWG:pX1 oligomerization rule
    ## Based on Leonard's reading of the supplement it seems that an 8-mer of the EWG:X1 complex
    ##  is what binds to gEN in order to get transcription
    ## Currently going to set it up as 8 pairs of EWG:pX1 bound together through EWG
        ## THIS DID NOT WORK, 16 molecules adding up is too much for BNG
    ## Changing it to have pEWG and pX1 form dimers and then add together

    Parameter('kf_ewg_olig', 1e8)
    Parameter('kr_ewg_olig', 1e9)
    Parameter('kf_g_ewg', 1e8)
    # kr = kf * (Kappa/H)^n / K1^(sum([n+n-1]))
    Parameter('kr_g_ewg', kf_g_ewg.value * \
              (K_WGen.value/(H_EWG.value*H_wg.value))**order_WGen * \
#              (K_WGen.value/(1))**order_WGen * \
              (kf_ewg_olig.value / kr_ewg_olig.value)**(sum([order_WGen + (order_WGen - 1)])))
    Parameter('ktr_en', 1)

    print (sum([order_WGen + (order_WGen - 1)]))
    print (kf_g_ewg.value)
    print (kr_g_ewg.value)
    K1 = kr_ewg_olig.value / kf_ewg_olig.value
    K2 = kr_g_ewg.value / kf_g_ewg.value
    print (K1**((sum([order_WGen + (order_WGen - 1)])))*K2)**(1./order_WGen)
    print (K_WGen.value / (H_EWG.value*H_wg.value))

    ## Doesn't currently work
    # reactants = []
    # products = []
    # for j in range(order_WGen):
    #     reactants.append(pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None))
    #     reactants.append(pX1(x1_l=None, x1_r=None, EWG=None))
    #     left = None if j == 0 else j - 1
    #     right = None if j == order_WGen - 1 else j
    #     x1bind = j+50
    #     products.append(pEWG(ewg_l=left, ewg_r=right, gEN=None, X1=x1bind))
    #     products.append(pX1(x1_l=None, x1_r=None, EWG=x1bind))
    # y = ReactionPattern([ComplexPattern([x], None) for x in reactants])
    # z = ComplexPattern(products, None)
    # #Rule('EWG_pX1_oligomerize', y | z, kf_ewg_olig, kr_ewg_olig)
    # Rule('EWG_pX1_oligomerize',  pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None)
    #      + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) + pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None) |
    #      PH(g=None)
    #      #pEWG(ewg_l=None, ewg_r=0, gEN=None, X1=50) % pX1(x1_l=None, x1_r=None, EWG=50) % pEWG(ewg_l=0, ewg_r=1, gEN=None, X1=51) % pX1(x1_l=None, x1_r=None, EWG=51) % pEWG(ewg_l=1, ewg_r=2, gEN=None, X1=52) % pX1(x1_l=None, x1_r=None, EWG=52) % pEWG(ewg_l=2, ewg_r=3, gEN=None, X1=53) % pX1(x1_l=None, x1_r=None, EWG=53) % pEWG(ewg_l=3, ewg_r=4, gEN=None, X1=54) % pX1(x1_l=None, x1_r=None, EWG=54) % pEWG(ewg_l=4, ewg_r=5, gEN=None, X1=55) % pX1(x1_l=None, x1_r=None, EWG=55) % pEWG(ewg_l=5, ewg_r=6, gEN=None, X1=56) % pX1(x1_l=None, x1_r=None, EWG=56) % pEWG(ewg_l=6, ewg_r=None, gEN=None, X1=57) % pX1(x1_l=None, x1_r=None, EWG=57)
    #      ,kf_ewg_olig,kr_ewg_olig)

# With more than 10 species trying to be added together, get a BngInterfaceError (no text)
#    File "/home/sbeik/miniconda3/envs/b2r/lib/python2.7/site-packages/pysb/bng.py", line 672, in generate_network
#       bngfile.execute()
#    File "/home/sbeik/miniconda3/envs/b2r/lib/python2.7/site-packages/pysb/bng.py", line 487, in execute
#       p_err.rstrip())
#pysb.bng.BngInterfaceError:

    Rule('EWG_X1_heterodimerize',pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=None) + pX1(x1_l=None, x1_r=None, EWG=None)
         | pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=0) % pX1(x1_l=None, x1_r=None, EWG=0), kf_ewg_olig,kr_ewg_olig)

    # pEWG:pX1 oligomerization rules (x + y <> z)
    x = pEWG(ewg_l=None, ewg_r=None, gEN=None, X1=50) % pX1(x1_l=None, x1_r=None, EWG=50)
    for i in range(order_WGen-1):
        reactants = []
        products = [pEWG(ewg_l=None, ewg_r=i, gEN=None, X1=50) % pX1(x1_l=None, x1_r=None, EWG=50)]
        for j in range(i+1):
            x1ewg_site = 50 + (j + 1)
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            reactants.append(pEWG(ewg_l=left, ewg_r=right, gEN=None, X1=x1ewg_site) % pX1(x1_l=None, x1_r=None, EWG=x1ewg_site))
            # product pattern
            left = i if j==0 else j-1
            products.append(pEWG(ewg_l=left, ewg_r=right, gEN=None, X1=x1ewg_site) % pX1(x1_l=None, x1_r=None, EWG=x1ewg_site))
            reactant_monomers = [y for z in [k.monomer_patterns for k in reactants] for y in z]
        if i == (order_WGen)-2:
            product_monomers = [EWG_X1_oligomer(gEN=None)]
        else:
            product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        y = ComplexPattern(reactant_monomers,None)
        z = ComplexPattern(product_monomers, None)
        Rule('EWG_X1_oligomerize_%dmer_to_%dmer' % ((i+1)*2, (i+1)*2+2), x + y | z, kf_ewg_olig,kr_ewg_olig)

    # ## Old version using the 16mer (as if the if i == (order_WGen-2) statement above isnt there)
    # # Promoter binding rule for EN now
    # olig_unbound = z.copy()
    # olig_tobind = z.copy()
    # g_unbound = gEN(ewg_X1=None)
    # g_bound = gEN(ewg_X1=100)
    # j = 0
    # olig_bound = [x for x in olig_tobind.monomer_patterns]
    # olig_bound[j].site_conditions['gEN']=100
    # olig_bound.append(g_bound)
    # bound_complex = ComplexPattern(olig_bound, None)
    # Rule('EWG_X1_heterooligomer_binds_EN_promoter', PH(g=None) + gEN(ewg_X1=None) | PH(g=0) % gEN(ewg_X1=0), kf_g_ewg, kr_g_ewg)

    # ## Old version using the output of the promoter-binding rule
    # # Transcription rule
    # #mp = MonomerPattern(mon_gEN, {s : ANY for s in mon_gEN.sites}, None)
    # gene_ready = bound_complex.copy()
    # Rule('EN_transcription', gene_ready >> gene_ready + en(), ktr_en)

    ## New version using species oligomer() as the 16mer
    # Promoter binding rule for EN now
    Rule('EWG_X1_heterooligomer_binds_EN_promoter', EWG_X1_oligomer(gEN=None) + gEN(ewg_X1=None)
         | EWG_X1_oligomer(gEN=0) % gEN(ewg_X1=0), kf_g_ewg, kr_g_ewg)

    ## New version, can use the bound complex outright
    # Transcription rule
    Rule('EN_transcription', EWG_X1_oligomer(gEN=0) % gEN(ewg_X1=0)
         >> EWG_X1_oligomer(gEN=0) % gEN(ewg_X1=0) + en(), ktr_en)

    # Degradation rule
    #Rule('mRNA_degradation', mEN() >> None, H_en)

    # Degradation rule
    Rule('en_degradation', en() >> None, H_en)
    #Rule('en_degradation', en() >> None, Parameter('endegtest',1))
#    print (model.rules)
