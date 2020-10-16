from pysb import *
from pysb.core import as_complex_pattern, MonomerPattern, ComplexPattern, \
    ReactionPattern, RuleExpression, as_reaction_pattern, ComponentDuplicateNameError
import sys
import copy

def homo_oligomerize(homo_olig,order,site_l,site_r,kf,kr):
    #x = copy.deepcopy(homo_olig)
    h_o = ComplexPattern([homo_olig],None)
    gb = h_o.copy()
    x = [k for k in gb.monomer_patterns][0]
    mname = homo_olig.monomer.name
    for i in range(order - 1):
        reactants = []
        #p = homo_olig.copy()
        p = h_o.copy()
        pp = [k for k in p.monomer_patterns][0]
        pp.site_conditions[site_r] = i
        products = [pp]
        for j in range(i + 1):
            left = None if j == 0 else j - 1
            right = None if j == i else j
            r = h_o.copy()
            rr = [k for k in r.monomer_patterns][0]
            rr.site_conditions[site_l] = left
            rr.site_conditions[site_r] = right
            reactants.append(rr)
            # product pattern
            left = i if j == 0 else j - 1
            q = h_o.copy()
            qq = [k for k in q.monomer_patterns][0]
            qq.site_conditions[site_l] = left
            qq.site_conditions[site_r] = right
            products.append(qq)
        y = ComplexPattern(reactants, None)
        z = ComplexPattern(products, None)
        rulename = mname+'_oligomerize_%dmer_to_%dmer'
        try:
            Rule(rulename % (i + 1, i + 2), x + y | z, kf, kr)
        except ComponentDuplicateNameError:
            pass
    return z

def hetero_oligomerize(dimer,order,site0_l,site0_r,site_0to1,site1_l,site1_r,site_1to0,kf,kr):
    x = dimer.copy()
    mname = '_'.join([i.monomer.name for i in x.monomer_patterns])
    for i in range(order-1):
        reactants = []
        p = dimer.copy()
        pr = [j for j in p.monomer_patterns]
        pr[0].site_conditions[site0_r] = i
        pr[0].site_conditions[site_0to1] = 50
        pr[1].site_conditions[site_1to0] = 50
        products = [p]
        for j in range(i+1):
            binding_site = 50 + (j+1)
            # reactant pattern
            left = None if j == 0 else j-1
            right = None if j == i else j
            r = dimer.copy()
            rr = [k for k in r.monomer_patterns]
            rr[0].site_conditions[site0_l] = left
            rr[0].site_conditions[site0_r] = right
            rr[0].site_conditions[site_0to1] = binding_site
            rr[1].site_conditions[site_1to0] = binding_site
            reactants.append(r)
            # product pattern
            left = i if j==0 else j-1
            q = dimer.copy()
            qr = [k for k in q.monomer_patterns]
            qr[0].site_conditions[site0_l] = left
            qr[0].site_conditions[site0_r] = right
            qr[0].site_conditions[site_0to1] = binding_site
            qr[1].site_conditions[site_1to0] = binding_site
            products.append(q)
            reactant_monomers = [y for z in [k.monomer_patterns for k in reactants] for y in z]
        #if i == (order)-2:
        #    product_monomers = [CI_X2_oligomer(e1=None)]
        #else:
        #    product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        product_monomers = [y for z in [k.monomer_patterns for k in products] for y in z]
        y = ComplexPattern(reactant_monomers,None)
        z = ComplexPattern(product_monomers, None)
        rulename = mname+'_oligomerize_%dmer_to_%dmer'
        Rule(rulename % ((i+1)*2, (i+1)*2+2), x + y | z, kf, kr)
    return z

def transcription_inhibition_by_olig_then_tx(inhibitor,inhib_genesite,gene,gene_inhibsite,rna,order,kfbind,krbind,ktr):
    iname = '_'.join([i.monomer.name for i in inhibitor.monomer_patterns])
    if len(set(iname.split('_'))) == 1:
        iname = list(set(iname.split('_')))[0]
    gname = gene.monomer.name
    olig_unbound = inhibitor.copy()
    olig_tobind = inhibitor.copy()
    #g_unbound = gX2(CN_oligomer=None)
    g = ComplexPattern([gene],None)
    g_unbound = g.copy()
    #g_unbound = gene.monomer.copy()
    #g_bound = gX2(CN_oligomer=50)
    g_bound = g.copy()
    gb = [x for x in g_bound.monomer_patterns]
    gb[0].site_conditions[gene_inhibsite] = 50
    j = order-1
    olig_bound = [x for x in olig_tobind.monomer_patterns]
    olig_bound[j].site_conditions[inhib_genesite]=50
    olig_bound.append(gb[0])
    bound_complex = ComplexPattern(olig_bound,None)
    bindrulename = iname+'_binds_'+gname+'_promoter'
    Rule(bindrulename, olig_unbound + g_unbound | bound_complex, kfbind, krbind)
    txrulename = gname+'_transcription'
    Rule(txrulename, g_unbound >> g_unbound + rna, ktr)


Model()

# en
Parameter('H_en', 1.0 / 56.195484)
Parameter('K_WGen', 0.016870342)
Parameter('nu_WGen', 8)  # 8.195488)
#Parameter('nu_WGen', nu_wg)  # 8.195488)
Parameter('K_CNen', 0.0015239986)
Parameter('nu_CNen', 6)  # 5.929026)
#Parameter('nu_CNen', nu_cn)  # 5.929026)
# wg
Parameter('H_wg', 1.0 / 65.74263)
Parameter('K_WGwg', 0.010839762)
Parameter('nu_WGwg', 5) #5.2690034)
#Parameter('nu_WGwg', nu_wg) #5.2690034)
Parameter('alpha_WGwg', 1.2893039)
Parameter('K_CIwg', 0.1358894)
Parameter('nu_CIwg', 3) #3.0448556)
#Parameter('nu_CIwg', nu_ci) #3.0448556)
Parameter('alpha_CIwg', 5.354691)
Parameter('K_CNwg', 0.08410451)
Parameter('nu_CNwg', 2) #2.2106833)
#Parameter('nu_CNwg', nu_cn) #2.2106833)
# ptc
Parameter('K_CIptc', 0.0012549695)
Parameter('nu_CIptc', 3.) # should it be 4??
#Parameter('nu_CIptc', 3.5400295)
Parameter('H_ptc', 1.0/62.299126)
Parameter('K_CNptc', 0.15335482)
Parameter('nu_CNptc', 9)
#Parameter('nu_CNptc', 9.436314)
# ci
Parameter('K_Bci', 0.016921775)
Parameter('nu_Bci', 1.0)
Parameter('H_ci', 1.0/88.05598)
Parameter('K_ENci', 0.05230557)
Parameter('nu_ENci', 3) #3.1201758)
Parameter('B', 1) # determines ci expression

### Extras needed for comparing reduced to mechanistic
Parameter('H_IWG', 1.0 / 95.70869)  # Parameter('H_IWG', 1.0 / 95.70869), Parameter('H_EWG', 1.0 / 47.52389)
Parameter('H_EWG', 1.0 / 47.52389)
Parameter('H_CI', 1.0 / 20.967995)  # Parameter('H_CI', 1.0 / 20.967995), Parameter('H_CN', 1.0 / 17.765957)
Parameter('H_CN', 1.0 / 17.765957)
Parameter('H_EN', 1.0/26.655169)
#########################################

########### MECHANISTIC MODEL ###########

order_WGen = int(8)  # 8
order_CNen = int(6)  # 6
order_WGwg = int(5)  # 5
order_CIwg = int(3)  # 3
order_CNwg = int(2)  # 2
order_CIptc = int(3)  # 3
order_CNptc = int(9)  # 9
order_ENci = int(3)  # 3
order_Bci = int(1.0) # given as 1

#Monomer('pCI', ['ci_l', 'ci_r', 'X2', 'e1', 'X3', 'gPTC'])  # binding sites relevant for wg,ptc only
Monomer('pCI', ['ci_l', 'ci_r', 'protein', 'gene'])  # binding sites relevant for wg,ptc only
#Monomer('pCN', ['cn_l', 'cn_r', 'gX1', 'gX2', 'gX3'])  # binding sites relevant for en,wg,ptc
Monomer('pCN', ['cn_l', 'cn_r', 'gene'])  # binding sites relevant for en,wg,ptc
Monomer('gX1', ['CN_oligomer'])
Monomer('mX1')
Monomer('pX1', ['x1_l', 'x1_r', 'EWG'])
Monomer('gX2', ['CN_oligomer'])
Monomer('mX2')
Monomer('pX2', ['x2_l', 'x2_r', 'CI'])
Monomer('gX3', ['CN_oligomer'])
Monomer('mX3')
Monomer('pX3', ['x3_l', 'x3_r', 'CI'])
Monomer('gX4', ['EN_oligomer'])
Monomer('mX4')
Monomer('pX4', ['x4_l', 'x4_r', 'B'])
#Monomer('pIWG',['iwg_l','iwg_r','e2'])
Monomer('pIWG',['iwg_l', 'iwg_r', 'gene'])
#Monomer('pEWG',['ewg_l','ewg_r','gEN','X1'])
Monomer('pEWG', ['ewg_l', 'ewg_r', 'protein', 'gene'])  # + ['cn%d' % i for i in range(order_CNen)])
Monomer('gEN', ['ewg_X1'])
Monomer('en')
Monomer('pEN',['en_l', 'en_r', 'gene'])
Monomer('pB', ['b_l','b_r','protein','gene'])
Monomer('mAB')
Monomer('pAB',['gWng'])
Monomer('e1', ['CI_X2_oligomer'])
Monomer('e2', ['IWG_oligomer'])
Monomer('gWng', ['pAB'])
Monomer('wng')
Monomer('gPTC', ['CI_X3'])
Monomer('ptc')
Monomer('gCI', ['X4'])
Monomer('ci')

Initial(pCI(ci_l=None, ci_r=None, protein=None, gene=None), Parameter('pCI_0'))
Initial(pCN(cn_l=None, cn_r=None, gene=None), Parameter('pCN_0'))
Initial(pIWG(iwg_l=None, iwg_r=None, gene=None), Parameter('pIWG_0'))
Initial(pEWG(ewg_l=None, ewg_r=None, gene=None, protein=None), Parameter('pEWG_0'))
Initial(pB(b_l=None, b_r=None, protein=None, gene=None), B)

Initial(gX1(CN_oligomer=None), Parameter('gX1_0', 1))
Initial(mX1(), Parameter('mX1_0'))
Initial(pX1(x1_l=None, x1_r=None, EWG=None), Parameter('pX1_0'))
Initial(gX2(CN_oligomer=None), Parameter('gX2_0',1))
Initial(mX2(), Parameter('mX2_0'))
Initial(pX2(x2_l=None, x2_r=None, CI=None), Parameter('pX2_0'))
Initial(gX3(CN_oligomer=None), Parameter('gX3_0',1))
Initial(mX3(), Parameter('mX3_0'))
Initial(pX3(x3_l=None, x3_r=None, CI=None), Parameter('pX3_0'))
Initial(gX4(EN_oligomer=None), Parameter('gX4_0',1))
Initial(mX4(), Parameter('mX4_0'))
Initial(pX4(x4_l=None, x4_r=None, B=None), Parameter('pX4_0'))

Initial(e1(CI_X2_oligomer=None), Parameter('e1_0',1))
Initial(e2(IWG_oligomer=None), Parameter('e2_0',1))
Initial(mAB(), Parameter('mAB_0'))
Initial(pAB(gWng=None), Parameter('AB_0'))

Initial(gEN(ewg_X1=None), Parameter('gEN_0',1))
Initial(en(), Parameter('en_0'))
Initial(pEN(en_l=None, en_r=None, gene=None), Parameter('EN_0'))

Initial(gWng(pAB=None), Parameter('gWng_0',1))
Initial(wng(), Parameter('wng_0'))

Initial(gPTC(CI_X3=None), Parameter('gPTC_0',1))
Initial(ptc(), Parameter('ptc_0'))

Initial(gCI(X4=None), Parameter('gCI_0',1))
Initial(ci(), Parameter('ci_0'))

#for n,ic in enumerate(model.initial_conditions):
#    print (n)
#    print (ic)

Observable('pCI_tot', pCI())
Observable('pCI_free', pCI(ci_l=None, ci_r=None, protein=None, gene=None))
Observable('pCI_X2_dimer', pCI(ci_l=None, ci_r=None, protein=0) % pX2(x2_l=None, x2_r=None, CI=0))
Observable('pCN_tot', pCN())
Observable('pCN_free', pCN(cn_l=None, cn_r=None, gene=None))

Observable('pX1_tot', pX1())
Observable('pX1_free', pX1(x1_l=None, x1_r=None, EWG=None))
Observable('pX2_tot', pX2())
Observable('pX2_free', pX2(x2_l=None, x2_r=None, CI=None))
Observable('pX3_tot', pX3())
Observable('pX3_free', pX3(x3_l=None, x3_r=None, CI=None))
Observable('pX4_tot', pX4())
Observable('pX4_free', pX4(x4_l=None, x4_r=None, B=None))

Observable('pIWG_tot', pIWG())
Observable('pIWG_free', pIWG(iwg_l=None, iwg_r=None, gene=None))
Observable('pEWG_tot', pEWG())
Observable('pEWG_free', pEWG(ewg_l=None, ewg_r=None, gene=None, protein=None))

Observable('mEN_obs', en())
Observable('pEN_tot', pEN())
Observable('pEN_free', pEN(en_l=None, en_r=None, gene=None))
Observable('mWng_obs', wng())
Observable('mPTC_obs', ptc())
Observable('mCI_obs', ci())
Observable('B_free',pB(b_l=None, b_r=None, protein=None, gene=None))


# Wingless

Parameter('kf_cn_olig', 1e8) # same for engrailed, wingless, patched
Parameter('kr_cn_olig', 1e9) #

cn_olig_forx2 = homo_oligomerize(pCN(cn_l=None,cn_r=None,gene=None),order_CNwg,'cn_l','cn_r',kf_cn_olig,kr_cn_olig)

Parameter('kf_cn_x2', 1e8)  # * 10 ** (order_CNen - 1))
Parameter('kr_cn_x2', kf_cn_x2.value * \
          (K_CNwg.value / (H_CN.value * H_ci.value)) ** order_CNwg * \
          (kf_cn_olig.value / kr_cn_olig.value) ** (order_CNwg - 1))
Parameter('ktr_x2', 1e8)

transcription_inhibition_by_olig_then_tx(cn_olig_forx2,'gene',gX2(CN_oligomer=None),'CN_oligomer',mX2(),order_CNwg,kf_cn_x2,kr_cn_x2,ktr_x2)

Parameter('ktl_x2', 1e9)
Parameter('kdeg_mx2', 1e9)
Parameter('kdeg_px2', 1e8)

Rule('mX2_degradation', mX2() >> None, kdeg_mx2)
Rule('X2_translation', mX2() >> mX2() + pX2(x2_l=None, x2_r=None, CI=None), ktl_x2)
Rule('pX2_degradation', pX2(x2_l=None, x2_r=None, CI=None) >> None, kdeg_px2)


Parameter('kf_ci_olig', 1e8)
Parameter('kr_ci_olig', 1e9)


Rule('CI_X2_heterodimerize', pCI(ci_l=None, ci_r=None, protein=None, gene=None) + pX2(x2_l=None, x2_r=None, CI=None)
     | pCI(ci_l=None, ci_r=None, protein=0, gene=None) % pX2(x2_l=None, x2_r=None, CI=0), kf_ci_olig, kr_ci_olig)

ci_x2_olig = hetero_oligomerize(pCI(ci_l=None, ci_r=None, protein=50, gene=None) % pX2(x2_l=None, x2_r=None, CI=50),order_CIwg,'ci_l','ci_r','protein','x2_l','x2_r','CI',kf_ci_olig,kr_ci_olig)

Parameter('kf_iwg_olig', 1e8)
Parameter('kr_iwg_olig', 1e9)

iwg_olig = homo_oligomerize(pIWG(iwg_l=None,iwg_r=None,gene=None),order_WGwg,'iwg_l','iwg_r',kf_iwg_olig,kr_iwg_olig)

Parameter('ktr_wg', 1)  # set the overall transcription rate

Parameter('kf_enh_1', 1e8)
Parameter('kr_enh_1', kf_enh_1.value * \
          (K_CIwg.value / (H_CI.value * H_ci.value)) ** order_CIwg * \
          (kf_ci_olig.value / kr_ci_olig.value) ** (sum([order_CIwg + (order_CIwg - 1)])))
Parameter('kf_enh_2', 1e8)
Parameter('kr_enh_2', kf_enh_2.value * \
          (K_WGwg.value / (H_IWG.value * H_wg.value)) ** order_WGwg * \
          (kf_iwg_olig.value / kr_iwg_olig.value) ** (order_WGwg - 1))
#
# alpha_1, alpha_2
Parameter('kdeg_mab', 1e9)
Parameter('ktl_ab', 1e8)
Parameter('kdeg_pab', 1e9)
Parameter('kf_ab', 1e9)
Parameter('kr_ab', 1e8)

Parameter('alpha_1_mass_axn',
          alpha_CIwg.value * (kr_ab.value * \
                              kdeg_mab.value * kdeg_pab.value) \
          / (kf_ab.value * ktl_ab.value)
          )

Parameter('alpha_2_mass_axn',
          alpha_WGwg.value * (kr_ab.value * \
                              kdeg_mab.value * kdeg_pab.value) \
          / (kf_ab.value * ktl_ab.value)
          )

ci_x2_olig_b = ci_x2_olig.copy()
ci_x2_olig_bound = [k for k in ci_x2_olig_b.monomer_patterns]
ci_x2_olig_bound[0].site_conditions['gene'] = 100 #actually an enhancer
ci_x2_olig_bound.append(e1(CI_X2_oligomer=100))
e1_bound_complex = ComplexPattern(ci_x2_olig_bound, None)
Rule('CI_X2_olig_binds_e1', ci_x2_olig + e1(CI_X2_oligomer=None)
    | e1_bound_complex, kf_enh_1, kr_enh_1)

iwg_olig_b = iwg_olig.copy()
iwg_olig_bound = [k for k in iwg_olig_b.monomer_patterns]
iwg_olig_bound[0].site_conditions['gene'] = 100 # actually an enhancer
iwg_olig_bound.append(e2(IWG_oligomer=100))
e2_bound_complex = ComplexPattern(iwg_olig_bound, None)
Rule('IWG_olig_binds_e2', iwg_olig + e2(IWG_oligomer=None)
     | e2_bound_complex, kf_enh_2, kr_enh_2)

Rule('mAB_production_by_e1', e1_bound_complex
     >> e1_bound_complex + mAB(), alpha_1_mass_axn)

Rule('mAB_production_by_e2', e2_bound_complex
     >> e2_bound_complex + mAB(), alpha_2_mass_axn)

Rule('mAB_degrades', mAB() >> None, kdeg_mab)

Rule('pAB_translation', mAB() >> mAB() + pAB(gWng=None), ktl_ab)

Rule('pAB_degrades', pAB(gWng=None) >> None, kdeg_pab)

Rule('B_binds_gWng', pAB(gWng=None) + gWng(pAB=None)
     | pAB(gWng=0) % gWng(pAB=0), kf_ab, kr_ab)

Rule('wg_transcription', pAB(gWng=0) % gWng(pAB=0)
     >> pAB(gWng=0) % gWng(pAB=0) + wng(), ktr_wg)

Rule('wng_degradation', wng() >> None, H_wg)

# Engrailed

# kf_cn_olig & kr same for wg and en
cn_olig_forx1 = homo_oligomerize(pCN(cn_l=None,cn_r=None,gene=None),order_CNen,'cn_l','cn_r',kf_cn_olig,kr_cn_olig)

Parameter('kf_cn_x1', 1e8)  # * 10 ** (order_CNen - 1))
Parameter('kr_cn_x1', kf_cn_x1.value * \
          (K_CNen.value / (H_CN.value * H_ci.value)) ** order_CNen * \
          (kf_cn_olig.value / kr_cn_olig.value) ** (order_CNen - 1))
Parameter('ktr_x1', 1e8)

transcription_inhibition_by_olig_then_tx(cn_olig_forx1,'gene',gX1(CN_oligomer=None),'CN_oligomer',mX1(),order_CNen,kf_cn_x1,kr_cn_x1,ktr_x1)

Parameter('ktl_x1', 1e9)
Parameter('kdeg_mx1', 1e9)
Parameter('kdeg_px1', 1e8)

Rule('mX1_degradation', mX1() >> None, kdeg_mx1)
Rule('X1_translation', mX1() >> mX1() + pX1(x1_l=None, x1_r=None, EWG=None), ktl_x1)
Rule('pX1_degradation', pX1(x1_l=None, x1_r=None, EWG=None) >> None, kdeg_px1)

Parameter('kf_ewg_olig', 1e8)
Parameter('kr_ewg_olig', 1e9)

Rule('EWG_X1_heterodimerize', pEWG(ewg_l=None, ewg_r=None, protein=None, gene=None) + pX1(x1_l=None, x1_r=None, EWG=None)
     | pEWG(ewg_l=None, ewg_r=None, gene=None, protein=0) % pX1(x1_l=None, x1_r=None, EWG=0), kf_ewg_olig, kr_ewg_olig)

ewg_x1_olig = hetero_oligomerize(pEWG(ewg_l=None, ewg_r=None, gene=None, protein=50) % pX1(x1_l=None, x1_r=None, EWG=50),order_WGen,'ewg_l','ewg_r','protein','x1_l','x1_r','EWG',kf_ewg_olig,kr_ewg_olig)

Parameter('kf_g_ewg', 1e8)
# kr = kf * (Kappa/H)^n / K1^(sum([n+n-1]))
Parameter('kr_g_ewg', kf_g_ewg.value * \
          (K_WGen.value/(H_EWG.value*H_wg.value))**order_WGen * \
#          (K_WGen.value/(1))**order_WGen * \
          (kf_ewg_olig.value / kr_ewg_olig.value)**(sum([order_WGen + (order_WGen - 1)])))
Parameter('ktr_en', 1)

ewg_x1_olig_b = ewg_x1_olig.copy()
ewg_x1_olig_bound = [k for k in ewg_x1_olig_b.monomer_patterns]
ewg_x1_olig_bound[0].site_conditions['gene'] = 100
ewg_x1_olig_bound.append(gEN(ewg_X1=100))
gEN_bound_complex = ComplexPattern(ewg_x1_olig_bound, None)
Rule('EWG_X1_olig_binds_EN_promoter', ewg_x1_olig + gEN(ewg_X1=None)
    | gEN_bound_complex, kf_g_ewg, kr_g_ewg)

Rule('EN_transcription', gEN_bound_complex
     >> gEN_bound_complex + en(), ktr_en)

Rule('en_degradation', en() >> None, H_en)

Rule('EN_translation', en() >> en() + pEN(en_l=None,en_r=None,gene=None), Parameter('ktl_en',5000))

Parameter('H_EN_altered',H_EN.value*5000)
#Rule('EN_degradation', pEN(en_l=None,en_r=None,gene=None) >> None, H_EN)
Rule('EN_degradation', pEN(en_l=None,en_r=None,gene=None) >> None, H_EN_altered)

# # Patched
#
# # kf_cn_olig & kr same for wg, en, ptc
# cn_olig_forx3 = homo_oligomerize(pCN(cn_l=None,cn_r=None,gene=None),order_CNptc,'cn_l','cn_r',kf_cn_olig,kr_cn_olig)
#
# Parameter('kf_cn_x3', 1e8)  # * 10 ** (order_CNen - 1))
# Parameter('kr_cn_x3', kf_cn_x3.value * \
#           (K_CNptc.value / (H_CN.value * H_ci.value)) ** order_CNptc * \
#           (kf_cn_olig.value / kr_cn_olig.value) ** (order_CNptc - 1))
# Parameter('ktr_x3', 1e8)
#
# transcription_inhibition_by_olig_then_tx(cn_olig_forx3,'gene',gX3(CN_oligomer=None),'CN_oligomer',mX3(),order_CNptc,kf_cn_x3,kr_cn_x3,ktr_x3)
#
# Parameter('ktl_x3', 1e9)
# Parameter('kdeg_mx3', 1e9)
# Parameter('kdeg_px3', 1e8)
#
# Rule('mX3_degradation', mX3() >> None, kdeg_mx3)
# Rule('X3_translation', mX3() >> mX3() + pX3(x3_l=None, x3_r=None, CI=None), ktl_x3)
# Rule('pX3_degradation', pX3(x3_l=None, x3_r=None, CI=None) >> None, kdeg_px3)
#
# # specified above and going to be the same
# #Parameter('kf_ci_olig', 1e8)
# #Parameter('kr_ci_olig', 1e9)
#
# Rule('CI_X3_heterodimerize', pCI(ci_l=None, ci_r=None, protein=None, gene=None) + pX3(x3_l=None, x3_r=None, CI=None)
#      | pCI(ci_l=None, ci_r=None, protein=0, gene=None) % pX3(x3_l=None, x3_r=None, CI=0), kf_ci_olig, kr_ci_olig)
#
# ci_x3_olig = hetero_oligomerize(pCI(ci_l=None, ci_r=None, gene=None, protein=50) % pX3(x3_l=None, x3_r=None, CI=50),order_CIptc,'ci_l','ci_r','protein','x3_l','x3_r','CI',kf_ci_olig,kr_ci_olig)
#
# Parameter('kf_g_ci', 1e8)
# # kr = kf * (Kappa/H)^n / K1^(sum([n+n-1]))
# Parameter('kr_g_ci', kf_g_ci.value * \
#           (K_CIptc.value/(H_CI.value*H_ci.value))**order_CIptc * \
#           (kf_ci_olig.value / kr_ci_olig.value)**(sum([order_CIptc + (order_CIptc - 1)])))
# Parameter('ktr_ptc', 1)
#
# ci_x3_olig_b = ci_x3_olig.copy()
# ci_x3_olig_bound = [k for k in ci_x3_olig_b.monomer_patterns]
# ci_x3_olig_bound[0].site_conditions['gene'] = 100
# ci_x3_olig_bound.append(gPTC(CI_X3=100))
# gPTC_bound_complex = ComplexPattern(ci_x3_olig_bound, None)
# Rule('CI_X3_olig_binds_PTC_promoter', ci_x3_olig + gPTC(CI_X3=None)
#      | gPTC_bound_complex, kf_g_ci, kr_g_ci)
#
# Rule('PTC_transcription', gPTC_bound_complex
#      >> gPTC_bound_complex + ptc(), ktr_ptc)
#
# Rule('ptc_degradation', ptc() >> None, H_ptc)
#
# #Cubitus Interruptus
#
# Parameter('kf_en_olig', 1e8) #
# Parameter('kr_en_olig', 1e9) #
#
# en_olig_forx4 = homo_oligomerize(pEN(en_l=None,en_r=None,gene=None),order_ENci,'en_l','en_r',kf_en_olig,kr_en_olig)
#
# Parameter('kf_en_x4', 1e8)  # * 10 ** (order_CNen - 1))
# Parameter('kr_en_x4', kf_en_x4.value * \
#           (K_ENci.value / (H_EN.value * H_en.value)) ** order_ENci * \
#           (kf_en_olig.value / kr_en_olig.value) ** (order_ENci - 1))
# Parameter('ktr_x4', 1e8)
#
# print (kf_en_x4.value)
# print (kr_en_x4.value)
# K1 = kr_en_olig.value / kf_en_olig.value
# K2 = kr_en_x4.value / kf_en_x4.value
# print (K1 ** (order_ENci - 1) * K2) ** (1. / order_ENci)
# print (K_ENci.value / (H_EN.value * H_en.value))
#
# transcription_inhibition_by_olig_then_tx(en_olig_forx4,'gene',gX4(EN_oligomer=None),'EN_oligomer',mX4(),order_ENci,kf_en_x4,kr_en_x4,ktr_x4)
#
# Parameter('ktl_x4', 1e9)
# Parameter('kdeg_mx4', 1e9)
# Parameter('kdeg_px4', 1e8)
#
# Rule('mX4_degradation', mX4() >> None, kdeg_mx4)
# Rule('X4_translation', mX4() >> mX4() + pX4(x4_l=None, x4_r=None, B=None), ktl_x4)
# Rule('pX4_degradation', pX4(x4_l=None, x4_r=None, B=None) >> None, kdeg_px4)
#
# Parameter('kf_b_x4', 1e8)
# Parameter('kr_b_x4', 1e9)
# Parameter('kf_g_b', 1e9)
# # # Kappa = (kr_g/kf_g)*(kr_b_x4/kf_b_x4)
# # # Kappa/(kr_b_x4/kf_b_x4) = kr_g/kf_g
# # # Kappa*kf_b_x4/kr_b_x4 = kr_g/kf_g
# # # kr_g = kf_g * Kappa*kf_b_x4/kr_b_x4
# Parameter('kr_g_b', kf_g_b.value * \
# #        kf_b_x4.value * (K_Bci.value/.4) / kr_b_x4.value)
#         kf_b_x4.value * (K_Bci.value*.4) / kr_b_x4.value)
# #          (K_Bci.value / B.value)) #order_Bci = 1
# Parameter('ktr_ci', 1)
#
# print (kf_g_b.value)
# print (kr_g_b.value)
#
# print ((kr_b_x4.value/kf_b_x4.value)*(kr_g_b.value/kf_g_b.value))
# print(K_Bci.value)
#
# Rule('X4_binds_B', pX4(x4_l=None, x4_r=None, B=None) + pB(b_l=None, b_r=None, gene=None, protein=None)
#      | pX4(x4_l=None, x4_r=None, B=0) % pB(b_l=None, b_r=None, gene=None, protein=0), kf_b_x4, kr_b_x4)
#
# Rule('X4_B_binds_CI_promoter', pX4(x4_l=None, x4_r=None, B=0) % pB(b_l=None, b_r=None, gene=None, protein=0) + gCI(X4=None)
#      | pX4(x4_l=None, x4_r=None, B=0) % pB(b_l=None, b_r=None, gene=1, protein=0) % gCI(X4=1), kf_g_b, kr_g_b)
#
# Rule('CI_transcription', pX4(x4_l=None, x4_r=None, B=0) % pB(b_l=None, b_r=None, gene=1, protein=0) % gCI(X4=1)
#      >> pX4(x4_l=None, x4_r=None, B=0) % pB(b_l=None, b_r=None, gene=1, protein=0) % gCI(X4=1) + ci(), ktr_ci)
#
# Rule('ci_degradation', ci() >> None, H_ci)

#for ic in model.initial_conditions:
#    print (ic)

print (model.rules)
