
from pysb import *
from SPN_macros import transcription_cooperative
from SPN_macros import transcription_cooperative_n
from SPN_macros import translation
from SPN_macros import perimeter_diffusion
from SPN_macros import face_to_face_diffusion
from SPN_macros import dimerization
from pysb.macros import bind
from SPN_macros import bind_multiple
from SPN_macros import exo_endo
from SPN_macros import degradation

Model()

# Monomers and Initials

Monomer('x1_g_1', ['CN'])
Monomer('x1_g_2', ['CN'])
Monomer('x1_g_3', ['CN'])
Monomer('x1_g_4', ['CN'])

Initial(x1_g_1(CN=None), Parameter('x1_g_1_0', 1))
Initial(x1_g_2(CN=None), Parameter('x1_g_2_0', 1))
Initial(x1_g_3(CN=None), Parameter('x1_g_3_0', 1))
Initial(x1_g_4(CN=None), Parameter('x1_g_4_0', 1))

Monomer('x1_1')
Monomer('x1_2')
Monomer('x1_3')
Monomer('x1_4')

Initial(x1_1(), Parameter('x1_1_0', 1))
Initial(x1_2(), Parameter('x1_2_0', 1))
Initial(x1_3(), Parameter('x1_3_0', 1))
Initial(x1_4(), Parameter('x1_4_0', 1))

Monomer('X1_1', ['EWG'])
Monomer('X1_2', ['EWG'])
Monomer('X1_3', ['EWG'])
Monomer('X1_4', ['EWG'])

Initial(X1_1(EWG=None), Parameter('X1_1_0', 1))
Initial(X1_2(EWG=None), Parameter('X1_2_0', 1))
Initial(X1_3(EWG=None), Parameter('X1_3_0', 1))
Initial(X1_4(EWG=None), Parameter('X1_4_0', 1))

Monomer('en_g_1', ['EWG_X1'])
Monomer('en_g_2', ['EWG_X1'])
Monomer('en_g_3', ['EWG_X1'])
Monomer('en_g_4', ['EWG_X1'])

Initial(en_g_1(EWG_X1=None), Parameter('en_g_1_0', 1))
Initial(en_g_2(EWG_X1=None), Parameter('en_g_2_0', 1))
Initial(en_g_3(EWG_X1=None), Parameter('en_g_3_0', 1))
Initial(en_g_4(EWG_X1=None), Parameter('en_g_4_0', 1))

Monomer('en_1')
Monomer('en_2')
Monomer('en_3')
Monomer('en_4')

Initial(en_1(), Parameter('en_1_0', 1))
Initial(en_2(), Parameter('en_2_0', 1))
Initial(en_3(), Parameter('en_3_0', 1))
Initial(en_4(), Parameter('en_4_0', 1))

Monomer('EN_1', ['x5_g', 'X7', 'hh_g'])
Monomer('EN_2', ['x5_g', 'X7', 'hh_g'])
Monomer('EN_3', ['x5_g', 'X7', 'hh_g'])
Monomer('EN_4', ['x5_g', 'X7', 'hh_g'])

Initial(EN_1(X7=None, x5_g=None, hh_g=None), Parameter('EN_1_0', 1))
Initial(EN_2(X7=None, x5_g=None, hh_g=None), Parameter('EN_2_0', 1))
Initial(EN_3(X7=None, x5_g=None, hh_g=None), Parameter('EN_3_0', 1))
Initial(EN_4(X7=None, x5_g=None, hh_g=None), Parameter('EN_4_0', 1))

Monomer('x2_g_1', ['CN'])
Monomer('x2_g_2', ['CN'])
Monomer('x2_g_3', ['CN'])
Monomer('x2_g_4', ['CN'])

Initial(x2_g_1(CN=None), Parameter('x2_g_1_0', 1))
Initial(x2_g_2(CN=None), Parameter('x2_g_2_0', 1))
Initial(x2_g_3(CN=None), Parameter('x2_g_3_0', 1))
Initial(x2_g_4(CN=None), Parameter('x2_g_4_0', 1))

Monomer('x2_1')
Monomer('x2_2')
Monomer('x2_3')
Monomer('x2_4')

Initial(x2_1(), Parameter('x2_1_0', 1))
Initial(x2_2(), Parameter('x2_2_0', 1))
Initial(x2_3(), Parameter('x2_3_0', 1))
Initial(x2_4(), Parameter('x2_4_0', 1))

Monomer('X2_1', ['wg_g'])
Monomer('X2_2', ['wg_g'])
Monomer('X2_3', ['wg_g'])
Monomer('X2_4', ['wg_g'])

Initial(X2_1(wg_g=None), Parameter('X2_1_0', 1))
Initial(X2_2(wg_g=None), Parameter('X2_2_0', 1))
Initial(X2_3(wg_g=None), Parameter('X2_3_0', 1))
Initial(X2_4(wg_g=None), Parameter('X2_4_0', 1))

Monomer('x3_g_1', ['IWG'])
Monomer('x3_g_2', ['IWG'])
Monomer('x3_g_3', ['IWG'])
Monomer('x3_g_4', ['IWG'])

Initial(x3_g_1(IWG=None), Parameter('x3_g_1_0', 1))
Initial(x3_g_2(IWG=None), Parameter('x3_g_2_0', 1))
Initial(x3_g_3(IWG=None), Parameter('x3_g_3_0', 1))
Initial(x3_g_4(IWG=None), Parameter('x3_g_4_0', 1))

Monomer('x3_1')
Monomer('x3_2')
Monomer('x3_3')
Monomer('x3_4')

Initial(x3_1(), Parameter('x3_1_0', 1))
Initial(x3_2(), Parameter('x3_2_0', 1))
Initial(x3_3(), Parameter('x3_3_0', 1))
Initial(x3_4(), Parameter('x3_4_0', 1))

Monomer('X3_1', ['wg_g'])
Monomer('X3_2', ['wg_g'])
Monomer('X3_3', ['wg_g'])
Monomer('X3_4', ['wg_g'])

Initial(X3_1(wg_g=None), Parameter('X3_1_0', 1))
Initial(X3_2(wg_g=None), Parameter('X3_2_0', 1))
Initial(X3_3(wg_g=None), Parameter('X3_3_0', 1))
Initial(X3_4(wg_g=None), Parameter('X3_4_0', 1))

Monomer('wg_g_1', ['X'])
Monomer('wg_g_2', ['X'])
Monomer('wg_g_3', ['X'])
Monomer('wg_g_4', ['X'])

Initial(wg_g_1(X=None), Parameter('wg_g_1_0', 1))
Initial(wg_g_2(X=None), Parameter('wg_g_2_0', 1))
Initial(wg_g_3(X=None), Parameter('wg_g_3_0', 1))
Initial(wg_g_4(X=None), Parameter('wg_g_4_0', 1))

Monomer('wg_1')
Monomer('wg_2')
Monomer('wg_3')
Monomer('wg_4')

Initial(wg_1(), Parameter('wg_1_0', 1))
Initial(wg_2(), Parameter('wg_2_0', 1))
Initial(wg_3(), Parameter('wg_3_0', 1))
Initial(wg_4(), Parameter('wg_4_0', 1))

Monomer('IWG_1', ['x3_g'])
Monomer('IWG_2', ['x3_g'])
Monomer('IWG_3', ['x3_g'])
Monomer('IWG_4', ['x3_g'])

Initial(IWG_1(x3_g=None), Parameter('IWG_1_0', 1))
Initial(IWG_2(x3_g=None), Parameter('IWG_2_0', 1))
Initial(IWG_3(x3_g=None), Parameter('IWG_3_0', 1))
Initial(IWG_4(x3_g=None), Parameter('IWG_4_0', 1))

Monomer('EWG_1_1', ['X1', 'en_g'])
Monomer('EWG_1_2', ['X1', 'en_g'])
Monomer('EWG_1_3', ['X1', 'en_g'])
Monomer('EWG_1_4', ['X1', 'en_g'])
Monomer('EWG_1_5', ['X1', 'en_g'])
Monomer('EWG_1_6', ['X1', 'en_g'])
Monomer('EWG_2_1', ['X1', 'en_g'])
Monomer('EWG_2_2', ['X1', 'en_g'])
Monomer('EWG_2_3', ['X1', 'en_g'])
Monomer('EWG_2_4', ['X1', 'en_g'])
Monomer('EWG_2_5', ['X1', 'en_g'])
Monomer('EWG_2_6', ['X1', 'en_g'])
Monomer('EWG_3_1', ['X1', 'en_g'])
Monomer('EWG_3_2', ['X1', 'en_g'])
Monomer('EWG_3_3', ['X1', 'en_g'])
Monomer('EWG_3_4', ['X1', 'en_g'])
Monomer('EWG_3_5', ['X1', 'en_g'])
Monomer('EWG_3_6', ['X1', 'en_g'])
Monomer('EWG_4_1', ['X1', 'en_g'])
Monomer('EWG_4_2', ['X1', 'en_g'])
Monomer('EWG_4_3', ['X1', 'en_g'])
Monomer('EWG_4_4', ['X1', 'en_g'])
Monomer('EWG_4_5', ['X1', 'en_g'])
Monomer('EWG_4_6', ['X1', 'en_g'])

Initial(EWG_1_1(X1=None, en_g=None), Parameter('EWG_1_1_0', 1))
Initial(EWG_1_2(X1=None, en_g=None), Parameter('EWG_1_2_0', 1))
Initial(EWG_1_3(X1=None, en_g=None), Parameter('EWG_1_3_0', 1))
Initial(EWG_1_4(X1=None, en_g=None), Parameter('EWG_1_4_0', 1))
Initial(EWG_1_5(X1=None, en_g=None), Parameter('EWG_1_5_0', 1))
Initial(EWG_1_6(X1=None, en_g=None), Parameter('EWG_1_6_0', 1))
Initial(EWG_2_1(X1=None, en_g=None), Parameter('EWG_2_1_0', 1))
Initial(EWG_2_2(X1=None, en_g=None), Parameter('EWG_2_2_0', 1))
Initial(EWG_2_3(X1=None, en_g=None), Parameter('EWG_2_3_0', 1))
Initial(EWG_2_4(X1=None, en_g=None), Parameter('EWG_2_4_0', 1))
Initial(EWG_2_5(X1=None, en_g=None), Parameter('EWG_2_5_0', 1))
Initial(EWG_2_6(X1=None, en_g=None), Parameter('EWG_2_6_0', 1))
Initial(EWG_3_1(X1=None, en_g=None), Parameter('EWG_3_1_0', 1))
Initial(EWG_3_2(X1=None, en_g=None), Parameter('EWG_3_2_0', 1))
Initial(EWG_3_3(X1=None, en_g=None), Parameter('EWG_3_3_0', 1))
Initial(EWG_3_4(X1=None, en_g=None), Parameter('EWG_3_4_0', 1))
Initial(EWG_3_5(X1=None, en_g=None), Parameter('EWG_3_5_0', 1))
Initial(EWG_3_6(X1=None, en_g=None), Parameter('EWG_3_6_0', 1))
Initial(EWG_4_1(X1=None, en_g=None), Parameter('EWG_4_1_0', 1))
Initial(EWG_4_2(X1=None, en_g=None), Parameter('EWG_4_2_0', 1))
Initial(EWG_4_3(X1=None, en_g=None), Parameter('EWG_4_3_0', 1))
Initial(EWG_4_4(X1=None, en_g=None), Parameter('EWG_4_4_0', 1))
Initial(EWG_4_5(X1=None, en_g=None), Parameter('EWG_4_5_0', 1))
Initial(EWG_4_6(X1=None, en_g=None), Parameter('EWG_4_6_0', 1))

Monomer('x4_g_1', ['CN'])
Monomer('x4_g_2', ['CN'])
Monomer('x4_g_3', ['CN'])
Monomer('x4_g_4', ['CN'])

Initial(x4_g_1(CN=None), Parameter('x4_g_1_0', 1))
Initial(x4_g_2(CN=None), Parameter('x4_g_2_0', 1))
Initial(x4_g_3(CN=None), Parameter('x4_g_3_0', 1))
Initial(x4_g_4(CN=None), Parameter('x4_g_4_0', 1))

Monomer('x4_1')
Monomer('x4_2')
Monomer('x4_3')
Monomer('x4_4')

Initial(x4_1(), Parameter('x4_1_0', 1))
Initial(x4_2(), Parameter('x4_2_0', 1))
Initial(x4_3(), Parameter('x4_3_0', 1))
Initial(x4_4(), Parameter('x4_4_0', 1))

Monomer('X4_1', ['CI'])
Monomer('X4_2', ['CI'])
Monomer('X4_3', ['CI'])
Monomer('X4_4', ['CI'])

Initial(X4_1(CI=None), Parameter('X4_1_0', 1))
Initial(X4_2(CI=None), Parameter('X4_2_0', 1))
Initial(X4_3(CI=None), Parameter('X4_3_0', 1))
Initial(X4_4(CI=None), Parameter('X4_4_0', 1))

Monomer('ptc_g_1', ['CI_X4'])
Monomer('ptc_g_2', ['CI_X4'])
Monomer('ptc_g_3', ['CI_X4'])
Monomer('ptc_g_4', ['CI_X4'])

Initial(ptc_g_1(CI_X4=None), Parameter('ptc_g_1_0', 1))
Initial(ptc_g_2(CI_X4=None), Parameter('ptc_g_2_0', 1))
Initial(ptc_g_3(CI_X4=None), Parameter('ptc_g_3_0', 1))
Initial(ptc_g_4(CI_X4=None), Parameter('ptc_g_4_0', 1))

Monomer('ptc_1')
Monomer('ptc_2')
Monomer('ptc_3')
Monomer('ptc_4')

Initial(ptc_1(), Parameter('ptc_1_0', 1))
Initial(ptc_2(), Parameter('ptc_2_0', 1))
Initial(ptc_3(), Parameter('ptc_3_0', 1))
Initial(ptc_4(), Parameter('ptc_4_0', 1))

Monomer('PTC_1_1', ['HH', 'x6_g'])
Monomer('PTC_1_2', ['HH', 'x6_g'])
Monomer('PTC_1_3', ['HH', 'x6_g'])
Monomer('PTC_1_4', ['HH', 'x6_g'])
Monomer('PTC_1_5', ['HH', 'x6_g'])
Monomer('PTC_1_6', ['HH', 'x6_g'])
Monomer('PTC_2_1', ['HH', 'x6_g'])
Monomer('PTC_2_2', ['HH', 'x6_g'])
Monomer('PTC_2_3', ['HH', 'x6_g'])
Monomer('PTC_2_4', ['HH', 'x6_g'])
Monomer('PTC_2_5', ['HH', 'x6_g'])
Monomer('PTC_2_6', ['HH', 'x6_g'])
Monomer('PTC_3_1', ['HH', 'x6_g'])
Monomer('PTC_3_2', ['HH', 'x6_g'])
Monomer('PTC_3_3', ['HH', 'x6_g'])
Monomer('PTC_3_4', ['HH', 'x6_g'])
Monomer('PTC_3_5', ['HH', 'x6_g'])
Monomer('PTC_3_6', ['HH', 'x6_g'])
Monomer('PTC_4_1', ['HH', 'x6_g'])
Monomer('PTC_4_2', ['HH', 'x6_g'])
Monomer('PTC_4_3', ['HH', 'x6_g'])
Monomer('PTC_4_4', ['HH', 'x6_g'])
Monomer('PTC_4_5', ['HH', 'x6_g'])
Monomer('PTC_4_6', ['HH', 'x6_g'])

Initial(PTC_1_1(HH=None, x6_g=None), Parameter('PTC_1_1_0', 1))
Initial(PTC_1_2(HH=None, x6_g=None), Parameter('PTC_1_2_0', 1))
Initial(PTC_1_3(HH=None, x6_g=None), Parameter('PTC_1_3_0', 1))
Initial(PTC_1_4(HH=None, x6_g=None), Parameter('PTC_1_4_0', 1))
Initial(PTC_1_5(HH=None, x6_g=None), Parameter('PTC_1_5_0', 1))
Initial(PTC_1_6(HH=None, x6_g=None), Parameter('PTC_1_6_0', 1))
Initial(PTC_2_1(HH=None, x6_g=None), Parameter('PTC_2_1_0', 1))
Initial(PTC_2_2(HH=None, x6_g=None), Parameter('PTC_2_2_0', 1))
Initial(PTC_2_3(HH=None, x6_g=None), Parameter('PTC_2_3_0', 1))
Initial(PTC_2_4(HH=None, x6_g=None), Parameter('PTC_2_4_0', 1))
Initial(PTC_2_5(HH=None, x6_g=None), Parameter('PTC_2_5_0', 1))
Initial(PTC_2_6(HH=None, x6_g=None), Parameter('PTC_2_6_0', 1))
Initial(PTC_3_1(HH=None, x6_g=None), Parameter('PTC_3_1_0', 1))
Initial(PTC_3_2(HH=None, x6_g=None), Parameter('PTC_3_2_0', 1))
Initial(PTC_3_3(HH=None, x6_g=None), Parameter('PTC_3_3_0', 1))
Initial(PTC_3_4(HH=None, x6_g=None), Parameter('PTC_3_4_0', 1))
Initial(PTC_3_5(HH=None, x6_g=None), Parameter('PTC_3_5_0', 1))
Initial(PTC_3_6(HH=None, x6_g=None), Parameter('PTC_3_6_0', 1))
Initial(PTC_4_1(HH=None, x6_g=None), Parameter('PTC_4_1_0', 1))
Initial(PTC_4_2(HH=None, x6_g=None), Parameter('PTC_4_2_0', 1))
Initial(PTC_4_3(HH=None, x6_g=None), Parameter('PTC_4_3_0', 1))
Initial(PTC_4_4(HH=None, x6_g=None), Parameter('PTC_4_4_0', 1))
Initial(PTC_4_5(HH=None, x6_g=None), Parameter('PTC_4_5_0', 1))
Initial(PTC_4_6(HH=None, x6_g=None), Parameter('PTC_4_6_0', 1))

Monomer('x5_g_1', ['EN'])
Monomer('x5_g_2', ['EN'])
Monomer('x5_g_3', ['EN'])
Monomer('x5_g_4', ['EN'])

Initial(x5_g_1(EN=None), Parameter('x5_g_1_0', 1))
Initial(x5_g_2(EN=None), Parameter('x5_g_2_0', 1))
Initial(x5_g_3(EN=None), Parameter('x5_g_3_0', 1))
Initial(x5_g_4(EN=None), Parameter('x5_g_4_0', 1))

Monomer('x5_1')
Monomer('x5_2')
Monomer('x5_3')
Monomer('x5_4')

Initial(x5_1(), Parameter('x5_1_0', 1))
Initial(x5_2(), Parameter('x5_2_0', 1))
Initial(x5_3(), Parameter('x5_3_0', 1))
Initial(x5_4(), Parameter('x5_4_0', 1))

Monomer('X5_1', ['B'])
Monomer('X5_2', ['B'])
Monomer('X5_3', ['B'])
Monomer('X5_4', ['B'])

Initial(X5_1(B=None), Parameter('X5_1_0', 1))
Initial(X5_2(B=None), Parameter('X5_2_0', 1))
Initial(X5_3(B=None), Parameter('X5_3_0', 1))
Initial(X5_4(B=None), Parameter('X5_4_0', 1))

Monomer('B_1', ['X5', 'ci_g'])
Monomer('B_2', ['X5', 'ci_g'])
Monomer('B_3', ['X5', 'ci_g'])
Monomer('B_4', ['X5', 'ci_g'])

Initial(B_1(X5=None, ci_g=None), Parameter('B_1_0', 1))
Initial(B_2(X5=None, ci_g=None), Parameter('B_2_0', 1))
Initial(B_3(X5=None, ci_g=None), Parameter('B_3_0', 1))
Initial(B_4(X5=None, ci_g=None), Parameter('B_4_0', 1))

Monomer('ci_g_1', ['B_X5'])
Monomer('ci_g_2', ['B_X5'])
Monomer('ci_g_3', ['B_X5'])
Monomer('ci_g_4', ['B_X5'])

Initial(ci_g_1(B_X5=None), Parameter('ci_g_1_0', 1))
Initial(ci_g_2(B_X5=None), Parameter('ci_g_2_0', 1))
Initial(ci_g_3(B_X5=None), Parameter('ci_g_3_0', 1))
Initial(ci_g_4(B_X5=None), Parameter('ci_g_4_0', 1))

Monomer('ci_1')
Monomer('ci_2')
Monomer('ci_3')
Monomer('ci_4')

Initial(ci_1(), Parameter('ci_1_0', 1))
Initial(ci_2(), Parameter('ci_2_0', 1))
Initial(ci_3(), Parameter('ci_3_0', 1))
Initial(ci_4(), Parameter('ci_4_0', 1))

Monomer('CI_1', ['X4', 'X6', 'ptc_g'])
Monomer('CI_2', ['X4', 'X6', 'ptc_g'])
Monomer('CI_3', ['X4', 'X6', 'ptc_g'])
Monomer('CI_4', ['X4', 'X6', 'ptc_g'])

Initial(CI_1(X4=None, X6=None, ptc_g=None), Parameter('CI_1_0', 1))
Initial(CI_2(X4=None, X6=None, ptc_g=None), Parameter('CI_2_0', 1))
Initial(CI_3(X4=None, X6=None, ptc_g=None), Parameter('CI_3_0', 1))
Initial(CI_4(X4=None, X6=None, ptc_g=None), Parameter('CI_4_0', 1))

Monomer('x6_g_1', ['PTC'])
Monomer('x6_g_2', ['PTC'])
Monomer('x6_g_3', ['PTC'])
Monomer('x6_g_4', ['PTC'])

Initial(x6_g_1(PTC=None), Parameter('x6_g_1_0', 1))
Initial(x6_g_2(PTC=None), Parameter('x6_g_2_0', 1))
Initial(x6_g_3(PTC=None), Parameter('x6_g_3_0', 1))
Initial(x6_g_4(PTC=None), Parameter('x6_g_4_0', 1))

Monomer('x6_1')
Monomer('x6_2')
Monomer('x6_3')
Monomer('x6_4')

Initial(x6_1(), Parameter('x6_1_0', 1))
Initial(x6_2(), Parameter('x6_2_0', 1))
Initial(x6_3(), Parameter('x6_3_0', 1))
Initial(x6_4(), Parameter('x6_4_0', 1))

Monomer('X6_1', ['CI'])
Monomer('X6_2', ['CI'])
Monomer('X6_3', ['CI'])
Monomer('X6_4', ['CI'])

Initial(X6_1(CI=None), Parameter('X6_1_0', 1))
Initial(X6_2(CI=None), Parameter('X6_2_0', 1))
Initial(X6_3(CI=None), Parameter('X6_3_0', 1))
Initial(X6_4(CI=None), Parameter('X6_4_0', 1))

Monomer('CN_1', ['x1_g', 'x2_g', 'x4_g', 'x7_g'])
Monomer('CN_2', ['x1_g', 'x2_g', 'x4_g', 'x7_g'])
Monomer('CN_3', ['x1_g', 'x2_g', 'x4_g', 'x7_g'])
Monomer('CN_4', ['x1_g', 'x2_g', 'x4_g', 'x7_g'])

Initial(CN_1(x1_g=None, x2_g=None, x4_g=None, x7_g=None), Parameter('CN_1_0', 1))
Initial(CN_2(x1_g=None, x2_g=None, x4_g=None, x7_g=None), Parameter('CN_2_0', 1))
Initial(CN_3(x1_g=None, x2_g=None, x4_g=None, x7_g=None), Parameter('CN_3_0', 1))
Initial(CN_4(x1_g=None, x2_g=None, x4_g=None, x7_g=None), Parameter('CN_4_0', 1))

Monomer('x7_g_1', ['CN'])
Monomer('x7_g_2', ['CN'])
Monomer('x7_g_3', ['CN'])
Monomer('x7_g_4', ['CN'])

Initial(x7_g_1(CN=None), Parameter('x7_g_1_0', 1))
Initial(x7_g_2(CN=None), Parameter('x7_g_2_0', 1))
Initial(x7_g_3(CN=None), Parameter('x7_g_3_0', 1))
Initial(x7_g_4(CN=None), Parameter('x7_g_4_0', 1))

Monomer('x7_1')
Monomer('x7_2')
Monomer('x7_3')
Monomer('x7_4')

Initial(x7_1(), Parameter('x7_1_0', 1))
Initial(x7_2(), Parameter('x7_2_0', 1))
Initial(x7_3(), Parameter('x7_3_0', 1))
Initial(x7_4(), Parameter('x7_4_0', 1))

Monomer('X7_1', ['EN'])
Monomer('X7_2', ['EN'])
Monomer('X7_3', ['EN'])
Monomer('X7_4', ['EN'])

Initial(X7_1(EN=None), Parameter('X7_1_0', 1))
Initial(X7_2(EN=None), Parameter('X7_2_0', 1))
Initial(X7_3(EN=None), Parameter('X7_3_0', 1))
Initial(X7_4(EN=None), Parameter('X7_4_0', 1))

Monomer('hh_g_1', ['EN_X7'])
Monomer('hh_g_2', ['EN_X7'])
Monomer('hh_g_3', ['EN_X7'])
Monomer('hh_g_4', ['EN_X7'])

Initial(hh_g_1(EN_X7=None), Parameter('hh_g_1_0', 1))
Initial(hh_g_2(EN_X7=None), Parameter('hh_g_2_0', 1))
Initial(hh_g_3(EN_X7=None), Parameter('hh_g_3_0', 1))
Initial(hh_g_4(EN_X7=None), Parameter('hh_g_4_0', 1))

Monomer('hh_1')
Monomer('hh_2')
Monomer('hh_3')
Monomer('hh_4')

Initial(hh_1(), Parameter('hh_1_0', 1))
Initial(hh_2(), Parameter('hh_2_0', 1))
Initial(hh_3(), Parameter('hh_3_0', 1))
Initial(hh_4(), Parameter('hh_4_0', 1))

Monomer('HH_1_1', ['PTC'])
Monomer('HH_1_2', ['PTC'])
Monomer('HH_1_3', ['PTC'])
Monomer('HH_1_4', ['PTC'])
Monomer('HH_1_5', ['PTC'])
Monomer('HH_1_6', ['PTC'])
Monomer('HH_2_1', ['PTC'])
Monomer('HH_2_2', ['PTC'])
Monomer('HH_2_3', ['PTC'])
Monomer('HH_2_4', ['PTC'])
Monomer('HH_2_5', ['PTC'])
Monomer('HH_2_6', ['PTC'])
Monomer('HH_3_1', ['PTC'])
Monomer('HH_3_2', ['PTC'])
Monomer('HH_3_3', ['PTC'])
Monomer('HH_3_4', ['PTC'])
Monomer('HH_3_5', ['PTC'])
Monomer('HH_3_6', ['PTC'])
Monomer('HH_4_1', ['PTC'])
Monomer('HH_4_2', ['PTC'])
Monomer('HH_4_3', ['PTC'])
Monomer('HH_4_4', ['PTC'])
Monomer('HH_4_5', ['PTC'])
Monomer('HH_4_6', ['PTC'])

Initial(HH_1_1(PTC=None), Parameter('HH_1_1_0', 1))
Initial(HH_1_2(PTC=None), Parameter('HH_1_2_0', 1))
Initial(HH_1_3(PTC=None), Parameter('HH_1_3_0', 1))
Initial(HH_1_4(PTC=None), Parameter('HH_1_4_0', 1))
Initial(HH_1_5(PTC=None), Parameter('HH_1_5_0', 1))
Initial(HH_1_6(PTC=None), Parameter('HH_1_6_0', 1))
Initial(HH_2_1(PTC=None), Parameter('HH_2_1_0', 1))
Initial(HH_2_2(PTC=None), Parameter('HH_2_2_0', 1))
Initial(HH_2_3(PTC=None), Parameter('HH_2_3_0', 1))
Initial(HH_2_4(PTC=None), Parameter('HH_2_4_0', 1))
Initial(HH_2_5(PTC=None), Parameter('HH_2_5_0', 1))
Initial(HH_2_6(PTC=None), Parameter('HH_2_6_0', 1))
Initial(HH_3_1(PTC=None), Parameter('HH_3_1_0', 1))
Initial(HH_3_2(PTC=None), Parameter('HH_3_2_0', 1))
Initial(HH_3_3(PTC=None), Parameter('HH_3_3_0', 1))
Initial(HH_3_4(PTC=None), Parameter('HH_3_4_0', 1))
Initial(HH_3_5(PTC=None), Parameter('HH_3_5_0', 1))
Initial(HH_3_6(PTC=None), Parameter('HH_3_6_0', 1))
Initial(HH_4_1(PTC=None), Parameter('HH_4_1_0', 1))
Initial(HH_4_2(PTC=None), Parameter('HH_4_2_0', 1))
Initial(HH_4_3(PTC=None), Parameter('HH_4_3_0', 1))
Initial(HH_4_4(PTC=None), Parameter('HH_4_4_0', 1))
Initial(HH_4_5(PTC=None), Parameter('HH_4_5_0', 1))
Initial(HH_4_6(PTC=None), Parameter('HH_4_6_0', 1))

Monomer('PH_1_1')
Monomer('PH_1_2')
Monomer('PH_1_3')
Monomer('PH_1_4')
Monomer('PH_1_5')
Monomer('PH_1_6')
Monomer('PH_2_1')
Monomer('PH_2_2')
Monomer('PH_2_3')
Monomer('PH_2_4')
Monomer('PH_2_5')
Monomer('PH_2_6')
Monomer('PH_3_1')
Monomer('PH_3_2')
Monomer('PH_3_3')
Monomer('PH_3_4')
Monomer('PH_3_5')
Monomer('PH_3_6')
Monomer('PH_4_1')
Monomer('PH_4_2')
Monomer('PH_4_3')
Monomer('PH_4_4')
Monomer('PH_4_5')
Monomer('PH_4_6')

Initial(PH_1_1(), Parameter('PH_1_1_0', 1))
Initial(PH_1_2(), Parameter('PH_1_2_0', 1))
Initial(PH_1_3(), Parameter('PH_1_3_0', 1))
Initial(PH_1_4(), Parameter('PH_1_4_0', 1))
Initial(PH_1_5(), Parameter('PH_1_5_0', 1))
Initial(PH_1_6(), Parameter('PH_1_6_0', 1))
Initial(PH_2_1(), Parameter('PH_2_1_0', 1))
Initial(PH_2_2(), Parameter('PH_2_2_0', 1))
Initial(PH_2_3(), Parameter('PH_2_3_0', 1))
Initial(PH_2_4(), Parameter('PH_2_4_0', 1))
Initial(PH_2_5(), Parameter('PH_2_5_0', 1))
Initial(PH_2_6(), Parameter('PH_2_6_0', 1))
Initial(PH_3_1(), Parameter('PH_3_1_0', 1))
Initial(PH_3_2(), Parameter('PH_3_2_0', 1))
Initial(PH_3_3(), Parameter('PH_3_3_0', 1))
Initial(PH_3_4(), Parameter('PH_3_4_0', 1))
Initial(PH_3_5(), Parameter('PH_3_5_0', 1))
Initial(PH_3_6(), Parameter('PH_3_6_0', 1))
Initial(PH_4_1(), Parameter('PH_4_1_0', 1))
Initial(PH_4_2(), Parameter('PH_4_2_0', 1))
Initial(PH_4_3(), Parameter('PH_4_3_0', 1))
Initial(PH_4_4(), Parameter('PH_4_4_0', 1))
Initial(PH_4_5(), Parameter('PH_4_5_0', 1))
Initial(PH_4_6(), Parameter('PH_4_6_0', 1))


# Parameters and Rules

# x1 transcription
Parameter('CN_x1_g_kf', 1)
Parameter('CN_x1_g_kr', 1)
Parameter('x1_ts', 1)
Parameter('x1_deg', 1)

# transcription_cooperative binds TFs to the given gene sequentially 

transcription_cooperative(x1_g_1, ['CN'], [CN_1], ['x1_g'], [0], [CN_x1_g_kf, CN_x1_g_kr], x1_1, x1_ts, x1_deg)
transcription_cooperative(x1_g_2, ['CN'], [CN_2], ['x1_g'], [0], [CN_x1_g_kf, CN_x1_g_kr], x1_2, x1_ts, x1_deg)
transcription_cooperative(x1_g_3, ['CN'], [CN_3], ['x1_g'], [0], [CN_x1_g_kf, CN_x1_g_kr], x1_3, x1_ts, x1_deg)
transcription_cooperative(x1_g_4, ['CN'], [CN_4], ['x1_g'], [0], [CN_x1_g_kf, CN_x1_g_kr], x1_4, x1_ts, x1_deg)

# X1 translation
Parameter('X1_tl', 1)
Parameter('X1_deg', 1)

translation(x1_1, X1_1, [X1_tl, X1_deg])
translation(x1_2, X1_2, [X1_tl, X1_deg])
translation(x1_3, X1_3, [X1_tl, X1_deg])
translation(x1_4, X1_4, [X1_tl, X1_deg])

# X1.EWG binding
Parameter('X1_EWG_kf')
Parameter('X1_EWG_kr')

# The intermediate species X1 binds external WG from the faces of surrounding cells.

bind_multiple(X1_1, 'EWG', [EWG_1_4, EWG_2_5, EWG_2_6, EWG_1_1, EWG_4_2, EWG_4_3], ['X1', 'X1', 'X1', 'X1', 'X1', 'X1'],  [X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr])
bind_multiple(X1_2, 'EWG', [EWG_2_4, EWG_3_5, EWG_3_6, EWG_2_1, EWG_1_2, EWG_1_3], ['X1', 'X1', 'X1', 'X1', 'X1', 'X1'],  [X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr])
bind_multiple(X1_3, 'EWG', [EWG_3_4, EWG_4_5, EWG_4_6, EWG_3_1, EWG_2_2, EWG_2_3], ['X1', 'X1', 'X1', 'X1', 'X1', 'X1'],  [X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr])
bind_multiple(X1_4, 'EWG', [EWG_4_4, EWG_1_5, EWG_1_6, EWG_4_1, EWG_3_2, EWG_3_3], ['X1', 'X1', 'X1', 'X1', 'X1', 'X1'],  [X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr, X1_EWG_kf, X1_EWG_kr])

# en transcription
Parameter('EWG_X1_en_g_kf', 1)
Parameter('EWG_X1_en_g_kr', 1)
Parameter('en_ts', 1)
Parameter('en_deg', 1)

# Bound X1.EWG act as a TF for en transcription.

en_g_1_binding_partners = [EWG_1_4(X1=50)%X1_1(EWG=50), EWG_2_5(X1=51)%X1_1(EWG=51), EWG_2_6(X1=52)%X1_1(EWG=52), EWG_1_1(X1=53)%X1_1(EWG=53), EWG_4_2(X1=54)%X1_1(EWG=54), EWG_4_3(X1=55)%X1_1(EWG=55)]
en_g_2_binding_partners = [EWG_2_4(X1=50)%X1_1(EWG=50), EWG_3_5(X1=51)%X1_1(EWG=51), EWG_3_6(X1=52)%X1_1(EWG=52), EWG_2_1(X1=53)%X1_1(EWG=53), EWG_1_2(X1=54)%X1_1(EWG=54), EWG_1_3(X1=55)%X1_1(EWG=55)]
en_g_3_binding_partners = [EWG_3_4(X1=50)%X1_1(EWG=50), EWG_4_5(X1=51)%X1_1(EWG=51), EWG_4_6(X1=52)%X1_1(EWG=52), EWG_3_1(X1=53)%X1_1(EWG=53), EWG_2_2(X1=54)%X1_1(EWG=54), EWG_2_3(X1=55)%X1_1(EWG=55)]
en_g_4_binding_partners = [EWG_4_4(X1=50)%X1_1(EWG=50), EWG_1_5(X1=51)%X1_1(EWG=51), EWG_1_6(X1=52)%X1_1(EWG=52), EWG_4_1(X1=53)%X1_1(EWG=53), EWG_3_2(X1=54)%X1_1(EWG=54), EWG_3_3(X1=55)%X1_1(EWG=55)]

# The single binding site on the en gene accepts multiple binding partners; in this case the banding partners are effectively identical. 

transcription_cooperative(en_g_1, ['EWG_X1'], [en_g_1_binding_partners], ['en_g'], [1], [EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf], en_1, en_ts, en_deg)
transcription_cooperative(en_g_2, ['EWG_X1'], [en_g_2_binding_partners], ['en_g'], [1], [EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf], en_2, en_ts, en_deg)
transcription_cooperative(en_g_3, ['EWG_X1'], [en_g_3_binding_partners], ['en_g'], [1], [EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf], en_3, en_ts, en_deg)
transcription_cooperative(en_g_4, ['EWG_X1'], [en_g_4_binding_partners], ['en_g'], [1], [EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf, EWG_X1_en_g_kf], en_4, en_ts, en_deg)

# EN translation
Parameter('EN_tl', 1)
Parameter('EN_deg', 1)

translation(en_1, EN_1, [EN_tl, EN_deg])
translation(en_2, EN_2, [EN_tl, EN_deg])
translation(en_3, EN_3, [EN_tl, EN_deg])
translation(en_4, EN_4, [EN_tl, EN_deg])

# x2 transcription
Parameter('CN_x2_g_kf', 1)
Parameter('CN_x2_g_kr', 1)
Parameter('x2_ts', 1)
Parameter('x2_deg', 1)

transcription_cooperative(x2_g_1, ['CN'], [CN_1], ['x2_g'], [0], [CN_x2_g_kf, CN_x2_g_kr], x2_1, x2_ts, x2_deg)
transcription_cooperative(x2_g_2, ['CN'], [CN_2], ['x2_g'], [0], [CN_x2_g_kf, CN_x2_g_kr], x2_2, x2_ts, x2_deg)
transcription_cooperative(x2_g_3, ['CN'], [CN_3], ['x2_g'], [0], [CN_x2_g_kf, CN_x2_g_kr], x2_3, x2_ts, x2_deg)
transcription_cooperative(x2_g_4, ['CN'], [CN_4], ['x2_g'], [0], [CN_x2_g_kf, CN_x2_g_kr], x2_4, x2_ts, x2_deg)

# X2 translation
Parameter('X2_tl', 1)
Parameter('X2_deg', 1)

translation(x2_1, X2_1, [X2_tl, X2_deg])
translation(x2_2, X2_2, [X2_tl, X2_deg])
translation(x2_3, X2_3, [X2_tl, X2_deg])
translation(x2_4, X2_4, [X2_tl, X2_deg])

# x3 transcription
Parameter('IWG_x3_g_kf', 1)
Parameter('IWG_x3_g_kr', 1)
Parameter('x3_ts', 1)
Parameter('x3_deg', 1)

transcription_cooperative(x3_g_1, ['IWG'], [IWG_1], ['x3_g'], [1], [IWG_x3_g_kf, IWG_x3_g_kr], x3_1, x3_ts, x3_deg)
transcription_cooperative(x3_g_2, ['IWG'], [IWG_2], ['x3_g'], [1], [IWG_x3_g_kf, IWG_x3_g_kr], x3_2, x3_ts, x3_deg)
transcription_cooperative(x3_g_3, ['IWG'], [IWG_3], ['x3_g'], [1], [IWG_x3_g_kf, IWG_x3_g_kr], x3_3, x3_ts, x3_deg)
transcription_cooperative(x3_g_4, ['IWG'], [IWG_4], ['x3_g'], [1], [IWG_x3_g_kf, IWG_x3_g_kr], x3_4, x3_ts, x3_deg)

# X3 translation
Parameter('X3_tl', 1)
Parameter('X3_deg', 1)

translation(x3_1, X3_1, [X3_tl, X3_deg])
translation(x3_2, X3_2, [X3_tl, X3_deg])
translation(x3_3, X3_3, [X3_tl, X3_deg])
translation(x3_4, X3_4, [X3_tl, X3_deg])

# wg translation
Parameter('X2_wg_g_kf', 1)
Parameter('X2_wg_g_kr', 1)
Parameter('X3_wg_g_kf', 1)
Parameter('X3_wg_g_kr', 1)
Parameter('wg_ts', 1)
Parameter('wg_deg', 1)

transcription_cooperative(wg_g_1, ['X'], [[X2_1, X3_1]], ['wg_g'], [1], [X2_wg_g_kf, X2_wg_g_kr, X3_wg_g_kf, X3_wg_g_kr], wg_1, wg_ts, wg_deg)
transcription_cooperative(wg_g_2, ['X'], [[X2_2, X3_2]], ['wg_g'], [1], [X2_wg_g_kf, X2_wg_g_kr, X3_wg_g_kf, X3_wg_g_kr], wg_2, wg_ts, wg_deg)
transcription_cooperative(wg_g_3, ['X'], [[X2_3, X3_3]], ['wg_g'], [1], [X2_wg_g_kf, X2_wg_g_kr, X3_wg_g_kf, X3_wg_g_kr], wg_3, wg_ts, wg_deg)
transcription_cooperative(wg_g_4, ['X'], [[X2_4, X3_4]], ['wg_g'], [1], [X2_wg_g_kf, X2_wg_g_kr, X3_wg_g_kf, X3_wg_g_kr], wg_4, wg_ts, wg_deg)

# IWG translation
Parameter('IWG_tl', 1)
Parameter('IWG_deg', 1)

translation(wg_1, IWG_1, [IWG_tl, IWG_deg])
translation(wg_2, IWG_2, [IWG_tl, IWG_deg])
translation(wg_3, IWG_3, [IWG_tl, IWG_deg])
translation(wg_4, IWG_4, [IWG_tl, IWG_deg])

# WG exo/endo
Parameter('WG_exo', 1)
Parameter('WG_endo', 1)

# IWG and EWG interconvert through endo- and exocytosis

exo_endo(IWG_1, [EWG_1_1, EWG_1_2, EWG_1_3, EWG_1_4, EWG_1_5, EWG_1_6], [WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo])
exo_endo(IWG_2, [EWG_2_1, EWG_2_2, EWG_2_3, EWG_2_4, EWG_2_5, EWG_2_6], [WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo])
exo_endo(IWG_3, [EWG_3_1, EWG_3_2, EWG_3_3, EWG_3_4, EWG_3_5, EWG_3_6], [WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo])
exo_endo(IWG_4, [EWG_4_1, EWG_4_2, EWG_4_3, EWG_4_4, EWG_4_5, EWG_4_6], [WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo, WG_exo, WG_endo])

# EWG face to face diffusion
Parameter('WG_f2f_diff', 1)

# EWG diffuses between faces that are in contact.

face_to_face_pairs = [[EWG_1_4, EWG_1_1], [EWG_2_5, EWG_1_2], [EWG_2_6, EWG_1_3], [EWG_1_1, EWG_1_4], [EWG_4_2, EWG_1_5], [EWG_4_3, EWG_1_6], [EWG_2_4, EWG_2_1], [EWG_3_5, EWG_2_2], [EWG_3_6, EWG_2_3], [EWG_2_1, EWG_2_4], [EWG_1_2, EWG_2_5], [EWG_1_3, EWG_2_6], \
                    [EWG_3_4, EWG_3_1], [EWG_4_5, EWG_3_2], [EWG_4_6, EWG_3_3], [EWG_3_1, EWG_3_4], [EWG_2_2, EWG_3_5], [EWG_2_3, EWG_3_6], [EWG_4_4, EWG_4_1], [EWG_1_5, EWG_4_2], [EWG_1_6, EWG_4_3], [EWG_4_1, EWG_4_4], [EWG_3_2, EWG_4_5], [EWG_3_3, EWG_4_6]]

face_to_face_diffusion(face_to_face_pairs, WG_f2f_diff)
 
# EWG perimter diffusion
Parameter('EWG_p_diff', 1)

# EWG diffuse around the perimeter of each cell

perimeter_diffusion([EWG_1_1,EWG_1_2,EWG_1_3,EWG_1_4,EWG_1_5,EWG_1_6], EWG_p_diff)
perimeter_diffusion([EWG_2_1,EWG_2_2,EWG_2_3,EWG_2_4,EWG_2_5,EWG_2_6], EWG_p_diff)
perimeter_diffusion([EWG_3_1,EWG_3_2,EWG_3_3,EWG_3_4,EWG_3_5,EWG_3_6], EWG_p_diff)
perimeter_diffusion([EWG_4_1,EWG_4_2,EWG_4_3,EWG_4_4,EWG_4_5,EWG_4_6], EWG_p_diff)
    
# EWG degradation
Parameter('EWG_deg', 1)

EWG_list = [EWG_1_1,EWG_1_2,EWG_1_3,EWG_1_4,EWG_1_5,EWG_1_6,EWG_2_1,EWG_2_2,EWG_2_3,EWG_2_4,EWG_2_5,EWG_2_6,EWG_3_1,EWG_3_2,EWG_3_3,EWG_3_4,EWG_3_5,EWG_3_6,EWG_4_1,EWG_4_2,EWG_4_3,EWG_4_4,EWG_4_5,EWG_4_6]
     
degradation(EWG_list, EWG_deg)
 
# x4 transcription
Parameter('CN_x4_g_kf', 1)
Parameter('CN_x4_g_kr', 1)
Parameter('x4_ts', 1)
Parameter('x4_deg', 1)
 
transcription_cooperative(x4_g_1, ['CN'], [CN_1], ['x4_g'], [0], [CN_x4_g_kf, CN_x4_g_kr], x4_1, x4_ts, x4_deg)
transcription_cooperative(x4_g_2, ['CN'], [CN_2], ['x4_g'], [0], [CN_x4_g_kf, CN_x4_g_kr], x4_2, x4_ts, x4_deg)
transcription_cooperative(x4_g_3, ['CN'], [CN_3], ['x4_g'], [0], [CN_x4_g_kf, CN_x4_g_kr], x4_3, x4_ts, x4_deg)
transcription_cooperative(x4_g_4, ['CN'], [CN_4], ['x4_g'], [0], [CN_x4_g_kf, CN_x4_g_kr], x4_4, x4_ts, x4_deg)
 
# X4 translation
Parameter('X4_tl', 1)
Parameter('X4_deg', 1)

translation(x4_1, X4_1, [X4_tl, X4_deg])
translation(x4_2, X4_2, [X4_tl, X4_deg])
translation(x4_3, X4_3, [X4_tl, X4_deg])
translation(x4_4, X4_4, [X4_tl, X4_deg])

# ptc transcription
Parameter('CI_X4_ptc_g_kf', 1)
Parameter('CI_X4_ptc_g_kr', 1)
Parameter('ptc_ts', 1)
Parameter('ptc_deg', 1)

transcription_cooperative(ptc_g_1, ['CI_X4'], [CI_1(X4=50)%X4_1(CI=50)], ['ptc_g'], [1], [CI_X4_ptc_g_kf, CI_X4_ptc_g_kf], ptc_1, ptc_ts, ptc_deg)
transcription_cooperative(ptc_g_2, ['CI_X4'], [CI_2(X4=50)%X4_2(CI=50)], ['ptc_g'], [1], [CI_X4_ptc_g_kf, CI_X4_ptc_g_kf], ptc_2, ptc_ts, ptc_deg)
transcription_cooperative(ptc_g_3, ['CI_X4'], [CI_3(X4=50)%X4_3(CI=50)], ['ptc_g'], [1], [CI_X4_ptc_g_kf, CI_X4_ptc_g_kf], ptc_3, ptc_ts, ptc_deg)
transcription_cooperative(ptc_g_4, ['CI_X4'], [CI_4(X4=50)%X4_4(CI=50)], ['ptc_g'], [1], [CI_X4_ptc_g_kf, CI_X4_ptc_g_kf], ptc_4, ptc_ts, ptc_deg)

# PTC translation
Parameter('PTC_tl', 1)
Parameter('PTC_deg', 1)

translation(ptc_1, PTC_1_1, [PTC_tl, PTC_deg])
translation(ptc_1, PTC_1_2, [PTC_tl, PTC_deg])
translation(ptc_1, PTC_1_3, [PTC_tl, PTC_deg])
translation(ptc_1, PTC_1_4, [PTC_tl, PTC_deg])
translation(ptc_1, PTC_1_5, [PTC_tl, PTC_deg])
translation(ptc_1, PTC_1_6, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_1, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_2, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_3, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_4, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_5, [PTC_tl, PTC_deg])
translation(ptc_2, PTC_2_6, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_1, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_2, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_3, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_4, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_5, [PTC_tl, PTC_deg])
translation(ptc_3, PTC_3_6, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_1, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_2, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_3, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_4, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_5, [PTC_tl, PTC_deg])
translation(ptc_4, PTC_4_6, [PTC_tl, PTC_deg])

# PTC perimter diffusion
Parameter('PTC_p_diff', 1)

perimeter_diffusion([PTC_1_1,PTC_1_2,PTC_1_3,PTC_1_4,PTC_1_5,PTC_1_6], PTC_p_diff)
perimeter_diffusion([PTC_2_1,PTC_2_2,PTC_2_3,PTC_2_4,PTC_2_5,PTC_2_6], PTC_p_diff)
perimeter_diffusion([PTC_3_1,PTC_3_2,PTC_3_3,PTC_3_4,PTC_3_5,PTC_3_6], PTC_p_diff)
perimeter_diffusion([PTC_4_1,PTC_4_2,PTC_4_3,PTC_4_4,PTC_4_5,PTC_4_6], PTC_p_diff)

# x5 transcription
Parameter('EN_x5_g_kf', 1)
Parameter('EN_x5_g_kr', 1)
Parameter('x5_ts', 1)
Parameter('x5_deg', 1)

transcription_cooperative(x5_g_1, ['EN'], [EN_1], ['x5_g'], [0], [EN_x5_g_kf, EN_x5_g_kr], x5_1, x5_ts, x5_deg)
transcription_cooperative(x5_g_2, ['EN'], [EN_2], ['x5_g'], [0], [EN_x5_g_kf, EN_x5_g_kr], x5_2, x5_ts, x5_deg)
transcription_cooperative(x5_g_3, ['EN'], [EN_3], ['x5_g'], [0], [EN_x5_g_kf, EN_x5_g_kr], x5_3, x5_ts, x5_deg)
transcription_cooperative(x5_g_4, ['EN'], [EN_4], ['x5_g'], [0], [EN_x5_g_kf, EN_x5_g_kr], x5_4, x5_ts, x5_deg)

# X5 translation
Parameter('X5_tl', 1)
Parameter('X5_deg', 1)

translation(x5_1, X5_1, [X5_tl, X5_deg])
translation(x5_2, X5_2, [X5_tl, X5_deg])
translation(x5_3, X5_3, [X5_tl, X5_deg])
translation(x5_4, X5_4, [X5_tl, X5_deg])

# ci transcription
Parameter('B_X5_ci_g_kf', 1)
Parameter('B_X5_ci_g_kr', 1)
Parameter('ci_ts', 1)
Parameter('ci_deg', 1)

transcription_cooperative(ci_g_1, ['B_X5'], [B_1(X5=50)%X5_1(B=50)], ['ci_g'], [1], [B_X5_ci_g_kf, B_X5_ci_g_kf], ci_1, ci_ts, ci_deg)
transcription_cooperative(ci_g_2, ['B_X5'], [B_2(X5=50)%X5_2(B=50)], ['ci_g'], [1], [B_X5_ci_g_kf, B_X5_ci_g_kf], ci_2, ci_ts, ci_deg)
transcription_cooperative(ci_g_3, ['B_X5'], [B_3(X5=50)%X5_3(B=50)], ['ci_g'], [1], [B_X5_ci_g_kf, B_X5_ci_g_kf], ci_3, ci_ts, ci_deg)
transcription_cooperative(ci_g_4, ['B_X5'], [B_4(X5=50)%X5_4(B=50)], ['ci_g'], [1], [B_X5_ci_g_kf, B_X5_ci_g_kf], ci_4, ci_ts, ci_deg)

# CI translation
Parameter('CI_tl', 1)
Parameter('CI_deg', 1)

translation(ci_1, CI_1, [CI_tl, CI_deg])
translation(ci_2, CI_2, [CI_tl, CI_deg])
translation(ci_3, CI_3, [CI_tl, CI_deg])
translation(ci_4, CI_4, [CI_tl, CI_deg])

# x6 transcription
Parameter('PTC_x6_g_kf', 1)
Parameter('PTC_x6_g_kr', 1)
Parameter('x6_ts', 1)
Parameter('x6_deg', 1)

x6_g_1_binding_partners = [PTC_1_1, PTC_1_2, PTC_1_3, PTC_1_4, PTC_1_5, PTC_1_6]
x6_g_2_binding_partners = [PTC_2_1, PTC_2_2, PTC_2_3, PTC_2_4, PTC_2_5, PTC_2_6]
x6_g_3_binding_partners = [PTC_3_1, PTC_3_2, PTC_3_3, PTC_3_4, PTC_3_5, PTC_3_6]
x6_g_4_binding_partners = [PTC_4_1, PTC_4_2, PTC_4_3, PTC_4_4, PTC_4_5, PTC_4_6]
 
transcription_cooperative(x6_g_1, ['PTC'], [x6_g_1_binding_partners], ['x6_g'], [1], [PTC_x6_g_kf, PTC_x6_g_kf], x6_1, x6_ts, x6_deg)
transcription_cooperative(x6_g_2, ['PTC'], [x6_g_2_binding_partners], ['x6_g'], [1], [PTC_x6_g_kf, PTC_x6_g_kf], x6_2, x6_ts, x6_deg)
transcription_cooperative(x6_g_3, ['PTC'], [x6_g_3_binding_partners], ['x6_g'], [1], [PTC_x6_g_kf, PTC_x6_g_kf], x6_3, x6_ts, x6_deg)
transcription_cooperative(x6_g_4, ['PTC'], [x6_g_4_binding_partners], ['x6_g'], [1], [PTC_x6_g_kf, PTC_x6_g_kf], x6_4, x6_ts, x6_deg)

# X6 translation
Parameter('X6_tl', 1)
Parameter('X6_deg', 1)

translation(x6_1, X6_1, [X6_tl, X6_deg])
translation(x6_2, X6_2, [X6_tl, X6_deg])
translation(x6_3, X6_3, [X6_tl, X6_deg])
translation(x6_4, X6_4, [X6_tl, X6_deg])

# CI-X6 binding
Parameter('CI_X6_kf', 1)
Parameter('CI_X6_kr', 1)

# X6 binds CI and induces cleavage to form CN

Rule('Bind_CI_X6_1', CI_1(X4=None, ptc_g=None, X6=None) + X6_1(CI=None) <> CI_1(X4=None, ptc_g=None, X6=1)%X6_1(CI=1), CI_X6_kf, CI_X6_kr)
Rule('Bind_CI_X6_2', CI_2(X4=None, ptc_g=None, X6=None) + X6_2(CI=None) <> CI_2(X4=None, ptc_g=None, X6=1)%X6_2(CI=1), CI_X6_kf, CI_X6_kr)
Rule('Bind_CI_X6_3', CI_3(X4=None, ptc_g=None, X6=None) + X6_3(CI=None) <> CI_3(X4=None, ptc_g=None, X6=1)%X6_3(CI=1), CI_X6_kf, CI_X6_kr)
Rule('Bind_CI_X6_4', CI_4(X4=None, ptc_g=None, X6=None) + X6_4(CI=None) <> CI_4(X4=None, ptc_g=None, X6=1)%X6_4(CI=1), CI_X6_kf, CI_X6_kr)
 
# CI cleavage
Parameter('Cleave_CI', 1) 

Rule('Cleave_CI_1', CI_1(X4=None, ptc_g=None, X6=1)%X6_1(CI=1) >> CN_1(x1_g=None, x2_g=None, x4_g=None, x7_g=None) + X6_1(CI=None), Cleave_CI)
Rule('Cleave_CI_2', CI_2(X4=None, ptc_g=None, X6=1)%X6_2(CI=1) >> CN_2(x1_g=None, x2_g=None, x4_g=None, x7_g=None) + X6_2(CI=None), Cleave_CI)
Rule('Cleave_CI_3', CI_3(X4=None, ptc_g=None, X6=1)%X6_3(CI=1) >> CN_3(x1_g=None, x2_g=None, x4_g=None, x7_g=None) + X6_3(CI=None), Cleave_CI)
Rule('Cleave_CI_4', CI_4(X4=None, ptc_g=None, X6=1)%X6_4(CI=1) >> CN_4(x1_g=None, x2_g=None, x4_g=None, x7_g=None) + X6_4(CI=None), Cleave_CI)

# CN degradation
Parameter('CN_deg', 1)

degradation([CN_1, CN_2, CN_3, CN_4], CN_deg)
 
# x7 transcription
Parameter('CN_x7_g_kf', 1)
Parameter('CN_x7_g_kr', 1)
Parameter('x7_ts', 1)
Parameter('x7_deg', 1) 

transcription_cooperative(x7_g_1, ['CN'], [CN_1], ['x7_g'], [0], [CN_x7_g_kf, CN_x7_g_kr], x7_1, x7_ts, x7_deg)
transcription_cooperative(x7_g_2, ['CN'], [CN_2], ['x7_g'], [0], [CN_x7_g_kf, CN_x7_g_kr], x7_2, x7_ts, x7_deg)
transcription_cooperative(x7_g_3, ['CN'], [CN_3], ['x7_g'], [0], [CN_x7_g_kf, CN_x7_g_kr], x7_3, x7_ts, x7_deg)
transcription_cooperative(x7_g_4, ['CN'], [CN_4], ['x7_g'], [0], [CN_x7_g_kf, CN_x7_g_kr], x7_4, x7_ts, x7_deg)

# X7 translation
Parameter('X7_tl', 1)
Parameter('X7_deg', 1)

translation(x7_1, X7_1, [X7_tl, X7_deg])
translation(x7_2, X7_2, [X7_tl, X7_deg])
translation(x7_3, X7_3, [X7_tl, X7_deg])
translation(x7_4, X7_4, [X7_tl, X7_deg])

# hh transcription
Parameter('EN_X7_hh_g_kf', 1)
Parameter('EN_X7_hh_g_kr', 1)
Parameter('hh_ts', 1)
Parameter('hh_deg', 1)

transcription_cooperative(hh_g_1, ['EN_X7'], [EN_1(X7=50)%X7_1(EN=50)], ['hh_g'], [1], [EN_X7_hh_g_kf, EN_X7_hh_g_kf], hh_1, hh_ts, hh_deg)
transcription_cooperative(hh_g_2, ['EN_X7'], [EN_2(X7=50)%X7_2(EN=50)], ['hh_g'], [1], [EN_X7_hh_g_kf, EN_X7_hh_g_kf], hh_2, hh_ts, hh_deg)
transcription_cooperative(hh_g_3, ['EN_X7'], [EN_3(X7=50)%X7_3(EN=50)], ['hh_g'], [1], [EN_X7_hh_g_kf, EN_X7_hh_g_kf], hh_3, hh_ts, hh_deg)
transcription_cooperative(hh_g_4, ['EN_X7'], [EN_4(X7=50)%X7_4(EN=50)], ['hh_g'], [1], [EN_X7_hh_g_kf, EN_X7_hh_g_kf], hh_4, hh_ts, hh_deg)

# HH translation
Parameter('HH_tl', 1)
Parameter('HH_deg', 1)

translation(hh_1, HH_1_1, [HH_tl, HH_deg])
translation(hh_1, HH_1_2, [HH_tl, HH_deg])
translation(hh_1, HH_1_3, [HH_tl, HH_deg])
translation(hh_1, HH_1_4, [HH_tl, HH_deg])
translation(hh_1, HH_1_5, [HH_tl, HH_deg])
translation(hh_1, HH_1_6, [HH_tl, HH_deg])
translation(hh_2, HH_2_1, [HH_tl, HH_deg])
translation(hh_2, HH_2_2, [HH_tl, HH_deg])
translation(hh_2, HH_2_3, [HH_tl, HH_deg])
translation(hh_2, HH_2_4, [HH_tl, HH_deg])
translation(hh_2, HH_2_5, [HH_tl, HH_deg])
translation(hh_2, HH_2_6, [HH_tl, HH_deg])
translation(hh_3, HH_3_1, [HH_tl, HH_deg])
translation(hh_3, HH_3_2, [HH_tl, HH_deg])
translation(hh_3, HH_3_3, [HH_tl, HH_deg])
translation(hh_3, HH_3_4, [HH_tl, HH_deg])
translation(hh_3, HH_3_5, [HH_tl, HH_deg])
translation(hh_3, HH_3_6, [HH_tl, HH_deg])
translation(hh_4, HH_4_1, [HH_tl, HH_deg])
translation(hh_4, HH_4_2, [HH_tl, HH_deg])
translation(hh_4, HH_4_3, [HH_tl, HH_deg])
translation(hh_4, HH_4_4, [HH_tl, HH_deg])
translation(hh_4, HH_4_5, [HH_tl, HH_deg])
translation(hh_4, HH_4_6, [HH_tl, HH_deg])

# HH perimter diffusion
Parameter('HH_p_diff', 1)

perimeter_diffusion([HH_1_1,HH_1_2,HH_1_3,HH_1_4,HH_1_5,HH_1_6], HH_p_diff)
perimeter_diffusion([HH_2_1,HH_2_2,HH_2_3,HH_2_4,HH_2_5,HH_2_6], HH_p_diff)
perimeter_diffusion([HH_3_1,HH_3_2,HH_3_3,HH_3_4,HH_3_5,HH_3_6], HH_p_diff)
perimeter_diffusion([HH_4_1,HH_4_2,HH_4_3,HH_4_4,HH_4_5,HH_4_6], HH_p_diff)

# HH-PTC binding to form PH
Parameter('PH_kf', 1)
Parameter('PH_kr', 1)
Parameter('PH_deg', 1)

# PTC and HH bind to form the dimer PH

dimerization(PTC_1_1, 'HH', HH_1_4, 'PTC', PH_1_1, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_1_2, 'HH', HH_2_5, 'PTC', PH_1_2, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_1_3, 'HH', HH_2_6, 'PTC', PH_1_3, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_1_4, 'HH', HH_1_1, 'PTC', PH_1_4, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_1_5, 'HH', HH_4_2, 'PTC', PH_1_5, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_1_6, 'HH', HH_4_3, 'PTC', PH_1_6, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_1, 'HH', HH_2_4, 'PTC', PH_2_1, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_2, 'HH', HH_3_5, 'PTC', PH_2_2, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_3, 'HH', HH_3_6, 'PTC', PH_2_3, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_4, 'HH', HH_2_1, 'PTC', PH_2_4, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_5, 'HH', HH_1_2, 'PTC', PH_2_5, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_2_6, 'HH', HH_1_3, 'PTC', PH_2_6, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_1, 'HH', HH_3_4, 'PTC', PH_3_1, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_2, 'HH', HH_4_5, 'PTC', PH_3_2, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_3, 'HH', HH_4_6, 'PTC', PH_3_3, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_4, 'HH', HH_3_1, 'PTC', PH_3_4, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_5, 'HH', HH_2_2, 'PTC', PH_3_5, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_3_6, 'HH', HH_2_3, 'PTC', PH_3_6, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_1, 'HH', HH_4_4, 'PTC', PH_4_1, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_2, 'HH', HH_1_5, 'PTC', PH_4_2, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_3, 'HH', HH_1_6, 'PTC', PH_4_3, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_4, 'HH', HH_4_1, 'PTC', PH_4_4, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_5, 'HH', HH_3_2, 'PTC', PH_4_5, [PH_kf, PH_kr, PH_deg])
dimerization(PTC_4_6, 'HH', HH_3_3, 'PTC', PH_4_6, [PH_kf, PH_kr, PH_deg])



