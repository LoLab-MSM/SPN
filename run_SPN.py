
from pysb.integrate import Solver
from pysb.bng import run_ssa
from sympy import sympify
from Model_SPN import model
from pysb.bng import generate_equations
from pysb.integrate import odesolve
import numpy as np
#from sys import argv
import matplotlib.pyplot as plt
import matplotlib.patches as patches

t=np.linspace(0, 500, 501)
#===============================================================================
# with open("Ingeneue_time", "r") as ins:
#     t = []
#     for line in ins:
#         t.append(float(line.rstrip('\n')))
#===============================================================================
#t = [float(line.rstrip('\n')) for line in open('Ingeneue_time')]

#generate_equations(model, verbose = True)
x = odesolve(model, t, verbose=True)
#x = run_ssa(model, t[-1], len(t)-1, verbose=True)
#print model.parameters
 #==============================================================================
 # for ode in model.odes:
 #     print ode
 #==============================================================================
#print x["hh_1_obs"]
EWG_1 = x["EWG_1_1_obs"][100] + x["EWG_1_2_obs"][100] + x["EWG_1_3_obs"][100] + x["EWG_1_4_obs"][100] + x["EWG_1_5_obs"][100] + x["EWG_1_6_obs"][100]
EWG_2 = x["EWG_2_1_obs"][100] + x["EWG_2_2_obs"][100] + x["EWG_2_3_obs"][100] + x["EWG_2_4_obs"][100] + x["EWG_2_5_obs"][100] + x["EWG_2_6_obs"][100]
EWG_3 = x["EWG_3_1_obs"][100] + x["EWG_3_2_obs"][100] + x["EWG_3_3_obs"][100] + x["EWG_3_4_obs"][100] + x["EWG_3_5_obs"][100] + x["EWG_3_6_obs"][100]
EWG_4 = x["EWG_4_1_obs"][100] + x["EWG_4_2_obs"][100] + x["EWG_4_3_obs"][100] + x["EWG_4_4_obs"][100] + x["EWG_4_5_obs"][100] + x["EWG_4_6_obs"][100]
print EWG_1
print EWG_2
print EWG_3
print EWG_4
PTC_1 = x["PTC_1_1_obs"][100] + x["PTC_1_2_obs"][100] + x["PTC_1_3_obs"][100] + x["PTC_1_4_obs"][100] + x["PTC_1_5_obs"][100] + x["PTC_1_6_obs"][100]
PTC_2 = x["PTC_2_1_obs"][100] + x["PTC_2_2_obs"][100] + x["PTC_2_3_obs"][100] + x["PTC_2_4_obs"][100] + x["PTC_2_5_obs"][100] + x["PTC_2_6_obs"][100]
PTC_3 = x["PTC_3_1_obs"][100] + x["PTC_3_2_obs"][100] + x["PTC_3_3_obs"][100] + x["PTC_3_4_obs"][100] + x["PTC_3_5_obs"][100] + x["PTC_3_6_obs"][100]
PTC_4 = x["PTC_4_1_obs"][100] + x["PTC_4_2_obs"][100] + x["PTC_4_3_obs"][100] + x["PTC_4_4_obs"][100] + x["PTC_4_5_obs"][100] + x["PTC_4_6_obs"][100]
print PTC_1
print PTC_2
print PTC_3
print PTC_4
HH_1 = x["HH_1_1_obs"][100] + x["HH_1_2_obs"][100] + x["HH_1_3_obs"][100] + x["HH_1_4_obs"][100] + x["HH_1_5_obs"][100] + x["HH_1_6_obs"][100]
HH_2 = x["HH_2_1_obs"][100] + x["HH_2_2_obs"][100] + x["HH_2_3_obs"][100] + x["HH_2_4_obs"][100] + x["HH_2_5_obs"][100] + x["HH_2_6_obs"][100]
HH_3 = x["HH_3_1_obs"][100] + x["HH_3_2_obs"][100] + x["HH_3_3_obs"][100] + x["HH_3_4_obs"][100] + x["HH_3_5_obs"][100] + x["HH_3_6_obs"][100]
HH_4 = x["HH_4_1_obs"][100] + x["HH_4_2_obs"][100] + x["HH_4_3_obs"][100] + x["HH_4_4_obs"][100] + x["HH_4_5_obs"][100] + x["HH_4_6_obs"][100]
#-------------------------------------------------------------------- print HH_1
#-------------------------------------------------------------------- print HH_2
#-------------------------------------------------------------------- print HH_3
print HH_4
PH_1 = x["PH_1_1_obs"][100] + x["PH_1_2_obs"][100] + x["PH_1_3_obs"][100] + x["PH_1_4_obs"][100] + x["PH_1_5_obs"][100] + x["PH_1_6_obs"][100]
PH_2 = x["PH_2_1_obs"][100] + x["PH_2_2_obs"][100] + x["PH_2_3_obs"][100] + x["PH_2_4_obs"][100] + x["PH_2_5_obs"][100] + x["PH_2_6_obs"][100]
PH_3 = x["PH_3_1_obs"][100] + x["PH_3_2_obs"][100] + x["PH_3_3_obs"][100] + x["PH_3_4_obs"][100] + x["PH_3_5_obs"][100] + x["PH_3_6_obs"][100]
PH_4 = x["PH_4_1_obs"][100] + x["PH_4_2_obs"][100] + x["PH_4_3_obs"][100] + x["PH_4_4_obs"][100] + x["PH_4_5_obs"][100] + x["PH_4_6_obs"][100]

print PH_1
print PH_2
print PH_3
print PH_4

#===============================================================================
# maxPH = max(PH_1, PH_2, PH_3, PH_4)
# 
# PH_1 = PH_1/maxPH
# PH_2 = PH_2/maxPH
# PH_3 = PH_3/maxPH
# PH_4 = PH_4/maxPH
#===============================================================================



print x["en_1_obs"][100]
print x["en_2_obs"][100]
print x["en_3_obs"][100]
print x["en_4_obs"][100]
print x["EN_1_obs"][100]
print x["EN_2_obs"][100]
print x["EN_3_obs"][100]
print x["EN_4_obs"][100]
print x["wg_1_obs"][100]
print x["wg_2_obs"][100]
print x["wg_3_obs"][100]
print x["wg_4_obs"][100]
print x["IWG_1_obs"][100]
print x["IWG_2_obs"][100]
print x["IWG_3_obs"][100]
print x["IWG_4_obs"][100]
print x["ptc_1_obs"][100]
print x["ptc_2_obs"][100]
print x["ptc_3_obs"][100]
print x["ptc_4_obs"][100]
print x["ci_1_obs"][100]
print x["ci_2_obs"][100]
print x["ci_3_obs"][100]
print x["ci_4_obs"][100]
print x["CI_1_obs"][100]
print x["CI_2_obs"][100]
print x["CI_3_obs"][100]
print x["CI_4_obs"][100]
print x["CN_1_obs"][100]
print x["CN_2_obs"][100]
print x["CN_3_obs"][100]
print x["CN_4_obs"][100]
print x["hh_1_obs"][100]
print x["hh_2_obs"][100]
print x["hh_3_obs"][100]
print x["hh_4_obs"][100]

fig1 = plt.figure(figsize=(5.1, 15))
ax1 = fig1.add_subplot(111, aspect='equal')

ax1.axes.get_xaxis().set_visible(False)
ax1.axes.get_yaxis().set_visible(False)
ax1.add_patch(
    patches.RegularPolygon(
        (2, 5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["en_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, 5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["en_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, 5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["en_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, 5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["en_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, 2.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["EN_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, 2.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["EN_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, 2.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["EN_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, 2.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,           # radius
        np.pi/6,
        alpha = min(x["EN_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, 0),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["wg_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, 0 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["wg_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, 0),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["wg_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, 0 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["wg_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -2.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["IWG_4_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -2.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["IWG_4_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -2.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["IWG_4_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -2.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["IWG_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(EWG_1, 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(EWG_2, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(EWG_3, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(EWG_4, 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -7.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ptc_4_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -7.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ptc_4_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -7.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ptc_4_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -7.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ptc_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -10),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PTC_1, 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -10 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PTC_2, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -10),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PTC_3, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -10 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PTC_4, 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -12.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ci_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -12.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ci_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -12.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ci_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -12.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["ci_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -15),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CI_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -15 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CI_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -15),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CI_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -15 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CI_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -17.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CN_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -17.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CN_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -17.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CN_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -17.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["CN_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -20),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["hh_1_obs"][100], 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -20 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["hh_2_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -20),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["hh_3_obs"][100], 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -20 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(x["hh_4_obs"][100], 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -22.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(HH_1, 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -22.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(HH_2, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -22.5),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(HH_3, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -22.5 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(HH_4, 1)
    )            
)

ax1.add_patch(
    patches.RegularPolygon(
        (2, -25),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PH_1, 1)
    )
)
ax1.add_patch(
    patches.RegularPolygon(
        (3.5, -25 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PH_2, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (5, -25),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PH_3, 1)
    )            
)
ax1.add_patch(
    patches.RegularPolygon(
        (6.5, -25 - np.cos(np.pi/6)),     # (x,y)
        6,              # number of vertices
        1,            # radius
        np.pi/6,
        alpha = min(PH_4, 1)
    )            
)

plt.annotate('en', xy=(-1.25, 5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('EN', xy=(-1.25, 2.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('wg', xy=(-1.25, 0 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('IWG', xy=(-1.25, -2.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('EWG', xy=(-1.25, -5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('ptc', xy=(-1.25, -7.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('PTC', xy=(-1.25, -10 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('ci', xy=(-1.25, -12.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('CI', xy=(-1.25, -15 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('CN', xy=(-1.25, -17.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('hh', xy=(-1.25, -20 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('HH', xy=(-1.25, -22.5 - np.cos(np.pi/6)/2), fontsize = 25)
plt.annotate('PH', xy=(-1.25, -25 - np.cos(np.pi/6)/2), fontsize = 25)


plt.axis([-2, 8, -27.5, 6.5])
plt.tight_layout()
plt.show()



#===============================================================================
# for i,sp in enumerate(model.species):
#     print i,sp
# print
# 
# for i,rxn in enumerate(model.reactions):
#     print i,rxn
# print
# for i,ode in enumerate(model.odes):
#     print i,ode
# print
# 
# for i,exp in enumerate(model.expressions):
#     print i,exp
#     print exp.expand_expr()
# print
#===============================================================================


#------------------------------------- exp2 = model.expressions["Hill_EWG_nT_1"]
#-------------------------------------------------------------------- print exp2
#------------------------------------------------------ print exp2.expand_expr()


