from pysb import *
from pysb.integrate import odesolve
import numpy as np
import matplotlib.pyplot as plt
import wg_boolean

# wg' := (CIA and SLP and not CIR) or (wg and (CIA or SLP) and not CIR)
# m' := (T1 and T2 and not R) or (m and (T1 or T2) and not R)
# m' := T1*T2*(1-R) + m*(T1+T2)*(1-R) = (T1*T2 + m*(T1+T2))*(1-R) = ((T1 and T2 or m and (T1 or T2)) and not R 

# print ( ( True and True or True and ( True or True ) ) and not True )
# quit()

from pysb.importers.boolean import model_from_boolean
model = model_from_boolean('wg.txt', mode='GSP')

print model.parameters
print
print model.rules

quit()

Model()

###############
T1_0 = 100
T2_0 = 100
R_0 = 0
P_0 = 0
###############

wg_boolean.monomers()
wg_boolean.initials(T1_0, T2_0, R_0, P_0)
wg_boolean.rules(version=3)
wg_boolean.observables()

# Monomer('G', ['t1', 't2', 'p', 'r', 'state'], {'state' : ['off', 'on']})
# # Monomer('G', ['t1', 't2', 'p', 'r'])
# Monomer('T1', ['g'])
# Monomer('T2', ['g'])
# Monomer('R', ['g'])
# Monomer('M')
# Monomer('P', ['g'])
# 
# Parameter('G_init', 1)
# Parameter('T1_init', T1_0)
# Parameter('T2_init', T2_0)
# Parameter('R_init', R_0)
# Parameter('P_init', P_0)
# 
# Initial(G(t1=None, t2=None, p=None, r=None, state='off'), G_init)
# # Initial(G(t1=None, t2=None, p=None, r=None), G_init)
# Initial(T1(g=None), T1_init)
# Initial(T2(g=None), T2_init)
# Initial(R(g=None), R_init)
# Initial(P(g=None), P_init)
# 
# Parameter('kgt1f', 10)
# Parameter('kgt1r', 1)
# Parameter('kgt2f', 10)
# Parameter('kgt2r', 1)
# Parameter('kgpf', 10)
# Parameter('kgpr', 1)
# Parameter('kgrf', 100)
# Parameter('kgrr', 1)
# Parameter('kg_act', 10)
# Parameter('kg_deact', 1)
# Parameter('k_transcribe', 1)
# Parameter('k_translate', 10)
# Parameter('kmrna_degrade', 1)
# Parameter('kp_degrade', 1)
# 
# Rule('T1_binds_G', G(t1=None) + T1(g=None) <> G(t1=1) % T1(g=1), kgt1f, kgt1r)
# Rule('T2_binds_G', G(t2=None) + T2(g=None) <> G(t2=1) % T2(g=1), kgt2f, kgt2r)
# Rule('P_binds_G', G(p=None) + P(g=None) <> G(p=1) % P(g=1), kgpf, kgpr)
# Rule('R_binds_G', G(r=None) + R(g=None) <> G(r=1) % R(g=1), kgrf, kgrr)
# #####
# Rule('GT1T2_activates', G(t1=ANY, t2=ANY, r=None, state='off') >> G(t1=ANY, t2=ANY, r=None, state='on'), kg_act)
# Rule('GT1P_activates', G(t1=ANY, p=ANY, r=None, state='off') >> G(t1=ANY, p=ANY, r=None, state='on'), kg_act)
# Rule('GT2P_activates', G(t2=ANY, p=ANY, r=None, state='off') >> G(t2=ANY, p=ANY, r=None, state='on'), kg_act)
# # Rule('G_deactivates', G(state='on') >> G(state='off'), kg_deact)
# Rule('G_deactivates', G(t1=None, t2=None, p=None, state='on') >> G(t1=None, t2=None, p=None, state='off'), kg_deact)
# #####
# 
# Rule('G_transcribes', G(state='on') >> G(state='on') + M(), k_transcribe)
# # Rule('GT1T2_transcribes', G(t1=ANY, t2=ANY) >> G(t1=ANY, t2=ANY) + M(), k_transcribe)
# # Rule('GT1P_transcribes', G(t1=ANY, p=ANY) >> G(t1=ANY, p=ANY) + M(), k_transcribe)
# # Rule('GT2P_transcribes', G(t2=ANY, p=ANY) >> G(t2=ANY, p=ANY) + M(), k_transcribe)
# Rule('mRNA_translates', M() >> M() + P(g=None), k_translate)
# Rule('mRNA_degrades', M() >> None, kmrna_degrade)
# Rule('P_degrades', P(g=None) >> None, kp_degrade)

################
Parameter('kt1_degrade', 0)
Rule('T1_degrades', T1(g=None) >> None, kt1_degrade)
Observable('T1_tot', T1())
################
Parameter('kt2_degrade', 0)
Rule('T2_degrades', T2(g=None) >> None, kt2_degrade)
Observable('T2_tot', T2())
################

# Observable('mRNA', M())

tspan = np.linspace(0, 1000, 1001)
x = odesolve(model, tspan, verbose=True) 

plt.figure('M')
plt.plot(tspan, x['mRNA'], lw=2)

plt.figure('T1')
plt.plot(tspan, x['T1_tot'], lw=2)

plt.figure('T2')
plt.plot(tspan, x['T2_tot'], lw=2)

##############
# for i,sp in enumerate(model.species):
#     print sp, x['__s%d' % i][-1]
##############

init_pops = []
for i in range(len(model.species)):
    init_pops.append(x[-1][i])

x = odesolve(model, tspan, y0 = init_pops, param_values = {'kt1_degrade' : 100}, verbose=True)

plt.figure('M')
plt.plot(tspan+tspan[-1], x['mRNA'], lw=2)

plt.figure('T1')
plt.plot(tspan+tspan[-1], x['T1_tot'], lw=2)

plt.figure('T2')
plt.plot(tspan+tspan[-1], x['T2_tot'], lw=2)

#####
init_pops = []
for i in range(len(model.species)):
    init_pops.append(x[-1][i])

x = odesolve(model, tspan, y0 = init_pops, param_values = {'kt1_degrade' : 100, 'kt2_degrade' : 100}, verbose=True)

plt.figure('M')
plt.plot(tspan+2*tspan[-1], x['mRNA'], lw=2)

plt.figure('T1')
plt.plot(tspan+2*tspan[-1], x['T1_tot'], lw=2)

plt.figure('T2')
plt.plot(tspan+2*tspan[-1], x['T2_tot'], lw=2)

plt.show()
   




