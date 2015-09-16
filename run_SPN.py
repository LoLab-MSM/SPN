
from pysb.integrate import Solver
from sympy import sympify
from Model_SPN import model
from pysb.bng import generate_equations
from pysb.integrate import odesolve
import numpy as np

generate_equations(model, verbose = True)
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
exp1 = model.expressions["EWG_nT_1"]
print exp1
print exp1.expand_expr()
print

exp2 = model.expressions["Hill_EWG_nT_1"]
print exp2
print exp2.expand_expr()


#---------------------------------------------------- t=np.linspace(0, 100, 101)
#------------------------------------------ x = odesolve(model, t, verbose=True)