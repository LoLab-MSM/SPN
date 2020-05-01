
from pysb import Model, Monomer, Parameter, Rule, Observable, Initial

Model()

Monomer('A', ['state'], {'state': ['off', 'on']})
Monomer('B', ['state'], {'state': ['off', 'on']})
Monomer('C', ['state'], {'state': ['off', 'on']})

Parameter('A_off_0', 1.0)
Parameter('A_on_0', 0.0)
Parameter('B_off_0', 1.0)
Parameter('B_on_0', 0.0)
Parameter('C_off_0', 0.0)
Parameter('C_on_0', 1.0)
Parameter('rate_A', 1.0)
Parameter('rate_B', 1.0)
Parameter('rate_C', 1.0)

# Observable('A_off_obs', A(state='off'))
Observable('A_on_obs', A(state='on'))
# Observable('B_off_obs', B(state='off'))
Observable('B_on_obs', B(state='on'))
# Observable('C_off_obs', C(state='off'))
Observable('C_on_obs', C(state='on'))

# Initials
Initial(A(state='off'), A_off_0)
Initial(A(state='on'), A_on_0)
Initial(B(state='off'), B_off_0)
Initial(B(state='on'), B_on_0)
Initial(C(state='off'), C_off_0)
Initial(C(state='on'), C_on_0)

Rule('A_rule0', A(state='off') + B(state='on') >> A(state='on') + B(state='on'), rate_A)
Rule('A_rule1', A(state='off') + B(state='off') + C(state='on') >> A(state='on') + B(state='off') + C(state='on'), rate_A)
Rule('A_rule2', A(state='on') + B(state='off') + C(state='off') >> A(state='off') + B(state='off') + C(state='off'), rate_A)
Rule('B_rule0', B(state='on') + A(state='on') >> B(state='off') + A(state='on'), rate_B)
Rule('B_rule1', B(state='off') + A(state='off') + C(state='on') >> B(state='on') + A(state='off') + C(state='on'), rate_B)
Rule('B_rule2', B(state='on') + A(state='off') + C(state='off') >> B(state='off') + A(state='off') + C(state='off'), rate_B)
Rule('C_rule0', C(state='on') + A(state='on') >> C(state='off') + A(state='on'), rate_C)
Rule('C_rule1', C(state='off') + A(state='on') >> C(state='on') + A(state='on'), rate_C)
Rule('C_rule2', C(state='on') + A(state='off') >> C(state='off') + A(state='off'), rate_C)
