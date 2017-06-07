
import numpy as np
from pysb.integrate import Solver
from copy import deepcopy
from pysb import ComponentSet

from test import model as model0

models = []
for each in model0.rules:
    models.append(deepcopy(model0))
    models[-1].rules = ComponentSet([each])

params = []
for each in model0.parameters:
    params.append(each.value)
tail = params[-len(model0.monomers):]
t = np.linspace(0, 10, 11)
steps = 10
trace = [deepcopy(params[:-len(model0.monomers)])]

for _ in range(steps):

    new_params = deepcopy(params[:-len(model0.monomers)])
    for each in models:
        s = Solver(each, t)
        s.run(params)
        m = list(s.y[-1])
        for i, each in enumerate(m):
            if each < 0.5:
                m[i] = 0.0
            else:
                m[i] = 1.0
        for i,each in enumerate(m):
            if each != params[i]:
                new_params[i] = each

    params = new_params + tail
    trace.append(deepcopy(new_params))

print ' A_F, A_T, B_F, B_T, C_F, C_T'
for i,each in enumerate(trace):
    print i, each
    
