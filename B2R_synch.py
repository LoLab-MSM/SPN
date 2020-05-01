import numpy as np
from pysb.integrate import Solver
from copy import deepcopy
from pysb import ComponentSet
from test import model as model0
import matplotlib.pyplot as plt
import re
from pysb.bng import run_ssa

monomer_names = [m.name for m in model0.monomers]
models = [deepcopy(model0) for m in monomer_names]
for i,m in enumerate(monomer_names):
    models[i].rules = ComponentSet([r for r in model0.rules if re.match(m, r.name)])

n_steps = 10
timecourse = []
param_values = np.array([p.value for p in model0.parameters])
for step in range(1, n_steps+1):
    print step
    if step == 1:
        timecourse.append([])
    timecourse.append([])
    for i,m in enumerate(models):
        x = run_ssa(m, t_end=100, n_steps=1, param_values=param_values, verbose=False, max_sim_steps=1)
        if step == 1:
            timecourse[0].append(x[m.observables[i].name][0])
        timecourse[step].append(x[m.observables[i].name][-1])
    param_values[range(1, 2*len(monomer_names), 2)] = timecourse[step]
    param_values[range(0, 2*len(monomer_names), 2)] = [1-t for t in timecourse[step]]
timecourse = np.array(timecourse).T

print timecourse

for i,tc in enumerate(timecourse):
    plt.plot(range(n_steps+1), tc, lw=2, label=monomer_names[i])
plt.legend(loc=0)

plt.show()

quit()

# ~~~ MAK ~~~
models = []
for each in model0.rules:
    models.append(deepcopy(model0))
    models[-1].rules = ComponentSet([each])

params = []
for each in model0.parameters:
    print each
    params.append(each.value)
tail = params[-len(model0.monomers):]
print tail

t = np.linspace(0, 10, 11)
steps = 10
trace = [deepcopy(params[:-len(model0.monomers)])]
print trace

for _ in range(steps):

    new_params = deepcopy(params[:-len(model0.monomers)])
    for each in models:
        s = Solver(each, t)
        s.run(param_values=params)
        m = list(s.y[-1])
        for i,each in enumerate(m):
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
    
print
print np.array(trace).shape
A_T = [t[1] for t in trace]
B_T = [t[3] for t in trace]
C_T = [t[5] for t in trace]

plt.plot(range(steps+1), A_T, lw=2, label='A_true')
plt.plot(range(steps+1), B_T, lw=2, label='B_true')
plt.plot(range(steps+1), C_T, lw=2, label='C_true')
plt.legend(loc=0)

plt.show()

