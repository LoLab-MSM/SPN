import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('SPN_Ingenue.txt', names = ['pset', 'time', 'cell', 'species', 'conc'])

sp_names = ['en', 'EN', 'wg', 'IWG', 'EWG', 'ptc', 'PTC', 'cid', 'CID', 'CN', 'hh', 'HH', 'PH']
n_cells = 4

time = [data['time'][i] for i in range(len(data['time'])) if 
        data['species'][i] == 0 and 
        data['cell'][i] == 0]

y = {}
for i,sp in enumerate(sp_names):
    for j in range(n_cells):
        y['%s_%d' % (sp,j)] = [data['conc'][k] for k in range(len(data['conc'])) if 
        data['species'][k] == i and 
        data['cell'][k] == j]

for name in sp_names:
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
    n = 0
    for i in range(2):
        for j in range(2):
            axes[i,j].plot(time, y['%s_%d' % (name,n)], lw=4, color='r', label='%s_%d' % (name, n))
            axes[i,j].set_ylim(-0.1 ,1.1)
            axes[i,j].legend(loc=0)
            if i==1:
                axes[i,j].set_xlabel('time')
            if j==0:
                axes[i,j].set_ylabel('scaled conc')            
            n += 1
    plt.tight_layout()

plt.show()
