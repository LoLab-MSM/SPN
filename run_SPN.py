from Model_SPN import model
from pysb.integrate import odesolve
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.cm import ScalarMappable

sp_names = ['en', 'EN', 'wg', 'IWG', 'EWG', 'ptc', 'PTC', 'ci', 'CI', 'CN', 'hh', 'HH', 'PH']
n_cells = 4

tspan = np.linspace(0, 1000, 501)
x = odesolve(model, tspan, verbose=True, integrator_options = {'atol' : 1e-6, 'rtol' : 1e-6})

EWG = [sum(x['EWG_%d_%d_obs' % (i+1,j+1)] for j in range(6)) for i in range(n_cells)]
PTC = [sum(x['PTC_%d_%d_obs' % (i+1,j+1)] for j in range(6)) for i in range(n_cells)]
HH =  [sum(x['HH_%d_%d_obs' % (i+1,j+1)] for j in range(6)) for i in range(n_cells)]
PH =  [sum(x['PH_%d_%d_obs' % (i+1,j+1)] for j in range(6)) for i in range(n_cells)]

names2_dict = {'EWG' : EWG, 'PTC' : PTC, 'HH'  : HH, 'PH'  : PH}

# create figure
fig1 = plt.figure(figsize=(20,20)) #figsize=(5.1, 15))
ax1 = fig1.add_subplot(111, aspect='equal')
ax1.axes.get_xaxis().set_visible(False)
ax1.axes.get_yaxis().set_visible(False)
ax1.axis([-3, 16, -27.5, 6.5])

cmap = plt.get_cmap('bwr') # blue-white-red colormap

# add cells (hexagons)
y_coord = 5
for name in sp_names:
    x_coord = 2
    for cell in range(n_cells):
        if name in names2_dict.keys():
            obs = names2_dict[name][cell]
        else:
            obs = x["%s_%d_obs" % (name, cell+1)]
        ax1.add_patch(
            patches.RegularPolygon(
                (x_coord, y_coord - (cell % 2)*np.cos(np.pi/6)),  # (x,y)
                6,          # number of vertices
                1,          # radius
                np.pi/6,
                ec='k',
                fc=cmap(obs[0])
            )
        )
        x_coord += 1.5
        ax1.annotate(name, xy=(-2.5, y_coord - np.cos(np.pi/6)/2), fontsize = 25)
    y_coord -= 2.5

# add current time 
time_text = ax1.text(8,4,'t = %g' % tspan[0], fontsize = 25, color='r')

# add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []
cbaxes = fig1.add_axes([0.58, 0.05, 0.03, 0.8]) # [left, bottom, width, height] -- fractions of plotting area between 0 and 1
cb = plt.colorbar(sm, cax=cbaxes)
cb.ax.tick_params(labelsize=20)
cb.ax.set_ylabel('Scaled concentration', fontsize=20, labelpad=-100)

plt.tight_layout()

def animate(t):
    time_text.set_text('t = %g' % tspan[t])
    i = 0
    for name in sp_names:
        for cell in range(n_cells):
            if name in names2_dict.keys():
                obs = names2_dict[name][cell]
            else:
                obs = x["%s_%d_obs" % (name, cell+1)]
            ax1.patches[i].set_fc(cmap(obs[t]))
            i += 1   
    return ax1.patches + [time_text]

ani = animation.FuncAnimation(fig1, animate, range(len(tspan)), interval=25, blit=True, repeat=False) 

print 'Saving animation...'
ani.save('anim.mp4', writer='ffmpeg')
# ani.save('anim.gif', writer='imagemagick')
print 'Done'

##### Plot time courses #####

# DATA FROM INGENUE 
data = np.genfromtxt('TEMP/SPN_Ingenue.txt', names = ['pset', 'time', 'cell', 'species', 'conc'])

y_time = [data['time'][i] for i in range(len(data['time'])) if 
        data['species'][i] == 0 and 
        data['cell'][i] == 0]

y = {}
for i,sp in enumerate(sp_names):
    for j in range(n_cells):
        y['%s_%d' % (sp,j+1)] = [data['conc'][k] for k in range(len(data['conc'])) if 
        data['species'][k] == i and 
        data['cell'][k] == j]

# Plot our results vs. Ingenue's
for name in sp_names:
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
    n = 1
    for i in range(2):
        for j in range(2):
            if name in names2_dict.keys():
                yvals = names2_dict[name][n-1]
            else:
                yvals = x['%s_%d_obs' % (name, n)]
            axes[i,j].plot(tspan, yvals, lw=4, color='r', label='%s_%d' % (name, n))
            axes[i,j].plot(y_time, y['%s_%d' % (name, n)], 'ob', mfc='None')
            axes[i,j].legend(loc=0)
            if i==1:
                axes[i,j].set_xlabel('time')
            if j==0:
                axes[i,j].set_ylabel('scaled conc')            
            n += 1
    plt.tight_layout()
#####

plt.show()
