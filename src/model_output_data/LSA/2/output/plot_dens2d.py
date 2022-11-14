import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

coldens_data = np.genfromtxt('column_density_day.out')
dens_2d_data = np.genfromtxt('density2d.out')

plt.plot(coldens_data[:,0], coldens_data[:,1])
plt.xlim(0, 32000)
plt.ylim(1e6, 1e14)
plt.yscale('log')
#plt.xscale('log')
plt.show()

v_min = 1
v_max = 1e6
my_cmap = mpl.cm.get_cmap('inferno')
my_cmap.set_bad('k')
my_norm = LogNorm()#vmin=v_min, vmax=v_max)

plt.imshow(dens_2d_data, cmap=my_cmap, interpolation='none', norm=my_norm)
plt.colorbar()
plt.show()
#plt.savefig('density2d.png')
