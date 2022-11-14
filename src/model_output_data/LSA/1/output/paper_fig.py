import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

G = 6.6743e-11
jev = 1.602176487e-19
M_mars = 6.4185e23
k = 1.380649e-23
k_ev = 8.617333262e-5
m = 1.6735575e-27
T = 200.0

fig = plt.figure()
fig.set_figheight(4)
fig.set_figwidth(14)

ax1 = plt.subplot2grid(shape=(2, 7), loc=(0, 0), colspan=2, rowspan=2)
ax2 = plt.subplot2grid(shape=(2, 7), loc=(0, 2))
ax3 = plt.subplot2grid(shape=(2, 7), loc=(0, 3))
ax4 = plt.subplot2grid(shape=(2, 7), loc=(0, 4))
ax5 = plt.subplot2grid(shape=(2, 7), loc=(0, 5))
ax6 = plt.subplot2grid(shape=(2, 7), loc=(1, 2))
ax7 = plt.subplot2grid(shape=(2, 7), loc=(1, 3))
ax8 = plt.subplot2grid(shape=(2, 7), loc=(1, 4))
ax9 = plt.subplot2grid(shape=(2, 7), loc=(1, 5))

dens_data_day = np.genfromtxt('./density1d_day.out')

v_min = 0.01
v_max = 100
my_cmap = mpl.cm.get_cmap('inferno')
my_cmap.set_bad('k')
my_norm = LogNorm(vmin=v_min, vmax=v_max)

EDF_day1 = np.genfromtxt('./EDF_day_145km.out')
im1 = ax6.imshow(EDF_day1.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day2 = np.genfromtxt('./EDF_day_300km.out')
im2 = ax7.imshow(EDF_day2.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day3 = np.genfromtxt('./EDF_day_600km.out')
im3 = ax8.imshow(EDF_day3.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
EDF_day4 = np.genfromtxt('./EDF_day_1200km.out')
im4 = ax9.imshow(EDF_day4.T, cmap=my_cmap, interpolation='none', origin='lower', extent=[0,10,-1,1], aspect='5', norm=my_norm)
#ax2.xlabel('$\mu$')
#ax2.ylabel('Energy [eV]')

markers_on = [145, 300, 600, 1200]

ax1.plot(dens_data_day[:,1], dens_data_day[:,0], markevery=markers_on, marker='D', mfc='r', mec='r')
ax1.xaxis.set_tick_params(top='on', direction='in', which='both')
ax1.yaxis.set_tick_params(right='on', direction='in', which='both')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylim(90, 5000)
ax1.set_xlim(1, 200)
ax1.set_ylabel('Altitude [km]')
ax1.set_xlabel('Density [$\mathrm{cm^{-3}}$]')
ax1.set_title('Dayside Hot H Density')


EDF_day_summed1 = np.sum(EDF_day1, axis=1)#np.zeros(201)
EDF_day_summed2 = np.sum(EDF_day2, axis=1)#np.zeros(201)
EDF_day_summed3 = np.sum(EDF_day3, axis=1)#np.zeros(201)
EDF_day_summed4 = np.sum(EDF_day4, axis=1)#np.zeros(201)

EDF_day_summed1 = EDF_day_summed1*(.01)
EDF_day_summed2 = EDF_day_summed2*(.01)
EDF_day_summed3 = EDF_day_summed3*(.01)
EDF_day_summed4 = EDF_day_summed4*(.01)


'''
for i in range(201):
    for k in range(100):
        EDF_day_summed1[i] += EDF_day1[i, k]*(.05*.01)*(-1.0+(k*.01))
        EDF_day_summed2[i] += EDF_day2[i, k]*(.05*.01)*(-1.0+(k*.01))
        EDF_day_summed3[i] += EDF_day3[i, k]*(.05*.01)*(-1.0+(k*.01))
        EDF_day_summed4[i] += EDF_day4[i, k]*(.05*.01)*(-1.0+(k*.01))
    for k in range(101):
        EDF_day_summed1[i] += EDF_day1[i, k+100]*(.05*.01)*(0.0+(k*.01))
        EDF_day_summed2[i] += EDF_day2[i, k+100]*(.05*.01)*(0.0+(k*.01))
        EDF_day_summed3[i] += EDF_day3[i, k+100]*(.05*.01)*(0.0+(k*.01))
        EDF_day_summed4[i] += EDF_day4[i, k+100]*(.05*.01)*(0.0+(k*.01))
'''

x = np.arange(0, 10.05, 0.05)
ax2.plot(x, EDF_day_summed1, label='EDF')
ax3.plot(x, EDF_day_summed2, label='EDF')
ax4.plot(x, EDF_day_summed3, label='EDF')
ax5.plot(x, EDF_day_summed4, label='EDF')
e = np.linspace(0, 10, 201)
p_e = 2*np.sqrt(e / np.pi) * (1 / (k_ev*T))**(3/2) * np.exp(-e/(k_ev*T))
sum = np.sum(p_e)
p_e1 = (p_e/sum) * np.sum(EDF_day_summed1)
#ax2.plot(e, p_e1, label='MB at 200K')
p_e2 = (p_e/sum) * np.sum(EDF_day_summed2)
#ax3.plot(e, p_e2, label='MB at 200K')
p_e3 = (p_e/sum) * np.sum(EDF_day_summed3)
#ax4.plot(e, p_e3, label='MB at 200K')
p_e4 = (p_e/sum) * np.sum(EDF_day_summed4)
#ax5.plot(e, p_e4, label='MB at 200K')

#ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0, 10)
ax2.set_ylim(0.1, 500)
#ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim(0, 10)
ax3.set_ylim(0.1, 500)
#ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlim(0, 10)
ax4.set_ylim(0.1, 500)
#ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_xlim(0, 10)
ax5.set_ylim(0.1, 500)

ax2.set_title('145 km')
ax3.set_title('300 km')
ax4.set_title('600 km')
ax5.set_title('1200 km')
ax2.set_ylabel('Dist. [$\mathrm{cm^{-3} eV^{-1}}$]')
ax6.set_ylabel('$\mu$')
ax6.set_xlabel('eV')
ax7.set_xlabel('eV')
ax8.set_xlabel('eV')
ax9.set_xlabel('eV')

plt.tight_layout()
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
#plt.setp(ax3.get_yticklabels(), visible=False)
#plt.setp(ax4.get_yticklabels(), visible=False)
#plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax7.get_yticklabels(), visible=False)
plt.setp(ax8.get_yticklabels(), visible=False)
plt.setp(ax9.get_yticklabels(), visible=False)

ax10 = plt.subplot2grid(shape=(2, 7), loc=(0, 6))
ax11 = plt.subplot2grid(shape=(2, 7), loc=(1, 6))
ax10.set_axis_off()
ax11.set_axis_off()
fig.colorbar(im1, ax=[ax10, ax11], location='left', label='Distribution [$\mathrm{cm^{-3} eV^{-1} \mu^{-1}}$]')
plt.savefig('fig1.png')
plt.close()


