import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

dens_data_day = np.genfromtxt('./density1d_day.out')
dens_data_night = np.genfromtxt('./density1d_night.out')
hannes_dens_data = np.genfromtxt('/home/rodney/Pictures/hannes_H_densities.csv', delimiter=',')
nagy_dens_data = np.genfromtxt('/home/rodney/Pictures/Nagy1990.csv', delimiter=',')

surface = 3.397e8

fig, ax = fig, ax = plt.subplots(1, 1)
plt.plot(dens_data_day[143:,0], dens_data_day[143:,1], label='This work')
plt.plot(hannes_dens_data[1:,0], hannes_dens_data[1:,1], label='Gr$\mathrm{\ddot{o}}$ller Poster')
plt.plot(nagy_dens_data[1:,0], nagy_dens_data[1:,1], label='Nagy 1990')
ax.xaxis.set_tick_params(top='on', direction='in', which='both')
ax.yaxis.set_tick_params(right='on', direction='in', which='both')
plt.yscale('log')
plt.xlim(0, 3000)
plt.ylim(0.01, 1000)
plt.xlabel('Altitude [km]')
plt.ylabel('Dayside density [$\mathrm{cm^{-3}}$]')
plt.title('Dayside Densities')
plt.legend()
plt.savefig('HotH_density_compare_day.svg')
#plt.show()
plt.close()

fig, ax = fig, ax = plt.subplots(1, 1)
plt.plot(dens_data_night[143:,0], dens_data_night[143:,1], label='This work')
plt.plot(hannes_dens_data[1:,0], hannes_dens_data[1:,1], label='Gr$\mathrm{\ddot{o}}$ller Poster')
plt.plot(nagy_dens_data[1:,0], nagy_dens_data[1:,1], label='Nagy 1990')
ax.xaxis.set_tick_params(top='on', direction='in', which='both')
ax.yaxis.set_tick_params(right='on', direction='in', which='both')
plt.yscale('log')
plt.xlim(0, 3000)
plt.ylim(0.1, 1000)
plt.xlabel('Altitude [km]')
plt.ylabel('Nightside density [$\mathrm{cm^{-3}}$]')
plt.title('Nightside Densities')
plt.legend()
plt.savefig('HotH_density_compare_night.svg')
#plt.show()
plt.close()

G = 6.6743e-11
jev = 1.602176487e-19
M_mars = 6.4185e23
k = 1.380649e-23
k_ev = 8.617333262e-5
m = 1.6735575e-27
T = 200.0

alt_list = [145, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000]
loss_rates = []
v_min = 0.01
v_max = 100
my_cmap = mpl.cm.get_cmap('inferno')
my_cmap.set_bad('r')
my_norm = LogNorm(vmin=v_min, vmax=v_max)

for j in alt_list:
    alt = 1e3*(j)
    r = 3397000.0 + alt
    v_esc = np.sqrt(2*G*M_mars/r)
    e_esc = (0.5*m*v_esc**2) / jev
    print('v_esc ' + str(j) + ' = ' + str(v_esc) + ' m/s')
    print('e_esc ' + str(j) + ' = ' + str(e_esc) + ' eV')
    
    EDF_day = np.genfromtxt('./EDF_day_' + str(j) + 'km.out')
    EDF_night = np.genfromtxt('./EDF_night_' + str(j) + 'km.out')
    EDF_day = np.clip(EDF_day, v_min, v_max)
    EDF_night = np.clip(EDF_night, v_min, v_max)

    surf = 2.0*np.pi * (surface + 1e5*(j+1))**2
    vol = (2.0/3.0)*np.pi * ((surface + 1e5*(j+1))**3 - (surface + 1e5*j)**3)
    
    plt.imshow(EDF_day, cmap=my_cmap, interpolation='none', origin='lower', extent=[-1,1,0,10], aspect='.2', norm=my_norm)
    plt.xlabel('$\mu$')
    plt.ylabel('Energy [eV]')
    #plt.xlim(700, 1000)
    #plt.ylim(0, 300)
    plt.colorbar(label='Distribution [$\mathrm{cm^{-3} eV^{-1} \mu^{-1}}$]')
    plt.savefig('EDF_day_map_' + str(j) + 'km.svg')
    plt.close()

    plt.imshow(EDF_night, cmap=my_cmap, interpolation='none', origin='lower', extent=[-1,1,0,10], aspect='.2', norm=my_norm)
    plt.xlabel('$\mu$')
    plt.ylabel('Energy [eV]')
    #plt.xlim(700, 1000)
    #plt.ylim(0, 300)
    plt.colorbar(label='Distribution [$\mathrm{cm^{-3} eV^{-1} \mu^{-1}}$]')
    plt.savefig('EDF_night_map_' + str(j) + 'km.svg')
    plt.close()
    
    EDF_day_summed = np.sum(EDF_day, axis=1)#np.zeros(201)
    EDF_night_summed = np.sum(EDF_night, axis=1)#np.zeros(201)
    EDF_day_summed = EDF_day_summed*0.01
    EDF_night_summed = EDF_night_summed*0.01

    '''
    for i in range(201):
        for k in range(100):
            EDF_day_summed[i] += EDF_day[i, k]*(.05*.01)*(-1.0+(k*.01))
            EDF_night_summed[i] += EDF_night[i, k]*(.05*.01)*(-1.0+(k*.01))
        for k in range(101):
            EDF_day_summed[i] += EDF_day[i, k+100]*(.05*.01)*(0.0+(k*.01))
            EDF_night_summed[i] += EDF_night[i, k+100]*(.05*.01)*(0.0+(k*.01))
    '''
    
    x = np.arange(0, 10.05, 0.05)
    plt.plot(x, EDF_day_summed, label='Dayside EDF at ' + str(j) + ' km')
    plt.plot(x, EDF_night_summed, label='Nightside EDF at ' + str(j) + ' km')

    e = np.linspace(0, 10, 201)
    p_e = 2*np.sqrt(e / np.pi) * (1 / (k_ev*T))**(3/2) * np.exp(-e/(k_ev*T))
    sum = np.sum(p_e)
    p_e = (p_e/sum) * np.sum(EDF_day_summed)
    plt.plot(e, p_e, label='MB at T = 200K')

    plt.xlabel('Energy [eV]')
    plt.ylabel('Distribution [$\mathrm{cm^{-3} eV^{-1}}$]')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 10)
    plt.ylim(0.01, 100)
    plt.title('Energy Distribution Function')
    plt.legend()
    plt.savefig('EDF_daynight' + str(j) + 'km.svg')
    #plt.show()
    plt.close()
    
    for i in range(201):
        EDF_day_summed[i] = EDF_day_summed[i] * (100.0*np.sqrt((2.0*i*0.05*jev)/m))
        EDF_night_summed[i] = EDF_night_summed[i] * (100.0*np.sqrt((2.0*i*0.05*jev)/m))
    
    sum_day = 0.0
    sum_night = 0.0
    for i in range(len(EDF_day_summed)):
        if (x[i]>e_esc):
            sum_day += 0.05*EDF_day_summed[i]
            sum_night += 0.05*EDF_night_summed[i]
    
    day = sum_day*surf
    night = sum_night*surf
    total_loss = day + night
    loss_rates.append(total_loss)

    print('day escape integrated over eV times surface: ' + str(day))
    print('night escape integrated over eV times surface: ' + str(night))    
    print('total loss rate: ' + str(total_loss) + '\n')
    
loss_rates2 = np.genfromtxt('./loss_rates.out')
    
plt.plot(alt_list, loss_rates, label='Loss rate from integrated EDFs')
plt.plot(alt_list, loss_rates2[0:,1], label='Loss rate directly from code')
#plt.axhline(y=2.10067e25, linestyle='--', label='Fraction lost at 60,000 km')
#plt.axhline(y=1.75384e25, linestyle='--', label='Fraction lost at 5,000 km')
#plt.axvline(x=14000, linestyle='--', label='14000 km')
#plt.axvline(x=16000, linestyle='--', label='16000 km')
plt.xlabel('Altitude [km]')
plt.ylabel('Loss rate [$\mathrm{s^{-1}}$]')
#plt.ylim(1.7e25, 1.8e25)
plt.legend()
plt.savefig('EDF_fraction_compare.svg')
#plt.show()
plt.close()
    
  
