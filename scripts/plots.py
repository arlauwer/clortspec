import numpy as np
import matplotlib.pyplot as plt

# column 1: wavelength; E (keV)
# column 2: total flux; F_E (1/s/cm2/keV)
# column 3: transparent flux; F_E (1/s/cm2/keV)
# column 4: direct primary flux; F_E (1/s/cm2/keV)
# column 5: scattered primary flux; F_E (1/s/cm2/keV)
# column 6: direct secondary flux; F_E (1/s/cm2/keV)
# column 7: scattered secondary flux; F_E (1/s/cm2/keV)
# column 8: transparent secondary flux; F_E (1/s/cm2/keV)
name = 'model_i90_sed.dat'

phot = np.loadtxt('phot/' + name)
ref = np.loadtxt('ref/' + name)
zone = np.loadtxt('zone/' + name)

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

def plot(ax, id, label, **kwargs):
	ax.plot(phot[:, 0], phot[:, id], label='phot ' + label, **kwargs)
	ax.plot(ref[:, 0], ref[:, id], label='ref ' + label, **kwargs)
	ax.plot(zone[:, 0], zone[:, id], label='zone ' + label, **kwargs)

plot(ax, 1, 'total', alpha=0.4)

ax.set_xlabel('wavelength (keV)')
ax.set_ylabel('flux (1/s/cm2/keV)')
ax.legend()

ax.set_xscale('log')
ax.set_yscale('log')

plt.show()