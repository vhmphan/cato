# Cosmic-ray Astrophysics TOols (CATO)

A package containing a function for ionization cross-sections.

**Installation**: You can install this package as follows
```sh
pip3 install git+https://github.com/vhmphan/cato.git
```

**Upgrade**: You can also upgrade this package as follows
```sh
pip3 install --upgrade git+https://github.com/vhmphan/cato.git
```

**Example code**: Here is an example for using the ionization cross section and gamma-ray cross section
```python
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
import cato.cato as ct


# Plotting ionization cross-sections 
E=np.logspace(1.0,11.0,100)

IonXS = ct.Cross_Section_Ion(E,Es)
sigma_ec=IonXS.func_sigma_ec()*1.0e17 # 1.0e-17 cm^2
sigma_p=IonXS.func_sigma_p()*1.0e17 # 1.0e-17 cm^2
sigma_e=IonXS.func_sigma_e()*1.0e17 # 1.0e-17 cm^2

fig=plt.figure(figsize=(10, 8))
ax=plt.subplot(111)

ax.plot(E,sigma_ec,'b--',linewidth=3,label=r'$\sigma_{\rm ec}$')
ax.plot(E,sigma_e,'g-.',linewidth=3,label=r'$\sigma^{\rm ion}_{\rm e}$')
ax.plot(E,sigma_p,'r-',linewidth=3,label=r'$\sigma^{\rm ion}_{\rm p}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_xlabel(r'$E {\rm (eV)}$',fontsize=fs)
ax.set_ylabel(r'$\sigma {\rm (10^{-17}\,cm^{2})}$',fontsize=fs)
for label_axd in (ax.get_xticklabels() + ax.get_yticklabels()):
    label_axd.set_fontsize(fs)
ax.set_xlim(1.0e1,1.0e11)
ax.set_ylim(1.0e-3,2.0e2)
ax.legend(loc='upper right', prop={"size":22})
ax.grid(linestyle='--')

plt.savefig("fg_ionizaton_cross-section.png")

# Load the data files for secondary ionization provided by the package
data_p=ct.load_data_file('phi_p.dat')
data_e=ct.load_data_file('phi_e.dat')

E_data_p, phi_data_p=data_p[:, 0], data_p[:, 1]
E_data_e, phi_data_e=data_e[:, 0], data_e[:, 1]

# Plotting secondary ionization functions
fig=plt.figure(figsize=(10, 8))
ax=plt.subplot(111)

ax.plot(E_data_p,phi_data_p,'r-',linewidth=3,label=r'${\rm Proton}$')
ax.plot(E_data_e,phi_data_e,'g-.',linewidth=3,label=r'${\rm Electron}$')

ax.set_xscale('log')
ax.legend()
ax.set_xlabel(r'$E {\rm (eV)}$',fontsize=fs)
ax.set_ylabel(r'$\sigma {\rm (10^{-17}\,cm^{2})}$',fontsize=fs)
for label_axd in (ax.get_xticklabels() + ax.get_yticklabels()):
    label_axd.set_fontsize(fs)
ax.set_xlim(1.0e2,1.0e11)
ax.set_ylim(0.0,1.2)
ax.legend(loc='upper right', prop={"size":22})
ax.grid(linestyle='--')

plt.savefig("fg_secondary-ionizaton.png")

# Testing gamma-ray cross-section by plotting the local gamma-ray emissivity
Tp=np.logspace(8,15,701)
Eg=Tp

RadXS=ct.Cross_Section_Rad(Tp,Eg)
d_sigma_g=RadXS.func_d_sigma_g()
eps_nucl=RadXS.func_enhancement()

phi_LOC=np.zeros_like(Eg)
for i in range(len(Eg)):
    phi_LOC[i]=trapezoid(eps_nucl*ct.func_jLocISM_p(Tp)*d_sigma_g[:,i],Tp) # eV^-1 s^-1

fs=22 # Fontsize for plot

# Plot gamma-ray emissivity
fig=plt.figure(figsize=(10, 8))
ax=plt.subplot(111)

ax.plot(Eg*1.0e-9,Eg**2*phi_LOC*1.0e-9,'g--',linewidth=5,label=r'${\rm Local\, Emissivity}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E_\gamma {\rm (GeV)}$',fontsize=fs)
ax.set_ylabel(r'$E_\gamma \phi(E_\gamma) {\rm (GeV\,s^{-1}\,sr^{-1})}$',fontsize=fs)
for label_ax in (ax.get_xticklabels() + ax.get_yticklabels()):
    label_ax.set_fontsize(fs)
ax.set_xlim(1.0e-1,1.0e6)
ax.set_ylim(1.0e-30,1.0e-26)
ax.legend(loc='upper right', prop={"size":fs})
ax.grid(linestyle='--')

plt.savefig("fg_LocISM_emissivity.png")
```

**References**: Ionization cross-sections are from [Padovani et al. 2009](https://ui.adsabs.harvard.edu/abs/2009A%26A...501..619P/abstract) and [Krause et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ICRC...34..518K/abstract). Gamma-ray cross-sections are from [Kafexhiu et al. 2014](https://ui.adsabs.harvard.edu/abs/2014PhRvD..90l3014K/abstract).
