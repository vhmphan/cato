# Cosmic-ray Astrophysics TOols (CATO)

A package containing a function for ionization cross-sections.

**Installation**: You can install this package as follows
    ```
    pip3 install git+https://github.com/vhmphan/cato.git
    ```

**Installation**: You can also upgrade this package as follows
    ```
    pip3 install git+https://github.com/vhmphan/cato.git
    ```

**Example code**: Here is an example for using the ionization cross section and gamma-ray cross section
    ```python
    import numpy as np
    from scipy.integrate import trapezoid
    import matplotlib.pyplot as plt
    import cato.cato as ct


    # Testing electron capture cross-section
    Ep_p = 1000 # eV 
    sigma_ec_value = ct.func_sigma_ec(Ep_p) # cm^2

    print(f"The electron capture cross-section for Ep_p = {Ep_p} eV is {sigma_ec_value} cm^2")


    # Testing gamma-ray cross-section by plotting the local gamma-ray emissivity
    Tp=np.logspace(8,15,701)
    Eg=Tp

    d_sigma_g=ct.func_d_sigma_g(Tp,Eg) # Minh's code
    eps_nucl=ct.func_enhancement(Tp) 
    phi_LOC=np.zeros_like(Eg)
    for i in range(len(Eg)):
        phi_LOC[i]=trapezoid(eps_nucl*ct.func_jLocISM_p(Tp)*d_sigma_g[:,i],Tp) # eV^-1 s^-1

    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    ax.plot(Eg*1.0e-9,Eg**2*phi_LOC*1.0e-9,'g--',linewidth=5,label=r'${\rm Local\, Emissivity}$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$E_\gamma {\rm (GeV)}$',fontsize=fs)
    ax.set_ylabel(r'$E_\gamma \phi(E_\gamma) {\rm (GeV\,s^{-1}\,sr^{-1})}$',fontsize=fs)
    for label_axd in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axd.set_fontsize(fs)
    ax.set_xlim(1.0e-1,1.0e6)
    ax.set_ylim(1.0e-30,1.0e-26)
    ax.legend(loc='upper right', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("fg_LocISM_emissivity.png")
    ```
    ```