# Cosmic-ray Astrophysics TOols (CATO)

A package containing a function for ionization cross-sections.

You can install this package as follows

**Shell Command**:
    ```markdown
    ```sh
    pip3 install git+https://github.com/vhmphan/cato.git
    ```
    ```


Here is an example for using the ionization cross section and gamma-ray cross section

**Python Code**:
    ```markdown
    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import cato.cato as ct


    # Testing electron capture cross-section
    Ep_p = 1000 # eV 
    sigma_ec_value = ct.func_sigma_ec(Ep_p) # cm^2

    print(f"The electron capture cross-section for Ep_p = {Ep_p} eV is {sigma_ec_value} cm^2")


    # Testing gamma-ray cross-section by plotting the local gamma-ray emissivity
    E=np.logspace(8,15,701)
    Eg=E

    d_sigma_g=ct.func_d_sigma_g(Tp,Eg) # Minh's code
    eps_nucl=ct.func_enhancement(Tp) 
    phi_LOC=np.zeros_like(Eg)
    for i in range(len(Eg)):
        phi_LOC[i]=trapezoid(eps_nucl*ct.func_jLocISM_p(E)*d_sigma_g[:,i],Tp) # eV^-1 s^-1

    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    ax.plot(Eg*1.0e-9,Eg**2*phi_LOC*1.0e-9,'g--',label='Local Emissivity')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$E_\gamma {\rm (GeV)}$',fontsize=fs)
    ax.set_ylabel(r'$E_\gamma \phi(E_\gamma) {\rm (GeV\,s^{-1})}$',fontsize=fs)
    for label_axd in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axd.set_fontsize(fs)
    ax.set_xlim(1.0e8,1.0e15)
    ax.legend(loc='lower left', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("fg_LocISM_emissivity.png")
    ```
    ```