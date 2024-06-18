import numpy as np
import scipy as sp

mp=938.272e6 # eV
me=0.510998e6 # eV

IH2=15.603 # eV
RH=13.6 # eV
w2=RH/(4.0*IH2)  	
alpha=0.87 # Parameter of the Rudd model (to be distinguish with the fine structure constant) 
alpha_f=1/137.035999074 # Fine structure constant
a0=0.52917721092*pow(10.0,-8) # cm
N=2.0 
Sp=4.0*np.pi*pow(a0,2)*N*pow(RH/IH2,2) 
lambda_pe=(2.0*mp+2.0*me)/me 
Se=4.0*np.pi*pow(a0,2)*N*pow(alpha_f,4) 


####################################################################################################################################
# Proton ionization cross-sections
####################################################################################################################################


# Electron capture cross-section
def func_sigma_ec(Ep_p):
    
    d_0=-52.793928 
    d_1=41.219156 
    d_2=-17.304947 
    d_3=3.292795 
    d_4=-0.238372 
    
    x=Ep_p/((mp/me)*IH2) 
    
    sigma_ec=pow(10.0,(d_0+d_1*np.log10(Ep_p)+d_2*pow(np.log10(Ep_p),2)+d_3*pow(np.log10(Ep_p),3)+d_4*pow(np.log10(Ep_p),4))) 
    
    return sigma_ec # cm^2


# Auxiliary functions for proton ionization differential cross-section 
def func_F1(v):

	A1=0.8 
	B1=2.9 
	C1=0.86 
	D1=1.48 
	E1=7.0 
	
	L1=C1*pow(v,D1)/(1.0+E1*pow(v,D1+4)) 
	H1=A1*np.log(1.0+pow(v,2))/(pow(v,2)+B1/pow(v,2.0)) 

	F1=L1+H1 
	return F1 


# Auxiliary function for proton ionization differential cross-section 
def func_F2(v):
	
	A2=1.06 
	B2=4.2 
	C2=1.39 
	D2=0.48 

	L2=C2*pow(v,D2) 
	H2=(A2/pow(v,2))+(B2/pow(v,4)) 

	F2=L2*H2/(L2+H2) 	
	return F2 


# Differential cross-section of Rudd model
def func_d_sigma_p_E(Ep_p, Ee_s):

	p=np.sqrt(pow(Ep_p+mp,2)-pow(mp,2)) 
	w=Ee_s/IH2 
		
	v_rel=p/np.sqrt(pow(p,2)+pow(mp,2)) 
	v=np.sqrt(2.0*mp*me/(IH2*(2.0*mp+2.0*me)))*v_rel 
	T=Ep_p/lambda_pe 
	wc=((4.0*Ep_p-2.0*np.sqrt(IH2*T))/IH2)-w2 

	d_sigma_E=(Sp/IH2)*(func_F1(v)+func_F2(v)*w)*pow(1.0+w,-3)/(1.0+np.exp(alpha*(w-wc)/v)) 

	return d_sigma_E 


# # Proton ionization cross-section
# def func_sigma_p(Ep_p):# Proton ionization cross-section
	
# 	p, sigma_p, w, dw, wmax, wc, T, Qmax, v, v_rel, gamma 
# 	p=np.sqrt(pow(Ep_p+mp,2)-pow(mp,2)) # eV
# 	w=0.0 

# 	v_rel=p/np.sqrt(pow(p,2)+pow(mp,2)) # Dimensionless -> p in eV and Ep_tot in eV
# 	gamma=1.0/(np.sqrt(1.0-pow(v_rel,2))) 
# 	Qmax=2.0*pow(gamma,2)*me*pow(v_rel,2)/(1.0+2.0*gamma*(me/mp)+pow(me/mp,2)) # Qmax is the maximum energy transferred-> eV 
# 	v=np.sqrt(2.0*mp*me/((2.0*mp+2.0*me)*I))*v_rel # Dimensionless -> reduced speed 

# 	T=Ep_p/lambda 
# 	wc=((4.0*Ep_p-2.0*np.sqrt(I*T))/I)-w2 
# 	wmax=Qmax/I 
		
# 	sigma_p=0.0 
#     w=0.0 

# 	while(w<wmax):				
# 		dw=min(0.001*w,wmax-w) 

# 		if(w==0.0):
# 			dw=0.001 
# 			
			
# 		sigma_p+=dw*Sp*(def func_F1(v)+def func_F2(v)*w)*pow(1.0+w,-3)/(1.0+exp(alpha*(w-wc)/v)) 	
# 		w+=dw 
# 

# 	return sigma_p 
#


# Proton ionization cross-section
def func_sigma_p(Ep_p): 

    A=0.314
    B=1.01
    C=3.24
    D=94.4 
    sigma0=6.51e-14 # cm^2 eV^2
    vrel=np.sqrt(pow(Ep_p+mp,2)-mp*mp)/(Ep_p+mp) 
    
	# Note that we use the formula as given in Padovani et al. A&A 2009 (Eq. 5) but we adopt the relativistic 
	# correction as suggested in Plante et al. NJP 2008 (Eq. 11, to replace x by its relativistic version).
    x=me*pow(vrel,2)/(2.0*IH2) 
    
    sigma_p=A*sigma0*2.0*pow(x,B)*np.log(1.0+D*x)/(C+pow(x,B+1)) 
    sigma_p*=1.0/pow(IH2,2) 

    return sigma_p # cm^2


####################################################################################################################################
# Electron ionization cross-section
####################################################################################################################################


# Auxiliary function for electron ionization differential cross-section 
def func_df_dw(w):

	c_df, d_df, e_df, f_df 
	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 

	df_dw=(c_df/pow(1.0+w,3))+(d_df/pow(1.0+w,4))+(e_df/pow(1.0+w,5))+(f_df/pow(1.0+w,6)) 
	return df_dw 


# Auxiliary function for electron ionization differential cross-section 
def func_D(t):

	c_df, d_df, e_df, f_df 
	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 

	D=0.0 
	u=(t+1.0)/2.0 
	D+=-c_df*(pow(u,-3)-1.0)/3.0 
	D+=-d_df*(pow(u,-4)-1.0)/4.0 
	D+=-e_df*(pow(u,-5)-1.0)/5.0 
	D+=-f_df*(pow(u,-6)-1.0)/6.0 
	D*=1.0/N 

	return D 


# Auxiliary function for electron ionization differential cross-section 
def func_d_sigma_e_E(Ee_p, Ee_s):# Differential cross-section of electron (RBEB binary-encounter-Bethe)

	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 
	N_i=(c_df/2.0)+(d_df/2.0)+(e_df/2.0)+(f_df/2.0) 	

	U_H2=39.603 # Average kinetic energy of orbital electron
	B_H2=IH2 # Binding energy
	u_primed=U_H2/me 
	b_primed=B_H2/me 
	beta_u=np.sqrt(1.0-(1.0/pow(1.0+u_primed,2))) 
	beta_b=np.sqrt(1.0-(1.0/pow(1.0+b_primed,2))) 

	t=Ee_p/IH2 
	t_primed=Ee_p/me 
	beta_t=np.sqrt(1.0-(1.0/pow(1.0+t_primed,2))) 	
	w=Ee_s/IH2 

	d_sigma_E=0.0 	
	d_sigma_E+=(((N_i/N)-2)/(t+1.0))*((1.0/(1+w))+(1.0/(t-w)))*((1.0+2.0*t_primed)/(pow(1.0+(t_primed/2.0),2))) 
	d_sigma_E+=(2.0-(N_i/N))*((1.0/pow(1.0+w,2))+(1.0/pow(t-w,2))+(pow(b_primed,2)/pow(1.0+(t_primed/2.0),2))) 
	d_sigma_E+=(1.0/(N*(1.0+w)))*func_df_dw(w)*(np.log(pow(beta_t,2)/(1.0-pow(beta_t,2)))) 
	d_sigma_E+=(1.0/(N*(1.0+w)))*func_df_dw(w)*(-pow(beta_t,2)-np.log(2.0*b_primed)) 
	d_sigma_E*=((Se/(2.0*b_primed))/(pow(beta_t,2)+pow(beta_u,2)+pow(beta_b,2))) 
	d_sigma_E*=(1.0/IH2) 	

	return d_sigma_E 


# Electron ionization cross-section
def func_sigma_e(Ee_p):

	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 
	N_i=(c_df/2.0)+(d_df/2.0)+(e_df/2.0)+(f_df/2.0) 

	U_H2=39.603 # Average kinetic energy of orbital electron
	B_H2=IH2 # Binding energy
	u_primed=U_H2/me 
	b_primed=B_H2/me 
	beta_u=np.sqrt(1.0-(1.0/pow(1.0+u_primed,2))) 
	beta_b=np.sqrt(1.0-(1.0/pow(1.0+b_primed,2))) 	

	t=Ee_p/IH2 
	t_primed=Ee_p/me 
	beta_t=np.sqrt(1.0-(1.0/pow(1.0+t_primed,2))) 	

	sigma_e=0.0 
	sigma_e+=func_D(t)*(np.log(pow(beta_t,2)/(1.0-pow(beta_t,2)))-pow(beta_t,2)-np.log(2.0*b_primed)) 
	sigma_e+=(2.0-(N_i/N))*(1.0-(1.0/t)-(np.log(t)/(t+1))*((1.0+2.0*t_primed)/pow(1.0+(t_primed/2.0),2))) 
	sigma_e+=(2.0-(N_i/N))*((pow(b_primed,2)/pow(1.0+(t_primed/2.0),2))*((t-1.0)/2.0)) 
	sigma_e*=((Se/(2.0*b_primed))/(pow(beta_t,2)+pow(beta_u,2)+pow(beta_b,2))) 
	
	return sigma_e 


####################################################################################################################################
# Cosmic-ray spectra in the local ISM
####################################################################################################################################


# Fit of the Voyager spectrum for protons 
def func_jLocISM_p(E):
    
    C=1.882e-9 
    alpha=0.129056 
    beta=2.82891 
    E0=624.5e6 
    
    f0=12.5e-10*pow(E/1.0e6,0.35)/((1.0+pow(E/80.0e6,1.3))*pow(1.0+pow(E/2.2e9,1.9/2.1),2.1)) 

    return f0 # eV^-1 cm^-2 s^-1 sr^-1


# Fit of the Voyager spectrum for protons 
def func_jLocISM_e(E):
    
    C=4.658e-7 
    alpha=-1.236 
    beta=2.033 
    E0=736.2e6
    
    f0=C*(E/1.0e6)**alpha*(1.0+(E/E0))**(-beta)

    return f0 # eV^-1 cm^-2 s^-1 sr^-1


####################################################################################################################################
# Function for interpolation on log-log scale
####################################################################################################################################

# Function for interpolation
def log_interp1d(xx, yy, kind='linear'):

    logx=np.log10(xx)
    logy=np.log10(yy)
    lin_interp=sp.interpolate.interp1d(logx,logy,kind=kind)
    log_interp=lambda zz: np.power(10.0,lin_interp(np.log10(zz)))

    return log_interp


meCGS=9.10938356e-28 # g
mpCGS=1.67262192e-24 # g
qeCGS=4.8032e-10 # CGS unit -> Electric charge of proton

kB=8.617333262145e-5 # eV/K
sigmaT=6.6524e-25 # cm^-2 -> Thompson cross-section
sigma_sb=3.5394474508e7 # erg cm^-2 s^-1 K^-4
hP=4.1357e-15 # eV s
Ktomec2=1.6863699549e-10 # -> To convert temperture T from K/kB to eV/me 

mpi=134.9766e6 # eV
Tpth=2.0*mpi+(pow(mpi,2)/(2.0*mp)) # eV


############################################################################################
# Gamma rays from proton-proton interaction -> Kafexhiu et al. 2014
############################################################################################


# Dimensionless Breit-Wigner distribution.
def func_fBW(sqrt_s):
# Eq. 4.

    M_res=1.1883*pow(10.0,9) # eV
    Gamma_res=0.2264*pow(10.0,9) # eV
    gamma=np.sqrt(pow(M_res,2)*(pow(M_res,2)+pow(Gamma_res,2))) 
    K=np.sqrt(8.0)*M_res*Gamma_res*gamma/(np.pi*np.sqrt(pow(M_res,2)+gamma)) 
    
    fBW=mp*K/(pow(pow(sqrt_s-mp,2)-pow(M_res,2),2)+pow(M_res*Gamma_res,2)) 
    
    return fBW


# Cross-section for p+p -> p+p+pi0.
def func_sigma_1pi(Tp):
# Eq. 2,
# Tp (eV) -> kinetic energy of CR proton.

    sigma_0=7.66e-3 
    sqrt_s=np.sqrt(2.0*mp*(Tp+2.0*mp)) 
    eta=np.sqrt(pow(pow(sqrt_s,2)-pow(mpi,2)-4.0*pow(mp,2),2)-16.0*pow(mpi*mp,2))/(2.0*mpi*sqrt_s) 
    
    sigma_1pi=sigma_0*pow(eta,1.95)*(1.0+eta+pow(eta,5))*pow(func_fBW(sqrt_s),1.86) 
    
    return sigma_1pi # m


# Cross-section for p+p -> p+p+2pi0.
def func_sigma_2pi(Tp):
# Eq. 5,
# Tp (eV) -> kinetic energy of CR proton.
       
    mask1=(Tp<0.56e9)
    mask2=(Tp>=0.56e9)

    sigma_2pi=np.zeros_like(Tp)

    sigma_2pi[mask1]=0.0 
    sigma_2pi[mask2]=5.7/(1.0+np.exp(-9.3*(Tp[mask2]*1.0e-9-1.4))) 
    
    return sigma_2pi # m


# Proton-proton total inelastic cross-section. 
def func_sigma_inel(Tp):
# Eq. 1,
# Tp (eV) -> kinetic energy of CR proton.

    sigma_inel=30.7-0.96*np.log(Tp/Tpth)+0.18*pow(np.log(Tp/Tpth),2) 
    sigma_inel*=pow(1.0-pow(Tpth/Tp,1.9),3) 
    sigma_inel[sigma_inel<0.0]=0.0
    
    return sigma_inel # m


# Average pi0 multiplicity.
def func_npi(Tp):
# Eq. 7 and Table 4 (GEANT 4 model),
# Tp (eV) -> kinetic energy of CR proton.

    a1=0.728 
    a2=0.596 
    a3=0.491 
    a4=0.2503 
    a5=0.117 
    
    Qp=(Tp-Tpth)/mp 
    xip=(Tp-3.0e9)/mp 

    mask1=(Tp>=1.0e9) & (Tp<5.0e9)
    mask2=(Tp>=5.0e9)

    npi=np.zeros_like(Tp)

    npi[mask1]=-6.0e-3+0.237*Qp[mask1]-0.023*pow(Qp[mask1],2) 
    npi[mask2]=a1*pow(xip[mask2],a4)*(1.0+np.exp(-a2*pow(xip[mask2],a5)))*(1.0-np.exp(-a3*pow(xip[mask2],0.25))) 
    
    return npi


# Pi0 production cross-section.
def func_sigma_pi(Tp):
# See paragraph above Table 4,
# Tp (eV) -> kinetic energy of CR proton.
    
    mask1=(Tp>=Tpth) & (Tp<2.0e9)
    mask2=(Tp>=2.0e9)

    sigma_pi=np.zeros_like(Tp)

    sigma_pi[mask1]=func_sigma_1pi(Tp[mask1])+func_sigma_2pi(Tp[mask1]) 
    sigma_pi[mask2]=func_sigma_inel(Tp[mask2])*func_npi(Tp[mask2]) 
    
    return sigma_pi # m


# Complementary function for the differential cross-section.
def func_Amax(Tp):
# Eq. 12 and Table 7 (GEANT 4 model),
# Tp (eV) -> kinetic energy of CR proton.

    sqrt_s=np.sqrt(2.0*mp*(Tp+2.0*mp)) 
    gamma_CM=(Tp+2.0*mp)/sqrt_s 
    beta_CM=np.sqrt(1.0-pow(gamma_CM,-2)) 
    Epi_CM=(pow(sqrt_s,2)-4.0*pow(mp,2)+pow(mpi,2))/(2.0*sqrt_s) 
    Ppi_CM=np.sqrt(pow(Epi_CM,2)-pow(mpi,2)) 
    Epi_max=gamma_CM*(Epi_CM+Ppi_CM*beta_CM) 
    Epi_min=gamma_CM*(Epi_CM-Ppi_CM*beta_CM) 
    theta_p=Tp/mp 

    mask1=(Tp<Tpth)
    mask2=(Tp>=Tpth) & (Tp<1.0e9)
    mask3=(Tp>=1.0e9) & (Tp<5.0e9)
    mask4=(Tp>=5.0e9)

    Amax=np.zeros_like(Tp)

    Amax[mask1]=0.0 
    Amax[mask2]=5.9*func_sigma_pi(Tp[mask2])/Epi_max[mask2] 
    Amax[mask3]=9.53*pow(theta_p[mask3],-0.52)*np.exp(0.054*pow(np.log(theta_p[mask3]),2))*func_sigma_pi(Tp[mask3])/mp 
    Amax[mask4]=9.13*pow(theta_p[mask4],-0.35)*np.exp(9.7e-3*pow(np.log(theta_p[mask4]),2))*func_sigma_pi(Tp[mask4])/mp 
    
    return Amax # mb/e


# Complementary function for the differential cross-section.
def func_alpha(Tp):
# Table 5, Eq. 14, and Eq. 15,
# Tp (eV) -> kinetic energy of CR proton.

    mask1=(Tp>=Tpth) & (Tp<=20.0e9)
    mask2=(Tp>20.0e9)

    alpha=np.zeros_like(Tp)

    alpha[mask1]=1.0
    alpha[mask2]=0.5

    return alpha


# Complementary function for the differential cross-section.
def func_beta(Tp):
# Table 5, Eq. 14, and Eq. 15,
# Tp (eV) -> kinetic energy of CR proton.

    q=(Tp-1.0e9)/mp 
    mu=1.25*pow(q,1.25)*np.exp(-1.25*q) 
    theta_p=Tp/mp 
    
    mask1=(Tp>=Tpth) & (Tp<=1.0e9)
    mask2=(Tp>1.0e9) & (Tp<=4.0e9)
    mask3=(Tp>4.0e9) & (Tp<=20.0e9)
    mask4=(Tp>20.0e9) & (Tp<=100.0e9)
    mask5=(Tp>100.0e9)

    beta=np.zeros_like(Tp)

    beta[mask1]=3.29-0.2*pow(theta_p[mask1],-1.5) 
    beta[mask2]=mu[mask2]+2.45 
    beta[mask3]=1.5*mu[mask3]+4.95 
    beta[mask4]=4.2 
    beta[mask5]=4.9 
    
    return beta

# Complementary function for the differential cross-section.
def func_gamma(Tp):
# Table 5, Eq. 14, and Eq. 15,
# Tp (eV) -> kinetic energy of CR proton.
    
    q=(Tp-1.0e9)/mp 
    mu=1.25*pow(q,1.25)*np.exp(-1.25*q) 
    
    mask1=(Tp>=Tpth) & (Tp<1.0e9)
    mask2=(Tp>=1.0e9) & (Tp<=4.0e9)
    mask3=(Tp>4.0e9) & (Tp<=20.0e9)
    mask4=(Tp>20.0e9)

    gamma=np.zeros_like(Tp)

    gamma[mask1]=0.0 
    gamma[mask2]=mu[mask2]+1.45 
    gamma[mask3]=mu[mask3]+1.5 
    gamma[mask4]=1.0 

    return gamma


# Complementary function for the differential cross-section.
def func_F(Tp, Eg):
# Eq. 11 and Table 5,    
# Tp (eV) -> kinetic energy of CR proton,
# Eg (eV) -> energy of gamma ray.

    Tp, Eg=np.meshgrid(Tp, Eg, indexing='ij')

    sqrt_s=np.sqrt(2.0*mp*(Tp+2.0*mp)) 
    gamma_CM=(Tp+2.0*mp)/sqrt_s 
    beta_CM=np.sqrt(1.0-pow(gamma_CM,-2)) 
    Epi_CM=(pow(sqrt_s,2)-4.0*pow(mp,2)+pow(mpi,2))/(2.0*sqrt_s) 
    Ppi_CM=np.sqrt(pow(Epi_CM,2)-pow(mpi,2)) 
    
    Epi_max_LAB=gamma_CM*(Epi_CM+Ppi_CM*beta_CM) 
    gammapi_LAB=Epi_max_LAB/mpi 
    betapi_LAB=np.sqrt(1.0-pow(gammapi_LAB,-2)) 
    Eg_max=mpi*gammapi_LAB*(1.0+betapi_LAB)/2.0 
    Eg_min=mpi*gammapi_LAB*(1.0-betapi_LAB)/2.0 
    
    mask=(Eg>=Eg_min) & (Eg<=Eg_max)
        
    Yg=Eg+pow(mpi,2)/(4.0*Eg) 
    Yg_max=Eg_max+pow(mpi,2)/(4.0*Eg_max) # Yg_max=mpi*gammapi_LAB
    Xg=(Yg-mpi)/(Yg_max-mpi) 
    C=3.0*mpi/Yg_max 

    F=np.where(mask, pow(1.0-pow(Xg,func_alpha(Tp)),func_beta(Tp))/pow(1.0+Xg/C,func_gamma(Tp)), 0)
    
    return F


# Complementary function for the nuclear enhancement factor.
def func_GTp(Tp):
# Eq. 19,
# Tp (eV) -> kinetic energy of CR proton.

    GTp=1.0+np.log(np.maximum(1.0,func_sigma_inel(Tp)/func_sigma_inel(1.0e12*Tp**0))) 
    
    return GTp


# Nuclear enhancement factor. 
def func_enhancement(Tp):
# Eq. 24,
# Tp (eV) -> kinetic energy of CR proton.

    eps_nucl=np.zeros_like(Tp)

    mask1=(Tp>=Tpth) & (Tp<1.0e9)
    mask2=(Tp>=1.0e9)

    eps_nucl[mask1]=1.7
    eps_nucl[mask2]=1.37+0.39*10.0*np.pi*func_GTp(Tp[mask2])/func_sigma_inel(Tp[mask2]) 

    return eps_nucl


# Differential cross-section of gamma-ray from pi0 decay.
def func_d_sigma_g(Tp, Eg):
# Eq. 8,  
# Tp (eV) -> kinetic energy of CR proton,
# Eg (eV) -> energy of gamma ray.

    Amax=np.zeros_like(Tp)

    mask=(Tp>=Tpth)

    Amax[mask]=func_Amax(Tp[mask])
    d_sigma_g=Amax[:,np.newaxis]*func_F(Tp,Eg)*1.0e-27

    return d_sigma_g # cm^2/eV


###################################################################################################
# Gamma absorption
###################################################################################################


# Cross section for gamma gamma interaction
def func_sigma_gg(Eg, Ebg):

    Eg, Ebg=np.meshgrid(Eg, Ebg, indexing='ij')

    s0=Eg*Ebg/(me*me) 
    sigma_gg=np.zeros_like(Eg) 

    mask=(s0>=1.0)
    
    sigma_gg[mask]=(s0[mask]+0.5*np.log(s0[mask])-(1.0/6.0)+(1.0/(2.0*s0[mask])))*np.log(np.sqrt(s0[mask])+np.sqrt(s0[mask]-1.0)) 
    sigma_gg[mask]+=-(s0[mask]+(4.0/9.0)-(1.0/(9.0*s0[mask])))*np.sqrt(1.0-(1.0/s0[mask])) 
    sigma_gg[mask]*=1.5*sigmaT/(s0[mask]*s0[mask]) 
    
    return sigma_gg # cm^2


# Photon density of the background photons from IR thermal dust emissions 
def func_fEtd(Urad, Trad, sigma, Ebg):

    zeta=sp.special.zeta(4.0+sigma) 

    fEtd=Urad/(pow(Trad,2)*sp.special.gamma(4.0+sigma)*zeta*(np.exp(Ebg/Trad)-1.0)) 
    fEtd*=pow(Ebg/Trad,2.0+sigma)

    return fEtd # eV^-1 cm^-3


