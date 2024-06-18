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

# The intensity of CR electrons fitted with Voyager and AMS data
def func_sptr_p_Voyager(E): # Voyager fit

    C=1.882e-9 
    alpha=0.129056 
    beta=2.82891 
    E0=624.5e6 
    
    f0=C*pow(E/1.0e6,alpha)*pow(1.0+(E/E0),-beta) 
    
    return f0 # eV^-1 cm^-2 s^-1 sr^-1

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

# Auxiliary functions for proton ionization differential cross-section and cross-section
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

# Electron ionization differential cross-section and cross-section
def func_df_dw(w):

	c_df, d_df, e_df, f_df 
	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 

	df_dw=(c_df/pow(1.0+w,3))+(d_df/pow(1.0+w,4))+(e_df/pow(1.0+w,5))+(f_df/pow(1.0+w,6)) 
	return df_dw 

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

def func_d_sigma_e_E(Ee_p, Ee_s):# Differential cross-section of electron (RBEB binary-encounter-Bethe)

	c_df, d_df, e_df, f_df, N_i 
	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 
	N_i=(c_df/2.0)+(d_df/2.0)+(e_df/2.0)+(f_df/2.0) 	

	U_H2, B_H2,  u_primed, b_primed, beta_u, beta_b 	
	U_H2=39.603 # Average kinetic energy of orbital electron
	B_H2=IH2 # Binding energy
	u_primed=U_H2/me 
	b_primed=B_H2/me 
	beta_u=np.sqrt(1.0-(1.0/pow(1.0+u_primed,2))) 
	beta_b=np.sqrt(1.0-(1.0/pow(1.0+b_primed,2))) 

	d_sigma_E, t, t_primed, w, beta_t 
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

def func_sigma_e(Ee_p):

	c_df, d_df, e_df, f_df, N_i 
	c_df=1.1262 
	d_df=6.382 
	e_df=-7.8055 
	f_df=2.144 
	N_i=(c_df/2.0)+(d_df/2.0)+(e_df/2.0)+(f_df/2.0) 

	U_H2, B_H2,  u_primed, b_primed, beta_u, beta_b 	
	U_H2=39.603 # Average kinetic energy of orbital electron
	B_H2=IH2 # Binding energy
	u_primed=U_H2/me 
	b_primed=B_H2/me 
	beta_u=np.sqrt(1.0-(1.0/pow(1.0+u_primed,2))) 
	beta_b=np.sqrt(1.0-(1.0/pow(1.0+b_primed,2))) 	

	sigma_e, t, t_primed, beta_t 
	t=Ee_p/IH2 
	t_primed=Ee_p/me 
	beta_t=np.sqrt(1.0-(1.0/pow(1.0+t_primed,2))) 	

	sigma_e=0.0 

	sigma_e+=func_D(t)*(np.log(pow(beta_t,2)/(1.0-pow(beta_t,2)))-pow(beta_t,2)-np.log(2.0*b_primed)) 
	sigma_e+=(2.0-(N_i/N))*(1.0-(1.0/t)-(np.log(t)/(t+1))*((1.0+2.0*t_primed)/pow(1.0+(t_primed/2.0),2))) 
	sigma_e+=(2.0-(N_i/N))*((pow(b_primed,2)/pow(1.0+(t_primed/2.0),2))*((t-1.0)/2.0)) 
	sigma_e*=((Se/(2.0*b_primed))/(pow(beta_t,2)+pow(beta_u,2)+pow(beta_b,2))) 
	
	return sigma_e 


# Function for interpolation
def log_interp1d(xx, yy, kind='linear'):

    logx=np.log10(xx)
    logy=np.log10(yy)
    lin_interp=sp.interpolate.interp1d(logx,logy,kind=kind)
    log_interp=lambda zz: np.power(10.0,lin_interp(np.log10(zz)))

    return log_interp


# Fit of the Voyager spectrum for protons 
def func_jLocISM_p(E):
    
    C=1.882e-9 
    alpha=0.129056 
    beta=2.82891 
    E0=624.5e6 
    
    f0=12.5e-10*pow(E/1.0e6,0.35)/((1.0+pow(E/80.0e6,1.3))*pow(1.0+pow(E/2.2e9,1.9/2.1),2.1)) 

    return f0 # eV^-1 cm^-2 s^-1 sr^-1

