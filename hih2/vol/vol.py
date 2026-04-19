import numpy as np
import astropy.units as u
import astropy.constants as c

# #########
#  * * * * 
# #########

def p24(nh, met, uv, dens_unit = u.cm**-3):
	"""
	Return volumetric molecular fraction from Polzin et al. (2024a)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		dens_unit (astropy units object): units attached to input nh
	"""

	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary

	d = 0.0199*met
	a = -34.7*uv**0.32 + 2.25*met**0.3
	b = -53.9*uv**0.31
	c = met/0.2
	ntr= a*np.log10(d) + b + c
	ntr[ntr < 0.1] = 0.1

	x = np.log(nh/ntr)

	r0 = 3.5e-17
	q = 30*u.Myr.to(u.s) * r0 * 0.2*(met/0.2)**1.3 * nh
	maxf = 1 - np.exp(-q)
	fm = 1/(1 + 2*(1 - maxf)/maxf)

	fmol = fm/(1 + np.exp(-(7.6*(met)**0.25)*x + np.log(fm/6e-4)))

	return fmol


def gk11(nh, met, uv, dens_unit = u.cm**-3):
	"""
	Return molecular fraction from Gnedin & Kravtsov (2011)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		dens_unit (astropy units object): units attached to input nh
	"""
	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary

	dmw = met

	nstar = 25 #cm^-3
	dstar = 1.5e-3*np.log(1 + (3*uv)**(1.7))
	alpha = 5*(uv/2)/(1 + (uv/2)**2)
	s = 0.04/(dstar + dmw)
	g = (1 + alpha*s + s**2)/(1 + s)
	lam = np.log(1 + g*dmw**(3/7)*(u/15)**(4/7))
	x = lam**(3/7) * np.log(dmw * nh/(lam*nstar))

	fmol = 1/(1 + np.exp(-4*x - 3*x**3))
	
	return fmol


def gd14(nh, met, scale, uv = None, dens_unit = u.cm**-3, scale_unit = u.cm):
	"""
	Return molecular fraction from Gnedin & Draine (2014)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3 by default)
		met (arr-like): metallicity in solar units
		scale (arr-like): scale of grid cell in cm
		uv (arr-like): normalized UV field strength
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary
	scale *= scale_unit.to(u.cm).value


	s = scale*u.cm.to(u.pc)/100
	dmw = met

	dstar = 0.17*(2 + s**5)/(1 + s**5)
	ustar = 9*dstar/s
	nstar = 14 * dstar**(1/2)/s #cm^-3
	g = (dmw**2 + dstar**2)**(1/2)
	lam = np.log(1 + (0.05/g + u)**(2/3) * g**(1/3)/ustar)
	nhalf = nstar*lam/g
	w = 0.8 + lam**(1/2)/s**(1/3)
	x = w*np.log(nh/nhalf)

	fmol = 1/(1 + np.exp(-x*(1 - 0.02*x + 0.001*x**2)))

	return fmol

def kmt09b(nh, z, l, uv = None, dens_unit = u.cm**-3, scale_unit = u.cm):
	"""
	Return molecular fraction from Krumholz et al. (2009b)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3)
		z (arr-like): metallicity in solar units
		u (arr-like): normalized UV field strength
		l (arr-like): scale of grid cell in cm
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary
	scale *= scale_unit.to(u.cm).value

	sig_comp = nh*l*c.m_p.to(uu.g).value #should be on 100 pc scales, but will see if this works

	chi = 0.77*(1 + 3.1*z**0.365)
	sig0 = sig_comp/(1*(uu.Msun/uu.pc**2).to(uu.g/uu.cm**2))
	s = np.log(1 + 0.6*chi)/(0.04*sig0 * z)
	delta = 0.0712*(0.1/s + 0.675)**(-2.8)

	fmol = 1 - (1 + ((3/4) * (s/(1 + delta)))**(-5))**(-1/5)

	return fmol

def k13(nh, z, uv, l, rho_sd = 1e-2, fc = 1, iter_ = True, niter = 10, dens_unit = u.cm**-3, scale_unit = u.cm):
	"""
	Return molecular fraction from Krumholz (2013)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3)
		z (arr-like): metallicity in solar units
		u (arr-like): normalized UV field strength
		l (arr-like): scale of grid cell in cm
		fc (int) : ~Unity on scales <~100 pc. Scales ~1 kpc -> fc = 5.
		iter (bool): If iter, iterate niter times to compute fmol.
		niter (int): Number of iterations since fmol necessary to compute fmol.
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary
	scale *= scale_unit.to(u.cm).value

	a = 5
	cw = 8*u.km/u.s
	fw = 0.5
	zeta_d = 0.33
	if type(rho_sd) == np.ndarray:
		## accounting for input from simulation directly
		rho_sd = rho_sd*u.Msun/(u.pc**3)
		# print(rho_sd)
	if type(rho_sd) != np.ndarray:
		#set single value between 1e-5 and 1e-1 Msun/pc**3 -- varies with galaxy structure, will use constant for now
		rho_sd = rho_sd#*uu.Msun/uu.pc**3 #rho_sd
	G = c.G.to(u.pc/u.Msun*(u.km/u.s)**2)
	## since generally low fH2 regime will take Sigma_HI ~ Sigma_H (see K13, Sec. 2.2)
	sig0 = (nh*l*c.m_p.to(u.g).value)*(u.g/u.cm**2).to(u.Msun/u.pc**2)*u.Msun/u.pc**2
	ncnm_2p = 23*u*(4.1/(1 + 3.1*z**0.365)) #cm-3
	if not iter_:
		ncnm_hydro = (np.sqrt(2*np.pi*G*zeta_d*fw*rho_sd/a) * (cw/(1.1*c.k_B*243*u.K))*sig0).to(u.cm**-3).value

		ncnm = np.maximum(ncnm_2p, ncnm_hydro)
		chi = 7.2*uv*(ncnm/10)**(-1)

		tau_c = 0.066*fc*z*sig0.value
		s = np.log(1 + 0.6*chi + 0.01*chi**2)/(0.6*tau_c)

		fmol = 1 - (3/4)*s/(1 + 0.25 * s)
		fmol[s >= 2] = 0

	if iter_:
		fmol = np.full_like(nh, 0.5)
		for i in range(niter):
			## for right now, not worrying about an error threshold -- ~following Diemer+18
			rh2 = fmol/(1 - fmol)
			sig_hi = (1 - fmol)*sig0 #assuming all neutral hydrogen
			pth = ((np.pi*G*sig_hi**2/(4*a))*
				  (1 + 2*rh2 + np.sqrt((1 + 2*rh2)**2 + 
						((32*zeta_d*a*fw*cw**2*rho_sd)/(np.pi*G*sig_hi**2)))))
			ncnm_hydro = (pth/(1.1*c.k_B*243*u.K)).to(u.cm**-3).value

			ncnm = np.maximum(ncnm_2p, ncnm_hydro)
			chi = 7.2*uv*(ncnm/10)**(-1)

			tau_c = 0.066*fc*z*sig0.value
			s = np.log(1 + 0.6*chi + 0.01*chi**2)/(0.6*tau_c)

			new_fmol = 1 - (3/4)*s/(1 + 0.25 * s)
			new_fmol[s >= 2] = 0

			fmol = 0.3*fmol + 0.7*new_fmol

	return fmol

def s14(nh, met, uv, scale, fc = 1, dens_unit = u.cm**-3, scale_unit = u.cm):
	"""
	Return molecular fraction from Sternberg et al. (2014)

	Following Diemer et al. (2018), Section C.5

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-3)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): scale of grid cell in cm
		fc (int): clumping factor
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-3).value # convert units if necessary
	scale *= scale_unit.to(u.cm).value

	g = 3e-5 * met * (9.9/(1 + 8.9*met))
	d0 = 5.8e-11*uv #s-1
	r = 3e-17 #cm-3 s-1
	ag = (d0*g)/(r*nh)
	nhi = 5.3e20 * ((1/met) * np.log((ag/fc)/2 + 1)) #cm-2
	nh_ = nh*scale

	fmol = 1 - nhi/nh_
	fmol[nhi > nh_] = 0

	return fmol