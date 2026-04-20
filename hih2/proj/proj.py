import numpy as np
import astropy.units as u
import astropy.constants as c

# #########
#  * * * * 
# #########

def p24(nh, met, uv, scale, dens_unit = u.cm**-2, scale_unit = u.pc):
	"""
	Return projected molecular fraction from Polzin et al. (2024)

	Parameters: 
		nh (arr-like): neutral hydrogen column density (cm^-2 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): spatial scale (pc by default)
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary
	scale *= scale_unit.to(u.pc)

	w = 0.27 - 1e-2 * (9.25*np.log10(met)**2 + 9.64*np.log10(met))
	ynorm = 21.96 - 0.19*np.log10(scale)
	y = ynorm * np.exp(-0.5*((np.log10(met) + 1.5)/6.84)**2)

	corr = (1 + 0.13*np.log10(0.1/met)*np.log10(scale/10))

	ntr = uv**w * 10**y * corr

	x = np.log(nh/ntr)

	r0 = 3.5e-17
	disk_scale = 150*u.pc.to(u.cm)
	q = 3*u.Myr.to(u.s) * r0 * (met/0.1)**1.3 * nh/disk_scale
	maxf = 1 - np.exp(-q)
	fm = 1/(1 + 2*(1 - maxf)/maxf)

	fmol = fm/(1 + np.exp(-(1 + 1.35*(1e-2/met)**0.25*(scale/10)**0.6 + 3.4*(met/0.6)**(0.02))*x + np.log(fm/1.65e-4)))
	
	return fmol


def gk11(nh, met, uv, scale = None, dens_unit = u.cm**-2, scale_unit = u.pc):
	"""
	Return molecular fraction from Gnedin & Kravtsov (2011)

	Parameters: 
		nh (arr-like): neutral hydrogen column density (cm^-2 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): spatial scale (pc by default; not required)
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary
	dmw = met

	nh_ = ((nh*u.cm**-2).to(u.pc**-2)*c.m_p.to(u.M_sun)).value #surface density, Msun/pc^2

	nstar = 25 #cm^-3
	dstar = 1.5e-3*np.log(1 + (3*uv)**(1.7))
	alpha = 5*(uv/2)/(1 + (uv/2)**2)
	s = 0.04/(dstar + dmw)
	g = (1 + alpha*s + s**2)/(1 + s)
	lam = np.log(1 + g*dmw**(3/7)*(uv/15)**(4/7))
	sig_hi = 40*(lam**(4/7)/dmw)/np.sqrt(1 + uv*dmw**2) #Eq 10 + 14

	fmol = 1 - sig_hi/nh_
	fmol[sig_hi > nh_] = 0

	return fmol


def gd14(nh, met, uv, scale, dens_unit = u.cm**-2, scale_unit = u.pc):
	"""
	Return molecular fraction from Gnedin & Draine (2014)

	Parameters: 
		nh (arr-like): neutral hydrogen number density (cm^-2 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): scale of grid cell in (pc by default)
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary
	scale *= scale_unit.to(u.pc)

	s = scale/100
	dmw = met

	nh_ = ((nh*u.cm**-2).to(u.pc**-2)*c.m_p.to(u.M_sun)).value #surface density, Msun/pc^2

	dstar = 0.17*(2 + s**5)/(1 + s**5)
	ustar = 9*dstar/s
	nstar = 14 * dstar**(1/2)/s #cm^-3
	g = (dmw**2 + dstar**2)**(1/2)
	a = 0.5 + 1/(1 + np.sqrt(uv*dmw**2/600))
	sig_r = (50/g)*np.sqrt(0.01 + uv)/(1 + 0.69*np.sqrt(0.01 + uv))

	rmol = (nh_/sig_r)**a


	fmol = rmol/(1 + rmol)

	return fmol


def kmt09b(nh, met, uv = None, scale = None, dens_unit = u.cm**-2, scale_unit = u.pc):
	"""
	Return molecular fraction from Krumholz et al. (2009b)

	Parameters: 
		nh (arr-like): neutral hydrogen column density (cm^-2)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength (not required)
		scale (arr-like): scale of grid cell (pc by default; not required)
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary

	sig_comp = nh*c.m_p.to(u.g).value

	chi = 0.77*(1 + 3.1*met**0.365)
	sig0 = sig_comp/(1*(u.Msun/u.pc**2).to(u.g/u.cm**2))
	s = np.log(1 + 0.6*chi)/(0.04*sig0 * met)
	delta = 0.0712*(0.1/s + 0.675)**(-2.8)

	fmol = 1 - (1 + ((3/4) * (s/(1 + delta)))**(-5))**(-1/5)

	return fmol

def k13(nh, met, uv, scale = None, rho_sd = 1e-2, fc = 1, iter_ = True, niter = 10, 
		dens_unit = u.cm**-2, scale_unit = u.pc, sddens_unit = u.Msun/u.pc**3):
	"""
	Return molecular fraction from Krumholz (2013)

	Parameters: 
		nh (arr-like): neutral hydrogen columns density (cm^-2 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): scale of grid cell (pc by default; not required)
		rho_sd (arr-like): density of stars and dark matter (Msun/pc^3 by default)
		fc (int): ~Unity on scales <~100 pc. Scales ~1 kpc -> fc = 5.
		iter (bool): If iter, iterate niter times to compute fmol.
		niter (int): Number of iterations since fmol necessary to compute fmol.
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
		sddens_unit (astropy units object): units attached to input rho_sd
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary
	rho_sd *= sddens_unit.to(u.Msun/u.pc**3)
	
	a = 5
	cw = 8*u.km/u.s
	fw = 0.5
	zeta_d = 0.33
	#can set single value between 1e-5 and 1e-1 Msun/pc**3 -- varies with galaxy structure though
	rho_sd = rho_sd*u.Msun/u.pc**3 #rho_sd
	G = c.G.to(u.pc/u.Msun*(u.km/u.s)**2)
	## since generally low fH2 regime will take Sigma_HI ~ Sigma_H (see K13, Sec. 2.2)
	sig0 = (nh*c.m_p.to(u.g).value)*(u.g/u.cm**2).to(u.Msun/u.pc**2)*u.Msun/u.pc**2
	ncnm_2p = 23*uv*(4.1/(1 + 3.1*met**0.365)) #cm-3
	if not iter_:
		ncnm_hydro = (np.sqrt(2*np.pi*G*zeta_d*fw*rho_sd/a) * (cw/(1.1*c.k_B*243*u.K))*sig0).to(u.cm**-3).value

		ncnm = np.maximum(ncnm_2p, ncnm_hydro)
		chi = 7.2*uv*(ncnm/10)**(-1)

		tau_c = 0.066*fc*met*sig0.value
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

			tau_c = 0.066*fc*met*sig0.value
			s = np.log(1 + 0.6*chi + 0.01*chi**2)/(0.6*tau_c)

			new_fmol = 1 - (3/4)*s/(1 + 0.25 * s)
			new_fmol[s >= 2] = 0

			fmol = 0.3*fmol + 0.7*new_fmol

	return fmol


def s14(nh, met, uv, scale = None, fc = 1, dens_unit = u.cm**-2, scale_unit = u.pc):
	"""
	Return molecular fraction from Sternberg et al. (2014)

	Following Diemer et al. (2018), Section C.5

	Parameters: 
		nh (arr-like): neutral hydrogen column density (cm^-2 by default)
		met (arr-like): metallicity in solar units
		uv (arr-like): normalized UV field strength
		scale (arr-like): scale of grid cell (pc by default; not required)
		fc (int): clumping factor
		dens_unit (astropy units object): units attached to input nh
		scale_unit (astropy units object): units attached to input scale
	"""
	nh *= dens_unit.to(u.cm**-2) # convert units if necessary
	nh_ = nh/(500*u.pc.to(u.cm)) #approximate height of the disk

	g = 3e-5 * met * (9.9/(1 + 8.9*met))
	d0 = 5.8e-11*uv #s-1
	r = 3e-17 #cm-3 s-1
	ag = (d0*g)/(r*nh_)
	nhi = 5.3e20 * ((1/met) * np.log((ag/fc)/2 + 1)) #cm-2

	fmol = 1 - nhi/nh
	fmol[nhi > nh] = 0

	return fmol