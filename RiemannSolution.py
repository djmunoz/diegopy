# Exact Riemann solver.  Based on Toro's FORTRAN code
# Joshua Suresh, Jan 2011.
from numpy import *


tolpre = 1.0*10**-8
nriter = 30

def sample_riemann_solution(rho_L,vel_L,press_L,rho_R,vel_R,press_R,gamma,time,positions,disc):
	
	timeout = time #output time
	dL = rho_L   #initial density on left state
	uL = vel_L   #initial velocity on left state
	pL = press_L #initial pressure on left state

	dR = rho_R   #initial density on right state
	uR = vel_R   #initial velocity on right state
	pR = press_R  #initial pressure on right state

	mpa = 1.0 #pressure normalizing constant
	
	diaph = disc #Initial Discontinuity Positio

	### Compute sound speeds ###
	cL = sqrt(gamma*pL/dL)
	cR = sqrt(gamma*pR/dR)

        ### The pressure positivity condition is tested for ###
	if ((2.0/(gamma - 1.0))*(cL+cR) < (uR-uL)):
	  ### The initial data is such that the vacuum is generated ###
	  ### Program stopped ###
	  print "The initial data is such that the vacuum is generated.  Program stopped."
	  # KILL PROGRAM HERE

	
        ### Exact solution for pressure and velocity in star region is found ###
	# CALL STARPU(PM, UM, MPA)
	(pm, um) = starpu(mpa,dL,uL,pL,dR,uR,pR,gamma)


	v_vec = []
	d_vec = []
	p_vec = []

	for i in arange(0,positions.shape[0]):
		s = (positions[i]-diaph)/timeout

                ### Solution at point (x,t) = (xpos-diaph, timeout) is found ###
		(ds,us,ps) = sample(pm,um,s,dL,uL,pL,dR,uR,pR,gamma)
		d_vec = append(d_vec,ds)
		p_vec = append(p_vec,ps)
		v_vec = append(v_vec,us)

	return d_vec,v_vec,p_vec

##############################################################################################
### starpu ###
### starpu's purpose: to compute the solution for pressure and velocity in the Star Region ###
def starpu(mpa,dL,uL,pL,dR,uR,pR,gamma):

	### Compute sound speeds ###
	cL = sqrt(gamma*pL/dL)
	cR = sqrt(gamma*pR/dR)
	
	pstart = guessp(dL,uL,pL,dR,uR,pR,gamma)

	pold = pstart
	udiff = uR - uL

	print "-----------------------------------"
	print "   Iteration number   Change   "
	print "-----------------------------------"

	for i in arange(nriter):
		(fL, fLd) = prefun(pold,dL,pL,cL,gamma)
		(fR, fRd) = prefun(pold,dR,pR,cR,gamma)

		p = pold - (fL + fR + udiff)/(fLd + fRd)
		change = 2.0 * abs((p - pold)/(p + pold))
		print i+1, change

		if change < tolpre:
			print "Divergence in Newton-Raphson iteration"
			break

		if p < 0.0: p = tolpre
		pold = p


	### Compute velocity in Star Region ###
	u = 0.5*(uL + uR + fR - fL)

	print "-----------------------------------"
	print "   Pressure   Velocity   "
	print "-----------------------------------"
	print p/mpa, u
	print "-----------------------------------"


	return (p,u) #return star values


##############################################################################################
### guessp ###
### guessp's purpose: to provide a guess value for pressure pm in the Star Region.  The    ###
### choice is made according to adaptive Riemann solver using the PVRS, TRRS, and TSRS     ###
### approximate Riemann solvers.  See Sect. 9.5 of Chapt. 9 of Ref. 1 in Toro for more     ###
def guessp(dL,uL,pL,dR,uR,pR,gamma):

        ### Compute sound speeds ###
	cL = sqrt(gamma*pL/dL)
	cR = sqrt(gamma*pR/dR)
	
	quser = 2.0
	
	cup = 0.25*(dL + dR)*(cL + cR)
	ppv = 0.5*(pL + pR) + 0.5*(uL - uR)*cup
	ppv = max(0.0, ppv)
	pmin = min(pL, pR)
	pmax = max(pL, pR)
	qmax = pmax/pmin

        ### Compute gamma-related constants ###
	g1,g2, g3, g4, g5, g6, g7, g8 = get_gamma_constants(gamma)
	
	
	if qmax<quser and pmin<ppv and ppv<pmax:
		### Select PVRS Riemann solver ###
		pm = ppv
	elif ppv<pmin:
		### Select Two-Rarefaction Riemann Solver ###
		pq = (pL/pR)**g1
		um = (pq*uL/cL + uR/cR + g4*(pq - 1.0))/(pq/cL + 1.0/cR)
		ptL = 1.0 + g7*(uL - um)/cL
		ptR = 1.0 + g7*(um - uR)/cR
		pm = 0.5*(pL*ptL**g3 + pR*ptR**g3)
	else :
		### Select Two-Shock Riemann Solver with PVRS as estimate ###
		geL = sqrt((g5/dL)/(g6*pL + ppv))
		geR = sqrt((g5/dR)/(g6*pR + ppv))
		pm = (geL*pL + geR*pR - (uR - uL))/(geL + geR)
		
	return pm

##############################################################################################
### prefun ###
### prefun's purpose: to evaluate the pressure functions fL and fR in exact Riemann solver ###
def prefun(p,dk,pk,ck,gamma):

	### Compute gamma-related constants ###
	g1,g2, g3, g4, g5, g6, g7, g8 = get_gamma_constants(gamma)
	

	
	if p <= pk:
		### Rarefaction Wave ###
		prat = p/pk
		f = g4*ck*(prat**g1 - 1.0)
		fd = (1.0/(dk*ck))*prat**(-g2)

	else :
		### Shock Wave ###
		ak = g5/dk
		bk = g6*pk
		qrt = sqrt(ak/(bk + p))
		f = (p - pk)*qrt
		fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt

	return (f,fd)

##############################################################################################
### sample ###
### sample's purpose: to sample the solution throughout the wave pattern.  Pressure pm and ###
### velocity um in the Star Region are known.  Sampling is performed in terms of the       ###
### "speed" S = x/t.  Sampled values are d,u,p
def sample(pm,um,s,dL,uL,pL,dR,uR,pR,gamma):

	### Compute sound speeds ###
	cL = sqrt(gamma*pL/dL)
	cR = sqrt(gamma*pR/dR)
	### Compute gamma-related constants ###
	g1,g2, g3, g4, g5, g6, g7, g8 = get_gamma_constants(gamma)

	
	if s <= um:
		### Sampling point lies to the left of the contact discontinuity ###
		if pm <= pL:
			### Left Rarefaction ###
			shL = uL - cL

			if s <= shL:
				### Sampled point is left data state ###
				d = dL
				u = uL
				p = pL

			else :
				cmL = cL*(pm/pL)**g1
				stL = um - cmL

				if s > stL:
					### Sampled point is Star Left state ###
					d = dL*(pm/pL)**(1.0/gamma)
					u = um
					p = pm

				else :
					### Sampled point is inside left fan ###
					u = g5*(cL + g7*uL + s)
					c = g5*(cL + g7*(uL - s))
					d = dL * (c/cL)**g4
					p = pL * (c/cL)**g3


		else :
			### Left shock ###
			pmL = pm/pL
			sL = uL - cL*sqrt(g2*pmL + g1)

			if s <= sL:
				### Sampled point is left data state ###
				d = dL
				u = uL
				p = pL

			else :
				### Sampled point is Star Left state ###

				d = dL*(pmL + g6)/(pmL*g6 + 1.0)
				u = um
				p = pm


	else :
		### Sampling point lies to the right of the contact discontinuity ###
		if pm > pR:
			### Right Shock ###

			pmR = pm/pR
			sR = uR + cR*sqrt(g2*pmR + g1)

			if s >= sR:
				### Sampled point is right data state ###
				d = dR
				u = uR
				p = pR

			else :
				### Sampled point is Star Right state ###
				d = dR*(pmR + g6)/(pmR*g6 + 1.0)
				u = um
				p = pm

		else :
			### Right Rarefaction ###
			shR = uR + cR

			if s >= shR:
				### Sampled point is right data state ###
				d = dR
				u = uR
				p = pR

			else :
				cmR = cR*(pm/pR)**g1
				stR = um + cmR

				if s <= stR:
					### Sampled point is Star Right state ###
					d = dR*(pm/pR)**(1.0/gamma)
					u = um
					p = pm

				else :
					### Sampled point is inside left fan ###
					u = g5*(-cR + g7*uR + s)
					c = g5*(cR - g7*(uR - s))
					d = dR*(c/cR)**g4
					p = pR*(c/cR)**g3

	return (d,u,p)

##############################################################################################
def get_gamma_constants(gamma):
	g1 = (gamma - 1.0)/(2.0*gamma)
	g2 = (gamma + 1.0)/(2.0*gamma)
	g3 = 2.0*gamma/(gamma - 1.0)
	g4 = 2.0/(gamma - 1.0)
	g5 = 2.0/(gamma + 1.0)
	g6 = (gamma - 1.0)/(gamma + 1.0)
	g7 = (gamma - 1.0)/2.0
	g8 = gamma - 1.0

	return 	g1,g2, g3, g4, g5, g6, g7, g8


##############################################################################################
def sample_riemann_isothermal_solution(rho_L,vel_L,rho_R,vel_R,csnd,time,positions,disc):
	
	timeout = time #output time
	dL = rho_L   #initial density on left state
	uL = vel_L   #initial velocity on left state
	pL = rho_L * csnd**2 #initial pressure on left state

	dR = rho_R   #initial density on right state
	uR = vel_R   #initial velocity on right state
	pR = rho_R * csnd**2  #initial pressure on right state

	mpa = 1.0 #pressure normalizing constant
	
	diaph = disc #Initial Discontinuity Positio


	### Exact solution for pressure and velocity in star region is found ###
	# CALL STARPU(PM, UM, MPA)
	(dm, um) = stardu(dL,uL,dR,uR,csnd)


	v_vec = []
	d_vec = []
	p_vec = []

	for i in arange(0,positions.shape[0]):
		s = (positions[i]-diaph)/timeout

                ### Solution at point (x,t) = (xpos-diaph, timeout) is found ###
		(ds,us) = sample_isothermal(dm,um,s,dL,uL,dR,uR,csnd)
		d_vec = append(d_vec,ds)
		p_vec = append(p_vec,ds * csnd**2)
		v_vec = append(v_vec,us)

	return d_vec,v_vec,p_vec


##############################################################################################
##############################################################################################
### stardu ###
### stardu's purpose: to compute the solution for density and velocity in the Star Region for
### strictly isothermal flow
def stardu(dL,uL,dR,uR,csnd):

	udiff = (uR - uL)/ csnd;
	
	if(udiff > 0):
		rho = sqrt(dL * dR * exp(-udiff))
	else:
		rho = 0.5 * (dL + dR)


	if( (dL <= 0) | (dR <= 0)):
		print "isothermal Riemann solver was called with zero or negative density\n"
		return
	
	udiff = uR - uL

	print "-----------------------------------"
	print "   Iteration number   Change   "
	print "-----------------------------------"

	for i in arange(nriter):
		(fL, fLd) = isothermal_function(rho,dL)
		(fR, fRd) = isothermal_function(rho,dR)

		rhoold = rho
		drho = -0.5 * (fL+ fR + udiff) / (fLd + fRd)
		if (abs(drho) > 0.25 * rho):
			drho = 0.25 * rho * abs(drho) / drho

		rho+=drho

		change = 2.0 * abs((rho - rhoold)/(rho + rhoold))
		print i+1, change

		if change < tolpre:
			print "Divergence in Newton-Raphson iteration"
			break

		if rho < 0.0: rho = tolpre



	### Compute velocity in Star Region ###
	u = 0.5*(uL + uR + csnd *(fR - fL))

	print "-----------------------------------"
	print "   Density   Velocity   "
	print "-----------------------------------"
	print rho, u
	print "-----------------------------------"


	return (rho,u) #return star values


##############################################################################################
### sample_isothermal ###
### sample's purpose: to sample the solution throughout the wave pattern.
###  Sampling is performed in terms of the       ###
### "speed" S = x/t.  Sampled values are d,u
def sample_isothermal(dm,um,s,dL,uL,dR,uR,csnd):

	if s <= um:
		### Sampling point lies to the left of the contact discontinuity ###

		if (dm <= dL):  ## left fan
			shL = uL - csnd

			if s <= shL:
				### Sampled point is left data state ###
				d = dL
				u = uL
			else :
				stL = um - csnd

				if s > stL:
					### Sampled point is Star Left state ###
					d = dm
					u = um
				else :
					### Sampled point is inside left fan ###
					u = um
					d = dL * exp(-(s + csnd - uL)/csnd)

		else :
			### Left shock ###
			sL = (dL * uL - dm * um)/(dL - dm)

			if s <= sL:
				### Sampled point is left data state ###
				d = dL
				u = uL
			else :
				### Sampled point is Star Left state ###

				d = dm
				u = um

	else :
		### Sampling point lies to the right of the contact discontinuity ###
		if dm > dR:
			### Right Shock ###

			sR = (dR * uR - dm * um)/(dR - dm)

			if s >= sR:
				### Sampled point is right data state ###
				d = dR
				u = uR
			else :
				### Sampled point is Star Right state ###
				d = dm
				u = um

		else :
			### Right Rarefaction ###
			shR = uR + csn

			if s >= shR:
				### Sampled point is right data state ###
				d = dR
				u = uR

			else :
				stR = um + csnd

				if s <= stR:
					### Sampled point is Star Right state ###
					d = dm
					u = um

				else :
					### Sampled point is inside left fan ###
					u = s - csnd
					d = dR * exp(( s- csnd - uR)/ csnd)

	return (d,u)

##############################################################################################
def isothermal_function(d,dk):

	if (d <= dk):
	        ### Rarefaction Wave ###	
		f = log(d/dk)
		fd = 1.0/dk
	else:
           	### Shock Wave ###
		f = (d - dk) / sqrt(d * dk)
		fd = 0.5 / d * (sqrt(d/dk) + sqrt(dk/d))

	return (f,fd)
