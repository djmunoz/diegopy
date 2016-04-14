import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import conversions as co
import readsnap as rs

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def PlotRhoVsTemp(base, num, arepo=0, BINS=200, Omegab=0.044, xmin=-3., xmax=7., ymin=2., ymax=8., format='.eps'):

        if arepo > 0:
                filename=base+'snap_arepo_'+str(num).zfill(3)
        else:
               filename=base+'snap_gadget_'+str(num).zfill(3)
        head=rs.snapshot_header(filename)
	mass=np.float64(rs.read_block(filename, "MASS", parttype=0, arepo=arepo))
        rho=np.float64(rs.read_block(filename, "RHO ", parttype=0, arepo=arepo))
        u=np.float64(rs.read_block(filename, "U   ", parttype=0, arepo=arepo))
        Nelec=np.float64(rs.read_block(filename, "NE  ", parttype=0, arepo=arepo))

        temp=co.GetTemp(u, Nelec, 5./3.)
        rho_b=Omegab*co.GetRhoCrit()
        rho/=rho_b

        print "z = ", head.redshift
        print "min/max of T [K]            = ", min(temp), max(temp)
        print "min/max of rho/rho_b        = ", min(rho), max(rho)

        rho=np.log10(rho)
        temp=np.log10(temp)

        if (xmin==0.) & (ymin==0.) & (xmax==0.) & (ymax==0.):
                print "range not specified -> adjusting min/max"
                xmin=min(rho)
                xmax=max(rho)
                ymin=min(temp)
                ymax=max(temp)

        Z,x,y=np.histogram2d(rho,temp, range=[[xmin,xmax],[ymin,ymax]], weights=mass, bins=BINS, normed=True)
        Z=np.log10(Z)
        
	Zmin=Z[Z>-np.inf].min()
	Zmax=Z.max()

	print "min/max of log10(histogram) = ", Zmin, Zmax

        fig = plt.figure(1, figsize=(10.0,10.0))

        ax = fig.add_subplot(1,1,1)

        im=ax.imshow(Z.T, vmin=Zmin, vmax=Zmax, origin='lower',interpolation='nearest', extent=[xmin, xmax, ymin, ymax], cmap=cm.get_cmap('jet'))
        ax.contour(Z.T,origin='lower',extent=[xmin, xmax, ymin, ymax], colors='black', vmin=Zmin, vmax=Zmax)
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))
        ax.set_xlabel(r'log $\rho/\rho_{b}$', fontsize=20)
        ax.set_ylabel('log T [K]', fontsize=20)

	plt.colorbar(im, shrink=0.5)
        if arepo > 0:
                plt.suptitle('Arepo   z='+str(round(head.redshift,2)))
                plt.savefig('Rho_vs_T_arepo_'+str(num).zfill(3)+format)
        else:
                plt.suptitle('Gadget   z='+str(round(head.redshift,2)))
                plt.savefig('Rho_vs_T_gagdet_'+str(num).zfill(3)+format)

 	fig.clf() 

def star_mass_dist(basename, num, arepo=1, bins=100, format=".eps"):	
	filename=basename+str(num).zfill(3)

        head=rs.snapshot_header(filename)
	mass = rs.read_block(filename,"MASS", parttype=4,arepo=arepo)
	print mass
	meanmass=mass.mean()

	print "mean star mass       = ", meanmass
        print "min/max of star mass = ", mass.min(),mass.max()


        fig = plt.figure(1, figsize=(10.0,10.0))
        ax = fig.add_subplot(1,1,1)

	# the histogram of the data
	n, bins, patches = plt.hist(np.log10(mass/meanmass), bins, normed=1, facecolor='blue')

	ax.set_xlabel('log[m/<m>]')
	ax.set_ylabel('df / dlog[m/<m>]')
	plt.suptitle('z='+str(round(head.redshift,2)))
	plt.savefig("star_mass_dist_"+str(num).zfill(3)+format)
        fig.clf()

def HaloProfiles(basename, num, centre, r200, rmin=0.05, rmax=10.0, bins=50, arepo=1, gamma=5./3., format=".eps"):	
	filename=basename+str(num).zfill(3)

        head=rs.snapshot_header(filename)
	mass_gas = rs.read_block(filename,"MASS", parttype=0,arepo=arepo).astype('float64')	
	mass_DM = rs.read_block(filename,"MASS", parttype=1,arepo=arepo).astype('float64')		
	pos_gas = rs.read_block(filename,"POS ", parttype=0,arepo=arepo).astype('float64')	
	pos_DM = rs.read_block(filename,"POS ", parttype=1,arepo=arepo).astype('float64')		
	u = rs.read_block(filename,"U   ", parttype=0,arepo=arepo).astype('float64')	
	rho = rs.read_block(filename,"RHO ", parttype=0,arepo=arepo).astype('float64')		
	Nele = rs.read_block(filename,"NE  ", parttype=0,arepo=arepo).astype('float64')			
	
	print "Centre     = ", centre
	print "R200       = ", r200
	print "rmin/rmax  = ", rmin, rmax 
	
	x=pos_gas[:,0] - centre[0]
	y=pos_gas[:,1] - centre[1]
	z=pos_gas[:,2] - centre[2]		
	r_gas=np.sqrt(x**2. + y**2. + z**2.) / r200
	

	x=pos_DM[:,0] - centre[0]
	y=pos_DM[:,1] - centre[1]
	z=pos_DM[:,2] - centre[2]		
	r_DM=np.sqrt(x**2. + y**2. + z**2.) / r200
	
	rmin=np.log10(rmin)
	rmax=np.log10(rmax)
	
    	dlog10=(rmax-rmin)/bins
	rho_DM_bin=np.zeros(bins)		
	rho_gas_bin=np.zeros(bins)	
	temp_bin=np.zeros(bins)
	entropy_bin=np.zeros(bins)

	rbinm=10.**((np.arange(bins)+0.5)*dlog10 + rmin)

	for n in range(0,bins):
		r1=10.**((n+0.)*dlog10 + rmin)
		r2=10.**((n+1.)*dlog10 + rmin)		
		index_gas=((r_gas>r1) & (r_gas<r2)).nonzero()
		index_DM=((r_DM>r1) & (r_DM<r2)).nonzero()		

		totmass_gas=mass_gas[index_gas].sum()
		totmass_DM=mass_DM[index_DM].sum()		

		rho_gas_bin[n]=totmass_gas/(4.*np.pi/3.*(r2**3.-r1**3.)*r200**3.)
		rho_DM_bin[n]=totmass_DM/(4.*np.pi/3.*(r2**3.-r1**3.)*r200**3.)		

		if (totmass_gas > 0.):
			entropy_bin[n]=np.average(co.GetEntropy(u[index_gas],rho[index_gas],gamma), weights=mass_gas[index_gas]) 
			temp_bin[n]=np.average(co.GetTemp(u[index_gas],Nele[index_gas],gamma), weights=mass_gas[index_gas])
								
	
        fig = plt.figure(1, figsize=(10.0,10.0))
        ax = fig.add_subplot(2,2,1)
	ax.set_xlabel('$r/r_{200}$')
	ax.set_ylabel(r'$\rho_{DM}$ [$h^2$ M$_\odot$ Kpc$^{-3}$]')
	ax.loglog()
	ax.plot(rbinm, 10.0**10.0*rho_DM_bin)
	ax.set_xlim((10.0**rmin,10.0**rmax))	

        ax = fig.add_subplot(2,2,2)
	ax.set_xlabel('$r/r_{200}$')
	ax.set_ylabel(r'$\rho_{gas}$ [$h^2$ M$_\odot$ Kpc$^{-3}$]')
	ax.loglog()
	ax.plot(rbinm, 10.0**10.0*rho_gas_bin)
	ax.set_xlim((10.0**rmin,10.0**rmax))	

        ax = fig.add_subplot(2,2,3)
	ax.set_xlabel('$r/r_{200}$')
	ax.set_ylabel('$T_{gas}$ [K]')
	ax.loglog()
	ax.plot(rbinm, temp_bin)
	ax.set_xlim((10.0**rmin,10.0**rmax))	
	
        ax = fig.add_subplot(2,2,4)
	ax.set_xlabel('$r/r_{200}$')
	ax.set_ylabel('Entropy')
	ax.loglog()
	ax.plot(rbinm, entropy_bin)
	ax.set_xlim((10.0**rmin,10.0**rmax))	

	plt.savefig("HaloProfiles_"+str(num).zfill(3)+format)
	#plt.show()
        #fig.clf()
