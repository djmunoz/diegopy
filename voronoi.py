from pylab import *
import matplotlib.pyplot as plt
import numpy as np

import conversions as co
import readsnap as rs


def visualize_mesh(basename, num, Lx, Ly, mi, ma, arepo=1, format=".eps"):	
	filename=basename+str(num).zfill(3)
        print filename
	rho = rs.read_block(filename,"RHO ", parttype=0, arepo=arepo)

        filename='voronoi_mesh_'+str(num).zfill(3)
        f = open(filename,'rb')

        ngas=np.fromfile(f,dtype=np.int32,count=1)
        nel=np.fromfile(f,dtype=np.int32,count=1)
        nedgepoints=np.fromfile(f,dtype=np.int32,count=1)
        nedges=np.fromfile(f,dtype=np.int32,count=ngas[0])
        nedges_offest=np.fromfile(f,dtype=np.int32,count=ngas[0])
        edgelist=np.fromfile(f,dtype=np.int32,count=nel[0])
        points=np.fromfile(f,dtype=np.dtype((np.float32,2)),count=nedgepoints[0])
        f.close()

        print "Number        = ", num
        print "Voronoi Cells = ", ngas[0]
        print "nel           = ", nel[0]
        print "Edge points   = ", nedgepoints[0]

        if (Lx>Ly):
	        fig = plt.figure(1, figsize=(Lx/Ly*2.5,2.5))
	else:
		fig = plt.figure(1, figsize=(10.0,10.0))
        ax = fig.add_subplot(1,1,1)

        print "min/max of rho = ", min(rho),max(rho)

        rho[rho > ma] = ma
        rho[rho < mi] = mi



        for i in range(0,ngas[0]):
		x = points[edgelist[nedges_offest[i]:nedges_offest[i]+nedges[i]],0]
		y = points[edgelist[nedges_offest[i]:nedges_offest[i]+nedges[i]],1]
                cmap=get_cmap("jet")
                co = cmap( (rho[i]-mi)/(ma-mi) )


		ax.fill(np.append(x,x[0]), np.append(y,y[0]), color=co, edgecolor='white')

	        if shape(np.nonzero((x < 0.0) | (x > Lx) | (y < 0.0) | (y > Ly )))[1] > 0:
		        for dx in range(-1,2):
			        for dy in range(-1,2):        
                      		 	ax.fill(np.append(x,x[0])+dx*Lx,np.append(y,y[0])+dy*Ly,color=co, edgecolor='white')
	                	
        
	ax.axis([0.0,Lx,0.0,Ly])
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.xaxis.set_major_locator(MultipleLocator(0.5))
	ax.yaxis.set_major_locator(MultipleLocator(0.5))

	plt.savefig("image_"+str(num).zfill(3)+format)
        clf()

    
def cell_mass_dist(basename, num, arepo=1, bins=100, format=".eps"):	
	filename=basename+str(num).zfill(3)

	mass = rs.read_block(filename,"MASS", parttype=0,arepo=arepo)

	meanmass=mass.mean()

	print "mean cell mass       = ", meanmass
        print "min/max of cell mass = ", mass.min(),mass.max()


        fig = plt.figure(1, figsize=(10.0,10.0))
        ax = fig.add_subplot(1,1,1)

	# the histogram of the data
	n, bins, patches = plt.hist(np.log10(mass/meanmass), bins, normed=1, facecolor='blue')

	ax.set_xlabel('log[m/<m>]')
	ax.set_ylabel('df / dlog[m/<m>]')

	plt.savefig("cell_mass_dist_"+str(num).zfill(3)+format)
        clf()

    
