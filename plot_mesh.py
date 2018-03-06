import matplotlib.pyplot as plt
import numpy as np
import readsnapHDF5 as rs
import os
from string import split
import sys

#path to output files
current_dir = split(os.getcwd(),"/")[-1]
snap_base="snap_"
directory='./'


if __name__ == '__main__':

    base = sys.argv[1]

    X0 = float(sys.argv[2])
    X1 = float(sys.argv[3])
    Y0 = float(sys.argv[4])
    Y1 = float(sys.argv[5])
        
    snap_list = np.arange(0,1000)

    time_list=[]
    for num in snap_list:
        filename=directory+base+snap_base+str(num).zfill(3)
	#open the snapshot header
	header = rs.snapshot_header(filename)
	time_list.append(header.time)  
        
        Lx = header.boxsize
        Ly = header.boxsize
	####################################################################
	#open mesh data
	filename=base+"voronoi_mesh_"+str(num).zfill(3)
	f = open(filename,'rb')
	ngas=np.fromfile(f,dtype=np.int32,count=1)
	nel=np.fromfile(f,dtype=np.int32,count=1)
	nedgepoints=np.fromfile(f,dtype=np.int32,count=1)
	nedges=np.fromfile(f,dtype=np.int32,count=ngas[0])
	nedges_offest=np.fromfile(f,dtype=np.int32,count=ngas[0])
	edgelist=np.fromfile(f,dtype=np.int32,count=nel[0])
	points=np.fromfile(f,dtype=np.dtype((np.float32,2)),count=nedgepoints[0])
	f.close()
	####################################################################
	#plot
	
	fig = plt.figure(1, figsize=(7.0,6.0))
	fig.subplots_adjust(top=0.98,bottom=0.08,right=0.95,left=0.05, hspace=0.0)
	
	########################
	ax=fig.add_subplot(1, 1, 1)
	
	
	########################
        #plot the mesh as polygons
	for i in range(0,ngas[0]):
		x = points[edgelist[nedges_offest[i]:nedges_offest[i]+nedges[i]],0]
	        y = points[edgelist[nedges_offest[i]:nedges_offest[i]+nedges[i]],1]
		ax.fill(np.append(x,x[0]), np.append(y,y[0]), edgecolor='k', lw=0.2, fill=False)
	        if np.shape(np.nonzero((x < 0.0) | (x > Lx) | (y < 0.0) | (y > Ly )))[1] > 0:
        		for dx in range(-1,2):
	               		    for dy in range(-1,2):
					    ax.fill(np.append(x,x[0])+dx*Lx,np.append(y,y[0])+dy*Ly,
                                                    lw=0.2,edgecolor='k', fill=False)
        ax.set_xlabel("$x$", fontsize=24)
	ax.set_ylabel("$y$", fontsize=24)
        ax.set_xlim(X0,X1)
        ax.set_ylim(Y0,Y1)
        ax.grid(True)
        ########################
	#write image (change output format by changing filename extension)
	fig.savefig("mesh_"+str(num).zfill(3)+".png")
	fig.clf()
