import numpy as np
import readsnapHDF5 as rs
import RayCast as rc


#absorption coefficient (=opacity kappa)
Ap=1.0
 
def ray_tracing(filename,nframe,target,cam_pos,cam_orient,frame,near,far,xbins,ybins,Box,basename):
	cam_pos = np.array(cam_pos)
	cam_orient = np.array(cam_orient)	
	target = np.array(target)

	print cam_pos,cam_orient,target
	
    	head=rs.snapshot_header(filename)
	pos=rs.read_block(filename, "POS ", parttype=0)
	mass=rs.read_block(filename, "MASS", parttype=0)
        hsmlfac=2.75*1.1
        hsml=hsmlfac*(3.0/4.0/np.pi*rs.read_block(filename, "VOL ", parttype=0))**(1.0/3.0)
	
	b,l,t,r = frame[0], frame[1], frame[2], frame[3]
	
        #translate box so that center is at Box/2.0
	centerx,centery,centerz = Box/2.0,Box/2.0,Box/2.0
        pos[:,0]=Box/2.0 + (pos[:,0]-centerx)
        pos[:,1]=Box/2.0 + (pos[:,1]-centery)
        pos[:,2]=Box/2.0 + (pos[:,2]-centerz)

        for k in range(0,3):
            ind=pos[:,k]>Box
            pos[ind,k]=pos[ind,k]-Box
            ind=pos[:,k]<0
            pos[ind,k]=pos[ind,k]+Box
            print "min/max coord: ",k,pos[:,k].min(), pos[:,k].max()

        #save original values
        x_orig=pos[:,0]
        y_orig=pos[:,1]
        z_orig=pos[:,2]
        hsml_orig=hsml
        mass_orig=mass
        
	#construct homog. transformation matrix
	PS=np.matrix([[(2*near)/(r-l),0,(r+l)/(r-l),0],[0,(2*near)/(t-b),(t+b)/(t-b),0],[0,0,-(far+near)/(far-near),-(2*far*near)/(far-near)],[0,0,-1,0]])
	nvec=-(cam_pos-target)*(-1.0)/np.sqrt((cam_pos[0]-target[0])**2.0 + (cam_pos[1]-target[1])**2.0 + (cam_pos[2]-target[2])**2.0)  #-1 
	temp=np.cross(cam_orient,nvec)
	rvec=temp/np.sqrt(temp[0]**2.0 + temp[1]**2.0 + temp[2]**2.0)
	uvec=np.cross(nvec, rvec) 
	R=np.matrix([[rvec[0],rvec[1],rvec[2],0],[uvec[0],uvec[1],uvec[2],0],[nvec[0],nvec[1],nvec[2],0],[0,0,0,1]])
	T=np.matrix([[1,0,0,-cam_pos[0]],[0,1,0,-cam_pos[1]],[0,0,1,-cam_pos[2]],[0,0,0,1]])

	PSRT=PS*R*T

	#PSRT tranformation: world coordinates -> camera coordinates
	x=PSRT[0,0]*x_orig + PSRT[0,1]*y_orig + PSRT[0,2]*z_orig + PSRT[0,3]*1
	y=PSRT[1,0]*x_orig + PSRT[1,1]*y_orig + PSRT[1,2]*z_orig + PSRT[1,3]*1
	z=PSRT[2,0]*x_orig + PSRT[2,1]*y_orig + PSRT[2,2]*z_orig + PSRT[2,3]*1
	w=PSRT[3,0]*x_orig + PSRT[3,1]*y_orig + PSRT[3,2]*z_orig + PSRT[3,3]*1

        hsml_x=PS[0,0]*hsml_orig + PS[0,1]*hsml_orig + PS[0,2]*hsml_orig + PS[0,3]*1
        hsml_y=PS[1,0]*hsml_orig + PS[1,1]*hsml_orig + PS[1,2]*hsml_orig + PS[1,3]*1
	
	#homog. scaling
	x/=w
	y/=w
	z/=w
	mass=mass_orig
	s=np.abs(w)
	hsml_x/=s
	hsml_y/=s
	hsml_o=hsml_orig

	#clipping in frustum (clip a bit larger for particle contributions outside of frustum)
	index=(np.abs(x) < 1.01)  & (np.abs(y) < 1.01) & (np.abs(z) < 1.01)
	x=x[index]
	y=y[index]
	z=z[index]
	hsml_x=hsml_x[index]
	hsml_y=hsml_y[index]
	hsml_o=hsml_o[index]
	mass=mass[index]

	print "Number of particles in frustum: ", mass.shape[0]
	
	#sort descending according to pseudo-depth
	index=np.argsort(z)[::-1]
	x=x[index]
	y=y[index]
	z=z[index]
	hsml_x=hsml_x[index]
	hsml_y=hsml_y[index]
        hsml_o=hsml_o[index]
	mass=mass[index]

	#avoid single pixel flickering
        pixfac=0.5
        hsml_x[hsml_x<pixfac*2.0/xbins]=0
        hsml_y[hsml_y<pixfac*2.0/ybins]=0

	#now ray-trace
	print "start render..."

        image=rc.Render(x, y, np.repeat(Ap, x.shape[0]).astype("float32"), hsml_x, hsml_y, hsml_o/(hsmlfac*(3.0/4.0/np.pi)**(1.0/3.0)), mass, xbins, ybins, hsmlfac)
	print "done."

	#save file
	fd=open(basename+str(nframe).zfill(4)+".dat", "wb")
	image.astype("float64").tofile(fd)
	fd.close()
	print image.shape	
	#clean image
	image=image*0.0	
