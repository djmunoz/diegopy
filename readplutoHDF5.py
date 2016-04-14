# import readsnapHDF5 as rs
# header = rs.snapshot_header("snap_063.0") 
# mass = rs.read_block("snap_063","MASS",parttype=5) # reads mass for particles of type 5

import numpy as np
import os
import sys
import math
import tables


############ 
#DATABLOCKS#
############
#descriptions of all datablocks -> add datablocks here!
#TAG:[HDF5_NAME,DIM]
datablocks = {"RHO ":["rho",1], 
	      "VELX":["vx1",1],
	      "VELY":["vx2",1],
	      "VELZ":["vx3",1],
              }
#####################################################################################################################
#                                                    READING ROUTINES			                            #
#####################################################################################################################


########################### 
#CLASS FOR SNAPSHOT HEADER#
###########################  
class snapshot_header:
	def __init__(self, *args, **kwargs):
		if (len(args) == 1):
			filename = args[0]
 			if os.path.exists(filename):
			        curfilename=filename
			elif os.path.exists(filename+".h5"):
				curfilename = filename+".h5"
			elif os.path.exists(filename+".0.h5"): 
				curfilename = filename+".0.h5"
			else:	
				print "[error] file not found : ", filename
			    	sys.exit()
    
                        timestep="Timestep_"+("%s"  % int(curfilename[curfilename.find("data.")+5:curfilename.find("data.")+9]))
			f=tables.openFile(curfilename)
                        self.timestep = timestep
			self.time = f.root._f_getChild(timestep)._v_attrs.Time 
			f.close()

		else:
			#read arguments
			self.time = kwargs.get("Time")

			#set default values
                        if (self.time == None):				
				self.time = np.array([0], dtype="float64")
			
		        if (self.time == None):				
				self.timestep = np.array(["Timestep_0"])

##############################
#READ ROUTINE FOR SINGLE FILE#
############################## 
def read_block_single_file(filename, block_name, fill_block_name="", verbose=False):

  if (verbose):
	  print "[single] reading file           : ", filename   	
	  print "[single] reading                : ", block_name
      
  head = snapshot_header(filename)
  timestep=head.timestep
  del head

  f=tables.openFile(filename)

  print timestep
  f.root._f_getChild(timestep)._f_getChild("vars")

  if (f.root._f_getChild(timestep).__contains__("vars")):
                
        #default: just read the block
        if (f.root._f_getChild(timestep)._f_getChild("vars").__contains__(block_name)):
                data=f.root._f_getChild(timestep)._f_getChild("vars")._f_getChild(block_name)[:]
                ret_val=data
          
  f.close()
        
  return ret_val

##############
#READ ROUTINE#
##############
def read_block(filename, block, fill_block="", verbose=False):
  if (verbose):
          print "reading block          : ", block

  curfilename=filename

  if os.path.exists(curfilename):
    multiple_files=False
  elif os.path.exists(filename+".0"+".h5"):
    curfilename = filename+".0"+".h5"
    multiple_files=True
  else:
    print "[error] file not found : ", filename
    sys.exit()

  head = snapshot_header(curfilename)

 
  if (datablocks.has_key(block)):
        block_name=datablocks[block][0]
        first=True
        if (verbose):
                print "Reading HDF5           : ", block_name
  else:
        print "[error] Block type ", block, "not known!"
        sys.exit()

  fill_block_name=""

  ret_val=read_block_single_file(curfilename, block_name, fill_block_name, verbose)

  return ret_val


#############
#LIST BLOCKS#
#############
def list_blocks(filename, parttype=-1, verbose=False):
  
  f=tables.openFile(filename)
  for parttype in range(0,5):
  	part_name='PartType'+str(parttype)
        if (f.root.__contains__(part_name)):
        	print "Parttype contains : ", parttype
		print "-------------------"
		iter = it=datablocks.__iter__()
		next = iter.next()
		while (1):
			if (verbose):
				print "check ", next, datablocks[next][0]
			if (f.root._f_getChild(part_name).__contains__(datablocks[next][0])):
  				print next, datablocks[next][0]
			try:
				next=iter.next()
			except StopIteration:
				break	
  f.close() 

#################
#CONTAINS BLOCKS#
#################
def contains_block(filename, tag, parttype=-1, verbose=False):
  
  contains_flag=False
  f=tables.openFile(filename)
  for parttype in range(0,5):
        part_name='PartType'+str(parttype)
        if (f.root.__contains__(part_name)):
                iter = it=datablocks.__iter__()
                next = iter.next()
                while (1):
                        if (verbose):
                                print "check ", next, datablocks[next][0]
                        if (f.root._f_getChild(part_name).__contains__(datablocks[next][0])):
                                if (next.find(tag)>-1):
					contains_flag=True	
                        try:
                                next=iter.next()
                        except StopIteration:
                                break
  f.close() 
  return contains_flag

############
#CHECK FILE#
############
def check_file(filename):
  f=tables.openFile(filename)
  f.close()
                                                                                                                                                  






#####################################################################################################################
#                                                    WRITING ROUTINES    		                            #
#####################################################################################################################



#######################
#OPEN FILE FOR WRITING#
#######################
def openfile(filename):
	f=tables.openFile(filename, mode = "w")	 
	return f

############
#CLOSE FILE#
############
def closefile(f):
	f.close()

##############################
#WRITE SNAPSHOT HEADER OBJECT#
##############################
def writeheader(f, header):	
    	group_header=f.createGroup(f.root, "Header")
	group_header._v_attrs.NumPart_ThisFile=header.npart
	group_header._v_attrs.NumPart_Total=header.nall
	group_header._v_attrs.NumPart_Total_HighWord=header.nall_highword
	group_header._v_attrs.MassTable=header.massarr
	group_header._v_attrs.Time=header.time
	group_header._v_attrs.Redshift=header.redshift
	group_header._v_attrs.BoxSize=header.boxsize
	group_header._v_attrs.NumFilesPerSnapshot=header.filenum						
	group_header._v_attrs.Omega0=header.omega0							
	group_header._v_attrs.OmegaLambda=header.omegaL							
	group_header._v_attrs.HubbleParam=header.hubble	
	group_header._v_attrs.Flag_Sfr=header.sfr	
	group_header._v_attrs.Flag_Cooling=header.cooling
	group_header._v_attrs.Flag_StellarAge=header.stellar_age		
	group_header._v_attrs.Flag_Metals=header.metals		
	group_header._v_attrs.Flag_Feedback=header.feedback				
	group_header._v_attrs.Flag_DoublePrecision=header.double			

###############
#WRITE ROUTINE#
###############
def write_block(f, block, parttype, data):
	part_name="PartType"+str(parttype)
	if (f.root.__contains__(part_name)==False):
	    	group=f.createGroup(f.root, part_name)
	else:
		group=f.root._f_getChild(part_name)	
	
	if (datablocks.has_key(block)):
        	block_name=datablocks[block][0]
	        dim2=datablocks[block][1]		
		if (group.__contains__(block_name)==False):
			table=f.createArray(group, block_name, data)
		else:
			print "I/O block already written"
	else:
		print "Unknown I/O block"		

	



