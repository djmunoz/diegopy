# code for reading Subfind's subhalo_tab files
# usage e.g.:
#
# import readsubf
# cat = readsubf.subfind_catalog("./m_10002_h_94_501_z3_csf/",63,masstab=True)
# print cat.nsubs
# print "largest halo x position = ",cat.sub_pos[0][0] 

import numpy as np
import os
import sys
 
class hsml_file:
  def __init__(self, basedir, snapnum, swap = False):
    self.filebase = basedir + "/hsmldir_" + str(snapnum).zfill(3) + "/hsml_" + str(snapnum).zfill(3) + "."
 
    print
    print "reading hsml file for snapshot",snapnum,"of",basedir
 
 
    filenum = 0
    doneflag = False
    skip = 0
    while not doneflag:
      curfile = self.filebase + str(filenum)
      
      if (not os.path.exists(curfile)):
        print "file not found:", curfile
        sys.exit()
      
      f = open(curfile,'rb')
              
      Nhsml = np.fromfile(f, dtype=np.uint32, count=1)[0]
      Nprevious = np.fromfile(f, dtype=np.uint32, count=1)[0]
      Ntotal = np.fromfile(f, dtype=np.uint64, count=1)[0]
      Ntask = np.fromfile(f, dtype=np.uint32, count=1)[0]
      
      if swap:
        Nhsml = Nhsml.byteswap()
        Nprevious = Nprevious.byteswap()
        Ntotal = Ntotal.byteswap()
        Ntask = Ntask.byteswap()

      if filenum == 0:
        self.Nhsml = Ntotal
	self.nfiles = Ntask

        self.Hsml = np.empty(Ntotal, dtype=np.float32)
        self.Density = np.empty(Ntotal, dtype=np.float32)
        self.Veldisp = np.empty(Ntotal, dtype=np.float32)
     
      if Nhsml > 0:
        #print skip, Nhsml, self.Nhsml, Ntotal
        locs = slice(skip, skip + Nhsml)
        self.Hsml[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
        self.Density[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
        self.Veldisp[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
        skip += Nhsml
        

      curpos = f.tell()
      f.seek(0,os.SEEK_END)
      if curpos != f.tell(): print "Warning: finished reading before EOF for file",filenum
      f.close()  
      #print 'finished with file number',filenum,"of",ntask
      filenum += 1
      if filenum == self.nfiles: doneflag = True
       
    if swap:
      self.Hsml.byteswap(True)
      self.Density.byteswap(True)
      self.Veldisp.byteswap(True)

    #print
    #print "number of groups =", self.ngroups
    #print "number of subgroups =", self.nsubs
    #if self.nsubs > 0:
    #  print "largest group of length",self.group_len[0],"has",self.group_nsubs[0],"subhalos"
    #  print 
