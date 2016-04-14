import numpy
import string
import os
import sys


G = 2.959122082855911E-4
c = 1.73144633E+2
Rsun = 4.6491E-3
Rjup = 4.6732617E-4
Msun = 1.0
Mjup = 9.54265748E-4



def write_mercury_files(pos,vel,mass,name,filename):
    
    f = open(filename,'w')
    f.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
    f.write(") Lines beginning with `)' are ignored.\n")
    f.write(")---------------------------------------------------------------------\n")
    f.write(" style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
    f.write(" epoch (in days) = 0.0\n")
    f.write(")---------------------------------------------------------------------\n")
    
    Nbody = pos.shape[0]
    for i in range(0,Nbody):
        #bodylabel = ' M%d'%(i)
        bodylabel = name[i]
        bodylabel+='     m=%23.17e r=3.d0 d=1.1'%(mass[i])
        f.write(bodylabel+"\n")
        
        f.write(' %24.17e %24.17e %24.17e\n'% (pos[i,0],pos[i,1],pos[i,2]))
        f.write(' %24.17e %24.17e %24.17e\n'% (vel[i,0],vel[i,1],vel[i,2]))
        f.write('  0. 0. 0.\n')
    
    f.close()











