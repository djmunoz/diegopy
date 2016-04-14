import numpy as np

def get_index_list(idarr1,idarr2):
 	"""returns indeces OF idarr2 of elements common in idarr1 and idarr2
	   import numpy as np
	   a=np.sort(np.array([2,8,33,70,100,20]))
           b=np.sort(np.array([2,5,8,90,10,20,30,11,2,33]))
           ind=get_index_list(a,b)
           print b[ind]

	"""
	narr1=idarr1.shape[0]
	narr2=idarr2.shape[0]
	ind_all=np.zeros(narr1,dtype="uint64")

	icountall=0
	icountarr1=0
	icountarr2=0
	lend=0

	iiarr2=np.argsort(idarr2)
	iiarr1=np.argsort(idarr1)

	icountall=0
	icountarr1=0
	icountarr2=0
	lend=0


	while (lend==0):
		if (idarr2[iiarr2[icountarr2]] == idarr1[iiarr1[icountarr1]]):
         		ind_all[icountall]=iiarr2[icountarr2]
         		icountall=icountall+1
         		icountarr2=icountarr2+1
         		icountarr1=icountarr1+1
		else:
         		if (idarr2[iiarr2[icountarr2]] < idarr1[iiarr1[icountarr1]]):
            			icountarr2=icountarr2+1 
         		else:
            			icountarr1=icountarr1+1        
      		if ((icountarr2 >= narr2) | (icountarr1 >= narr1)):
      			lend=1

	if (icountall > 0): 
   		ind_all=ind_all[0:icountall] 
	else:
		ind_all=-1
	
	return ind_all







import numpy as np
from numpy import random

rottable3=np.array([
  [36, 28, 25, 27, 10, 10, 25, 27],
  [29, 11, 24, 24, 37, 11, 26, 26],
  [8, 8, 25, 27, 30, 38, 25, 27],
  [9, 39, 24, 24, 9, 31, 26, 26],
  [40, 24, 44, 32, 40, 6, 44, 6],
  [25, 7, 33, 7, 41, 41, 45, 45],
  [4, 42, 4, 46, 26, 42, 34, 46],
  [43, 43, 47, 47, 5, 27, 5, 35],
  [33, 35, 36, 28, 33, 35, 2, 2],
  [32, 32, 29, 3, 34, 34, 37, 3],
  [33, 35, 0, 0, 33, 35, 30, 38],
  [32, 32, 1, 39, 34, 34, 1, 31],
  [24, 42, 32, 46, 14, 42, 14, 46],
  [43, 43, 47, 47, 25, 15, 33, 15],
  [40, 12, 44, 12, 40, 26, 44, 34],
  [13, 27, 13, 35, 41, 41, 45, 45],
  [28, 41, 28, 22, 38, 43, 38, 22],
  [42, 40, 23, 23, 29, 39, 29, 39],
  [41, 36, 20, 36, 43, 30, 20, 30],
  [37, 31, 37, 31, 42, 40, 21, 21],
  [28, 18, 28, 45, 38, 18, 38, 47],
  [19, 19, 46, 44, 29, 39, 29, 39],
  [16, 36, 45, 36, 16, 30, 47, 30],
  [37, 31, 37, 31, 17, 17, 46, 44],
  [12, 4, 1, 3, 34, 34, 1, 3],
  [5, 35, 0, 0, 13, 35, 2, 2],
  [32, 32, 1, 3, 6, 14, 1, 3],
  [33, 15, 0, 0, 33, 7, 2, 2],
  [16, 0, 20, 8, 16, 30, 20, 30],
  [1, 31, 9, 31, 17, 17, 21, 21],
  [28, 18, 28, 22, 2, 18, 10, 22],
  [19, 19, 23, 23, 29, 3, 29, 11],
  [9, 11, 12, 4, 9, 11, 26, 26],
  [8, 8, 5, 27, 10, 10, 13, 27],
  [9, 11, 24, 24, 9, 11, 6, 14],
  [8, 8, 25, 15, 10, 10, 25, 7],
  [0, 18, 8, 22, 38, 18, 38, 22],
  [19, 19, 23, 23, 1, 39, 9, 39],
  [16, 36, 20, 36, 16, 2, 20, 10],
  [37, 3, 37, 11, 17, 17, 21, 21],
  [4, 17, 4, 46, 14, 19, 14, 46],
  [18, 16, 47, 47, 5, 15, 5, 15],
  [17, 12, 44, 12, 19, 6, 44, 6],
  [13, 7, 13, 7, 18, 16, 45, 45],
  [4, 42, 4, 21, 14, 42, 14, 23],
  [43, 43, 22, 20, 5, 15, 5, 15],
  [40, 12, 21, 12, 40, 6, 23, 6],
  [13, 7, 13, 7, 41, 41, 22, 20]
])

subpix3= np.array([
  [0, 7, 1, 6, 3, 4, 2, 5],
  [7, 4, 6, 5, 0, 3, 1, 2],
  [4, 3, 5, 2, 7, 0, 6, 1],
  [3, 0, 2, 1, 4, 7, 5, 6],
  [1, 0, 6, 7, 2, 3, 5, 4],
  [0, 3, 7, 4, 1, 2, 6, 5],
  [3, 2, 4, 5, 0, 1, 7, 6],
  [2, 1, 5, 6, 3, 0, 4, 7],
  [6, 1, 7, 0, 5, 2, 4, 3],
  [1, 2, 0, 3, 6, 5, 7, 4],
  [2, 5, 3, 4, 1, 6, 0, 7],
  [5, 6, 4, 7, 2, 1, 3, 0],
  [7, 6, 0, 1, 4, 5, 3, 2],
  [6, 5, 1, 2, 7, 4, 0, 3],
  [5, 4, 2, 3, 6, 7, 1, 0],
  [4, 7, 3, 0, 5, 6, 2, 1],
  [6, 7, 5, 4, 1, 0, 2, 3],
  [7, 0, 4, 3, 6, 1, 5, 2],
  [0, 1, 3, 2, 7, 6, 4, 5],
  [1, 6, 2, 5, 0, 7, 3, 4],
  [2, 3, 1, 0, 5, 4, 6, 7],
  [3, 4, 0, 7, 2, 5, 1, 6],
  [4, 5, 7, 6, 3, 2, 0, 1],
  [5, 2, 6, 1, 4, 3, 7, 0],
  [7, 0, 6, 1, 4, 3, 5, 2],
  [0, 3, 1, 2, 7, 4, 6, 5],
  [3, 4, 2, 5, 0, 7, 1, 6],
  [4, 7, 5, 6, 3, 0, 2, 1],
  [6, 7, 1, 0, 5, 4, 2, 3],
  [7, 4, 0, 3, 6, 5, 1, 2],
  [4, 5, 3, 2, 7, 6, 0, 1],
  [5, 6, 2, 1, 4, 7, 3, 0],
  [1, 6, 0, 7, 2, 5, 3, 4],
  [6, 5, 7, 4, 1, 2, 0, 3],
  [5, 2, 4, 3, 6, 1, 7, 0],
  [2, 1, 3, 0, 5, 6, 4, 7],
  [0, 1, 7, 6, 3, 2, 4, 5],
  [1, 2, 6, 5, 0, 3, 7, 4],
  [2, 3, 5, 4, 1, 0, 6, 7],
  [3, 0, 4, 7, 2, 1, 5, 6],
  [1, 0, 2, 3, 6, 7, 5, 4],
  [0, 7, 3, 4, 1, 6, 2, 5],
  [7, 6, 4, 5, 0, 1, 3, 2],
  [6, 1, 5, 2, 7, 0, 4, 3],
  [5, 4, 6, 7, 2, 3, 1, 0],
  [4, 3, 7, 0, 5, 2, 6, 1],
  [3, 2, 0, 1, 4, 5, 7, 6],
  [2, 5, 1, 6, 3, 4, 0, 7]
])

# This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
#  with x,y,z in the range between 0 and 2^bits-1.
def peano_hilbert_key(xa,ya,za, bits):
  	key = np.zeros(len(xa),dtype="uint32");
        for i in range(0,len(xa)):
        	x=xa[i]
        	y=ya[i]
        	z=za[i]               
       		rotation = 0;
		mask = 1 << (bits - 1)
		while (mask >0):
        		if (x & mask):
	                	px=4
			else:     
	                	px=0
	        	if (y & mask):
        	        	py=2
			else:     
        	        	py=0
	        	if (z & mask):
        	        	pz=1
			else:     
        	        	pz=0
			pix =  px | py | pz                                           
		        key[i] <<= 3
		        key[i] |= subpix3[rotation][pix]
		        rotation = rottable3[rotation][pix]
                
        	        mask >>= 1

	return key
#  for(mask = 1 << (bits - 1); mask > 0; mask >>= 1)
#    {
#      unsigned char pix = ((x & mask) ? 4 : 0) | ((y & mask) ? 2 : 0) | ((z & mask) ? 1 : 0);##
#
##      key <<= 3;
#      key |= subpix3[rotation][pix];
#      rotation = rottable3[rotation][pix];
#    }

#  return key;



