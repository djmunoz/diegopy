import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plotxy(x,y):
	fig = plt.figure(1, figsize=(10.0,10.0))
	plt.plot(x,y)
	plt.show()

def scatterxy(x,y,s=1,c='black',xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0):
	fig = plt.figure(1, figsize=(10.0,10.0))
	ax=fig.add_subplot(1,1,1)
        ax.scatter(x,y,s=s, c=c)
	if (xmin!=0.0) | (xmax!=0.0) | (ymin!=0.0) | (ymax!=0.0):
		ax.set_xlim([xmin,xmax])
		ax.set_ylim([ymin,ymax])
        plt.show()

def scatterxyz(x,y,z,s=1,c='black'):
        fig = plt.figure(1, figsize=(10.0,10.0))
	ax = Axes3D(fig)
    	ax.scatter(x, y, z, s=s, c=c)
	plt.show()

