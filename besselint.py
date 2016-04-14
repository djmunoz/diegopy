import inspect
import pylab as m
import numpy as np
from scipy.special import jn
from scipy.linalg import solve as linearsolve

filename =  inspect.getfile(inspect.currentframe())


def bessel0int(func, a, b, r, args=(), maxpoints=12,abstol=1.49e-8,reltol=1.49e-8):
    
    def u(x,k): return (x - 0.5 * (a + b) - abstol)**(k - 1) 
    def uprime(x,k): return (k - 1) * (x - 0.5 * (a + b) - abstol)**(k - 2)
    vecU=np.vectorize(u)
    vecUprime=np.vectorize(uprime)

    default_points = 5
    vecfunc= np.vectorize(func)

    points = default_points - 1
    old_sol = 1.0
    tol = 1.0e-6
    
    while(True):
        n = points
        nodes = np.linspace(a,b,n)
        fvals = vecfunc(nodes)
        gvals = np.zeros(n)
        vals = np.append(fvals,gvals)
                
        matrix = np.zeros((2*n,2*n))
        for k in range(0,n):
            matrix[0:n,k] =  vecUprime(nodes,k)
            matrix[0:n,k+n] =  r * vecU(nodes,k)
            matrix[n:2*n,k] =  -r * vecU(nodes,k)
            matrix[n:2*n,k+n] =  vecUprime(nodes,k) - vecU(nodes,k) / nodes

        coeff = linearsolve(matrix,vals)

        sol = + (coeff[0:n]*vecU(b,np.arange(1,n+1))).sum()*jn(0,r*b)  \
              - (coeff[0:n]*vecU(a,np.arange(1,n+1))).sum()*jn(0,r*a)  \
              + (coeff[n:2*n]*vecU(b,np.arange(1,n+1))).sum()*jn(1,r*b)  \
              - (coeff[n:2*n]*vecU(a,np.arange(1,n+1))).sum()*jn(1,r*a)  \

        points = points + 1
        if (points > default_points):
            abserr,relerr = np.abs(old_sol-sol),np.abs(old_sol-sol)/np.abs(old_sol)
            if (np.abs(sol) > abstol/reltol):
                if (abserr < abstol): break
            else:
                if (relerr < reltol): break

        if (points > maxpoints):
            print 'File','"'+filename+'"'+', line ',inspect.currentframe().f_back.f_lineno 
            print 'AccuracyWarning: maxpoints(',maxpoints,') exceeded. Latest difference =', abserr ,'AccuracyWarning)'
            break
            
        old_sol = sol
    return sol

