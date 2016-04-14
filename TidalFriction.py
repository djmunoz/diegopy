import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
import mercury
import numpy.random as ran

G = mercury.G
clight = mercury.c
Rsun = mercury.Rsun
Mjup = mercury.Mjup
Rjup = mercury.Rjup


def binary_evolution_elements(y,t,m1,m2,m3,a_out,e_out,R1,R2):
    #time-dependent variables
    ein = y[0]
    ain = y[1]
    I = y[2]
    omega = y[3]
    Node = y[4]
    Omega1x = y[5]
    Omega1y = y[6]
    Omega1z = y[7]
    Omega2x = y[8]
    Omega2y = y[9]
    Omega2z = y[10]


    #some local definitions for convenience
    einx = (np.cos(Node)*np.cos(omega) - np.sin(Node) * np.sin(omega) * np.cos(I))  # * ein
    einy = (np.sin(Node)*np.cos(omega) + np.cos(Node) * np.sin(omega) * np.cos(I)) # * ein
    einz = np.sin(omega) * np.sin(I) # * ein

    hin= np.sqrt(ain * G * (m1 + m2) * (1 - ein**2))
    hinx =  np.sin(Node) * np.sin(I) #* hin
    hiny = -np.cos(Node) * np.sin(I) #* hin
    hinz =  np.cos(I) #* hin

    qinx = -np.sin(omega) * np.cos(Node) - np.cos(omega) * np.sin(Node) * np.cos(I)
    qiny = -np.sin(omega) * np.sin(Node) + np.cos(omega) * np.cos(Node) * np.cos(I)
    qinz = np.cos(omega) * np.sin(I)

    #qinx = hiny * einz - hinz * einy
    #qiny = hinz * einx - hinx * einz
    #qinz = hinx * einy - hiny * einx
    #qin = np.sqrt(qinx**2 + qiny**2 + qinz**2)

    Omega1_e = (Omega1x * einx + Omega1y * einy + Omega1z * einz)# / ein
    Omega1_h = (Omega1x * hinx + Omega1y * hiny + Omega1z * hinz)# / hin
    Omega1_q = (Omega1x * qinx + Omega1y * qiny + Omega1z * qinz)# / qin
    
    Omega2_e = (Omega2x * einx + Omega2y * einy + Omega2z * einz)# / ein
    Omega2_h = (Omega2x * hinx + Omega2y * hiny + Omega2z * hinz)# / hin
    Omega2_q = (Omega2x * qinx + Omega2y * qiny + Omega2z * qinz)#/ qin

    Pin = 2 * np.pi * np.sqrt(ain * ain * ain /G / (m1 + m2))
    Pout = 2 * np.pi * np.sqrt(a_out * a_out * a_out /G / (m1 + m2 + m3))
    mu = m1 * m2 / (m1 +  m2)

    tau = 2* Pout**2 /3 /np.pi / Pin * (m1 + m2 + m3)/m3 *  (1 - e_out**2)**1.5

    C = 1.0/3/tau / np.sqrt(1 - ein**2)

    houtx = 0
    houty = 0
    houtz = 1

    Seq = C  * (-3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * qinx + houty * qiny + houtz * qinz) )#/ ein/ qin)

    See = C  * (1 - 3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * einx + houty * einy + houtz * einz))#  / ein/ ein)

    Sqq = C  * (1 - 3 * (houtx * qinx + houty * qiny + houtz * qinz) * 
                (houtx * qinx + houty * qiny + houtz * qinz))# / qin/ qin)

    Sqh = C  * (-3 * (houtx * qinx + houty * qiny + houtz * qinz) * 
                (houtx * hinx + houty * hiny + houtz * hinz))# / qin/ hin)

    Seh = C  * (-3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * hinx + houty * hiny + houtz * hinz))# / ein/ hin)
    
  
    Qstar = 0.028
    Qplanet = 0.51
    rgstar = 0.08
    rgplanet = 0.25
    tvstar =   2.0e4 # in days, about 55 years
    tvplanet = 0.365242  # in days, about 0.001 years

    Q1, Q2 = Qstar, Qplanet

    tv1, tv2 = tvstar, tvplanet
    
    rg1, rg2 = rgstar, rgplanet

    I1 = rg1 * m1 * R1**2
    I2 = rg2 * m2 * R2**2
    
    k1 = 0.5 * Q1/(1 - Q1)
    k2 = 0.5 * Q2/(1 - Q2)

    tf1 = tv1/9 * (ain/R1)**8 * m1**2 / ((m1 + m2)*m2) / (1 + 2 * k1)**2 

    tf2 = tv2/9 * (ain/R2)**8 * m2**2 / ((m1 + m2)*m1) / (1 + 2 * k2)**2 

    #print tf1/Pin, tf2/Pin

    V1 = 9.0 / tf1 * ((1 + 15.0/4 * ein**2 + 15.0/8 * ein**4 + 5.0/64 * ein**6)/(1 - ein**2)**6.5 - \
                      11.0/18 * Omega1_h* Pin/2/np.pi * (1 + 1.5 * ein**2 + 0.125 * ein**4)/( 1 - ein**2)**5)
    
    W1 = 1.0 / tf1 * ((1 + 15.0/2 * ein**2 + 45.0/8 * ein**4 + 5.0/16 * ein**6)/(1 - ein**2)**6.5 - \
                      Omega1_h* Pin/2/np.pi * (1 + 3 * ein**2 + 0.375 * ein**4)/( 1 - ein**2)**5)    
    
    X1 = Pin/2/np.pi *(-m2 * k1 * (R1/ain)**5 / mu  *  Omega1_h *  Omega1_e / (1 - ein**2)**2 \
                       -Omega1_q /2/tf1 * (1 + 4.5 * ein**2 + 0.625 * ein**4) / (1 - ein**2)**5)
    
    Y1 = Pin/2/np.pi *(-m2 * k1 * (R1/ain)**5 / mu  *  Omega1_h *  Omega1_q / (1 - ein**2)**2 \
                       + Omega1_e /2/tf1 * (1 + 1.5 * ein**2 + 0.125 * ein**4) / (1 - ein**2)**5)                    
    
    Z1 = Pin/2/np.pi * m2 * k1 * (R1/ain)**5 /mu * (0.5*(2 * Omega1_h**2 - Omega1_e**2 - Omega1_q**2) / (1 - ein**2)**2 \
                                                    +15*G * m2 / ain**3 *(1 + 1.5*ein**2 + 0.125*ein**4) / (1 - ein**2)**5)


    V2 = 9.0 / tf2 * ((1 + 15.0/4 * ein**2 + 15.0/8 * ein**4 + 5.0/64 * ein**6)/(1 - ein**2)**6.5 - \
                      11.0/18 * Omega2_h* Pin/2/np.pi * (1 + 1.5 * ein**2 + 0.125 * ein**4)/( 1 - ein**2)**5)

    W2 = 1.0 / tf2 * ((1 + 15.0/2 * ein**2 + 45.0/8 * ein**4 + 5.0/16 * ein**6)/(1 - ein**2)**6.5 - \
                      Omega2_h* Pin/2/np.pi * (1 + 3 * ein**2 + 0.375 * ein**4)/( 1 - ein**2)**5)    
    
    X2 = Pin/2/np.pi *(-m2 * k2 * (R2/ain)**5 / mu  *  Omega2_h *  Omega2_e / (1 - ein**2)**2 \
                       -Omega2_q /2/tf2 * (1 + 4.5 * ein**2 + 0.625 * ein**4) / (1 - ein**2)**5)
    
    Y2 = Pin/2/np.pi *(-m2 * k2 * (R2/ain)**5 / mu  *  Omega2_h *  Omega2_q / (1 - ein**2)**2 \
                       +Omega2_e /2/tf2 * (1 + 1.5 * ein**2 + 0.125 * ein**4) / (1 - ein**2)**5)                    
    
    Z2 = Pin/2/np.pi * m2 * k2 * (R2/ain)**5 /mu * (0.5*(2 * Omega2_h**2 - Omega2_e**2 - Omega2_q**2) / (1 - ein**2)**2 \
                                                    +15*G * m2 / ain**3 *(1 + 1.5*ein**2 + 0.125*ein**4) / (1 - ein**2)**5)


    ZGR = 3 * G * (m1 + m2) * 2 * np.pi / Pin / ain / clight / clight / (1 - ein**2)
    
    #X1, X2, Y1, Y2, V1, V2, W1, W2, Z1, Z2, ZGR = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    
    #11 differential equations for 11 time-dependent quantities

    dein_dt = ein * (-V1 -V2 - 5*(1 - ein**2) * Seq)
    
    dain_dt = ain * (-2*W1 - 2*W2 - 2*ein**2/(1-ein**2) * (V1 + V2))

    dI_dt = -np.sin(omega) * (Y1 + Y2 + (1-ein**2) * Sqh) + np.cos(omega) * (X1 + X2 + (4*ein**2 + 1) * Seh)

    domega_dt = (Z1 + Z2 + ZGR + (1 - ein**2)*(4*See - Sqq))  \
                - np.cos(omega) * np.cos(I)/np.sin(I) * (Y1 + Y2 + (1 - ein**2) * Sqh) \
                - np.sin(omega) * np.cos(I)/np.sin(I) * (X1 + X2 + (4 * ein**2 + 1) * Seh)
    
    dNode_dt = np.cos(omega) /np.sin(I) * (Y1 + Y2 + (1 - ein**2) * Sqh) \
               + np.sin(omega) / np.sin(I) * (X1 + X2 + (4 * ein**2 + 1) * Seh)
    
    dOmega1x_dt = mu * hin / I1 * (- Y1 * einx + X1 * qinx + W1 * hinx)

    dOmega1y_dt = mu * hin / I1 * (- Y1 * einy + X1 * qiny + W1 * hiny)

    dOmega1z_dt = mu * hin / I1 * (- Y1 * einz + X1 * qinz + W1 * hinz)

    dOmega2x_dt = mu * hin / I2 * (- Y2 * einx + X2 * qinx + W2 * hinx)

    dOmega2y_dt = mu * hin / I2 * (- Y2 * einy + X2 * qiny + W2 * hiny)

    dOmega2z_dt = mu * hin / I2 * (- Y2 * einz + X2 * qinz + W2 * hinz)

    return [dein_dt, dain_dt, dI_dt, domega_dt, dNode_dt, dOmega1x_dt, dOmega1y_dt, dOmega1z_dt, dOmega2x_dt, dOmega2y_dt, dOmega2z_dt]


def kozai_tidal_friction(m1,m2,m3,a,e,I,omega,Node,Spin1,Spin2,R1,R2,a_out,e_out,time):
    
    
    #initial conditions
    ein_init = e
    ain_init = a
    I_init = I 
    omega_init = omega 
    Node_init  = Node 
    Omega1_init =  2*np.pi/Spin1
    Omega2_init =  2*np.pi/Spin2

    #project the stellar spins (assumed to have zero obliquity at initial time)
    Omega1x_init =  Omega1_init * np.sin(Node_init) * np.sin(I_init)
    Omega1y_init = -Omega1_init * np.cos(Node_init) * np.sin(I_init) 
    Omega1z_init =  Omega1_init * np.cos(I_init)
    Omega2x_init =  Omega2_init * np.sin(Node_init) * np.sin(I_init)
    Omega2y_init = -Omega2_init * np.cos(Node_init) * np.sin(I_init) 
    Omega2z_init =  Omega2_init * np.cos(I_init)


    ics = [ein_init, ain_init, I_init, omega_init, Node_init, Omega1x_init, Omega1y_init, Omega1z_init, Omega2x_init, Omega2y_init, Omega2z_init]

    #solve ODE
    t_years = time
    #t_years=np.append(0.0,np.logspace(6.0,9.3,500))

    t=t_years*365.25


    sol=integ.odeint(binary_evolution_elements,ics,t,args=(m1,m2,m3,a_out,e_out,R1,R2),mxords=16,mxordn=15,hmin=10.0)
    
    return sol
