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

#parameters
m1 = 1.1
m2 = 7.8 * Mjup
m3 = 1.1 #0.00001
R1 = Rsun
R2 = Rjup
a_out = 1000.0
e_out = 0.5


def secular_ode(y,t):
    #time-dependent variables
    einx = y[0]
    einy = y[1]
    einz = y[2]
    hinx = y[3]
    hiny = y[4]
    hinz = y[5]
    Omega1x = y[6]
    Omega1y = y[7]
    Omega1z = y[8]
    Omega2x = y[9]
    Omega2y = y[10]
    Omega2z = y[11]
    houtx = y[12]
    houty = y[13]
    houtz = y[14]


    #some local definitions for convenience
    qinx = hiny * einz - hinz * einy
    qiny = hinz * einx - hinx * einz
    qinz = hinx * einy - hiny * einx
    

    ein = np.sqrt(einx**2 + einy**2 + einz**2)
    hin = np.sqrt(hinx**2 + hiny**2 + hinz**2)
    qin = np.sqrt(qinx**2 + qiny**2 + qinz**2)
    hout = np.sqrt(houtx**2 + houty**2 + houtz**2)

    Omega1_e = (Omega1x * einx + Omega1y * einy + Omega1z * einz) / ein
    Omega1_h = (Omega1x * hinx + Omega1y * hiny + Omega1z * hinz) / hin
    Omega1_q = (Omega1x * qinx + Omega1y * qiny + Omega1z * qinz) / qin
    
    Omega2_e = (Omega2x * einx + Omega2y * einy + Omega2z * einz) / ein
    Omega2_h = (Omega2x * hinx + Omega2y * hiny + Omega2z * hinz) / hin
    Omega2_q = (Omega2x * qinx + Omega2y * qiny + Omega2z * qinz) / qin


    ain = hin**2 / G /(m1 + m2) / (1 - ein**2)
    aout = hout**2 / G /(m1 + m2 + m3)  / (1 - e_out**2)
    Pin = 2 * np.pi * np.sqrt(ain * ain * ain /G / (m1 + m2))
    Pout = 2 * np.pi * np.sqrt(a_out * a_out * a_out /G / (m1 + m2 + m3))
    mu = m1 * m2 / (m1 +  m2)

    tau = 2* Pout**2 /3 /np.pi / Pin * (m1 + m2 + m3)/m3 *  (1 - e_out**2)**1.5

    C = 1.0/3/tau / np.sqrt(1 - ein**2)

    Seq = C  * (-3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * qinx + houty * qiny + houtz * qinz) /hout**2 / ein/ qin)

    See = C  * (1 - 3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * einx + houty * einy + houtz * einz) /hout**2 / ein/ ein)

    Sqq = C  * (1 - 3 * (houtx * qinx + houty * qiny + houtz * qinz) * 
                (houtx * qinx + houty * qiny + houtz * qinz)/hout**2 / qin/ qin)

    Sqh = C  * (-3 * (houtx * qinx + houty * qiny + houtz * qinz) * 
                (houtx * hinx + houty * hiny + houtz * hinz)/hout**2 / qin/ hin)

    Seh = C  * (-3 * (houtx * einx + houty * einy + houtz * einz) * 
                (houtx * hinx + houty * hiny + houtz * hinz) /hout**2 / ein/ hin)
    
  
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

    print ain

    V1 = 9.0 / tf1 * ((1 + 15.0/4 * ein**2 + 15.0/8 * ein**4 + 5.0/64 * ein**6) / \
                      (1 - ein**2)**6.5 - \
                      11.0/18 * Omega1_h* Pin/2/np.pi * (1 + 1.5 * ein**2 + 0.125 * ein**4)/( 1 - ein**2)**5)
    
    W1 = 1.0 / tf1 * ((1 + 15.0/2 * ein**2 + 45.0/8 * ein**4 + 5.0/16 * ein**6)/ \
                      (1 - ein**2)**6.5 - \
                      Omega1_h* Pin/2/np.pi * (1 + 3 * ein**2 + 0.375 * ein**4)/( 1 - ein**2)**5)    
    
    X1 = Pin/2/np.pi *(-m2 * k1 * (R1/ain)**5 / mu  *  Omega1_h *  Omega1_e / (1 - ein**2)**2 \
                       -Omega1_q /2/tf1 * (1 + 4.5 * ein**2 + 0.625 * ein**4) / (1 - ein**2)**5)
    
    Y1 = Pin/2/np.pi *(-m2 * k1 * (R1/ain)**5 / mu  *  Omega1_h *  Omega1_q / (1 - ein**2)**2 \
                       +Omega1_e /2/tf1 * (1 + 1.5 * ein**2 + 0.125 * ein**4) / (1 - ein**2)**5)                    
    
    Z1 = Pin/2/np.pi * m2 * k1 * (R1/ain)**5 /mu * (0.5*(2 * Omega1_h**2 - Omega1_e**2 - Omega1_q**2) / (1 - ein**2)**2 \
                                                    +15*G * m2 / ain**3 *(1 + 1.5*ein**2 + 0.125*ein**4) / (1 - ein**2)**5)


    V2 = 9.0 / tf2 * ((1 + 15.0/4 * ein**2 + 15.0/8 * ein**4 + 5.0/64 * ein**6) / \
                      (1 - ein**2)**6.5 - \
                      11.0/18 * Omega2_h* Pin/2/np.pi * (1 + 1.5 * ein**2 + 0.125 * ein**4)/( 1 - ein**2)**5)

    W2 = 1.0 / tf2 * ((1 + 15.0/2 * ein**2 + 45.0/8 * ein**4 + 5.0/16 * ein**6)/ \
                      (1 - ein**2)**6.5 - \
                      Omega2_h* Pin/2/np.pi * (1 + 3 * ein**2 + 0.375 * ein**4)/( 1 - ein**2)**5)    
    
    X2 = Pin/2/np.pi *(-m2 * k2 * (R2/ain)**5 / mu  *  Omega2_h *  Omega2_e / (1 - ein**2)**2 \
                       -Omega2_q /2/tf2 * (1 + 4.5 * ein**2 + 0.625 * ein**4) / (1 - ein**2)**5)
    
    Y2 = Pin/2/np.pi *(-m2 * k2 * (R2/ain)**5 / mu  *  Omega2_h *  Omega2_q / (1 - ein**2)**2 \
                       +Omega2_e /2/tf2 * (1 + 1.5 * ein**2 + 0.125 * ein**4) / (1 - ein**2)**5)                    
    
    Z2 = Pin/2/np.pi * m2 * k2 * (R2/ain)**5 /mu * (0.5*(2 * Omega2_h**2 - Omega2_e**2 - Omega2_q**2) / (1 - ein**2)**2 \
                                                    +15*G * m2 / ain**3 *(1 + 1.5*ein**2 + 0.125*ein**4) / (1 - ein**2)**5)


    ZGR = 3 * G * (m1 + m2) * 2 * np.pi / Pin / ain / clight / (1 - ein**2)

    
    #print X1,Y1,Z1, V1,W1
    
    #Fifteen differential equations for 12 time-dependent quantities

    deinx_dt = ein * ((Z1 + Z2 + ZGR) * qinx/qin - (Y1 + Y2) * hinx/hin - (V1 + V2) * einx/ein  \
                    - (1 - ein**2) * (5 * Seq * einx/ein - (4 * See - Sqq) * qinx/qin +  Sqh * hinx/hin))
    
    deiny_dt = ein * ((Z1 + Z2 + ZGR) * qiny/qin - (Y1 + Y2) * hiny/hin - (V1 + V2) * einy/ein  \
                      - (1 - ein**2) * (5 * Seq * einy/ein - (4 * See - Sqq) * qiny/qin +  Sqh * hiny/hin))
    
    deinz_dt = ein * ((Z1 + Z2+ ZGR) * qinz/qin - (Y1 + Y2) * hinz/hin - (V1 + V2) * einz/ein  \
                      - (1 - ein**2) * (5 * Seq * einz/ein - (4 * See - Sqq) * qinz/qin +  Sqh * hinz/hin))
    

    dhinx_dt = hin * ((Y1 + Y2) * einx/ein - (X1 + X2) * qinx/qin - (W1 + W2) * hinx/hin \
                      + (1 - ein**2) * Sqh * einx/ein - (4 * ein**2 + 1) * Seh * qinx/qin \
                      +5 * ein**2 * Seq * hinx/hin)    
    
    dhiny_dt = hin * ((Y1 + Y2) * einy/ein - (X1 + X2) * qiny/qin - (W1 + W2) * hiny/hin \
                      + (1 - ein**2) * Sqh * einy/ein - (4 * ein**2 + 1) * Seh * qiny/qin \
                      +5 * ein**2 * Seq * hiny/hin)    

    dhinz_dt = hin * ((Y1 + Y2) * einz/ein - (X1 + X2) * qinz/qin - (W1 + W2) * hinz/hin \
                      + (1 - ein**2) * Sqh * einz/ein - (4 * ein**2 + 1) * Seh * qinz/qin \
                      +5 * ein**2 * Seq * hinz/hin) 
    
    
    dOmega1x_dt = mu * hin / I1 * (- Y1 * einx/ein + X1 * qinx/qin + W1 * hinx/hin)

    dOmega1y_dt = mu * hin / I1 * (- Y1 * einy/ein + X1 * qiny/qin + W1 * hiny/hin)

    dOmega1z_dt = mu * hin / I1 * (- Y1 * einz/ein + X1 * qinz/qin + W1 * hinz/hin)

    dOmega2x_dt = mu * hin / I2 * (- Y2 * einx/ein + X2 * qinx/qin + W2 * hinx/hin)

    dOmega2y_dt = mu * hin / I2 * (- Y2 * einy/ein + X2 * qiny/qin + W2 * hiny/hin)

    dOmega2z_dt = mu * hin / I2 * (- Y2 * einz/ein + X2 * qinz/qin + W2 * hinz/hin)


    dhoutx_dt = 0

    dhouty_dt = 0

    dhoutz_dt = 0


    return [deinx_dt, deiny_dt, deinz_dt, dhinx_dt, dhiny_dt, dhinz_dt, dOmega1x_dt, dOmega1y_dt, dOmega1z_dt, dOmega2x_dt, dOmega2y_dt, dOmega2z_dt,dhoutx_dt,dhouty_dt,dhoutz_dt]




#initial conditions
a_init = 5.0
e_init = 0.1
omega_init = 45.0 * np.pi/180.0
I_init = 85.6 * np.pi/180.0
Node_init  = 0.0 * np.pi/180.0

einx_init = e_init * (np.cos(Node_init)*np.cos(omega_init) - np.sin(Node_init) * np.sin(omega_init) * np.cos(I_init))
einy_init = e_init * (np.sin(Node_init)*np.cos(omega_init) + np.cos(Node_init) * np.sin(omega_init) * np.cos(I_init))
einz_init = e_init * np.sin(omega_init) * np.sin(I_init)

hin_init= np.sqrt(a_init * G * (m1 + m2) * (1 - e_init**2))
hinx_init =  hin_init * np.sin(Node_init) * np.sin(I_init)
hiny_init = -hin_init * np.cos(Node_init) * np.sin(I_init)
hinz_init =  hin_init * np.cos(I_init)
 
Omega1x_init = 0.0
Omega1y_init = 0.0
Omega1z_init = 2*np.pi/20.0
Omega2x_init = 0.0
Omega2y_init = 0.0
Omega2z_init = 2*np.pi/0.4

hout = np.sqrt(a_out * G * (m1 + m2 + m3) * (1 - e_out**2))
houtx_init = 0.0
houty_init = 0.0
houtz_init = hout

print hout, hinz_init

ics = [einx_init, einy_init, einz_init, hinx_init, hiny_init, hinz_init, Omega1x_init, Omega1y_init, Omega1z_init, Omega2x_init, Omega2y_init, Omega2z_init, houtx_init, houty_init, houtz_init]

#solve ODE
t = np.linspace(0,3.65e13,500)

sol=integ.odeint(secular_ode,ics,t)#,mxords=15,hmin=1e-4,hmax=100)

einx_sol = sol[:,0]
einy_sol = sol[:,1]
einz_sol = sol[:,2]
hinx_sol = sol[:,3]
hiny_sol = sol[:,4]
hinz_sol = sol[:,5]
Omega1x_sol = sol[:,6]
Omega1y_sol = sol[:,7]
Omega1z_sol = sol[:,8]
Omega2x_sol = sol[:,9]
Omega2y_sol = sol[:,10]
Omega2z_sol = sol[:,11]
houtx_sol = sol[:,12]
houty_sol = sol[:,13]
houtz_sol = sol[:,14]



ein_sol = np.sqrt(einx_sol**2 + einy_sol**2 + einz_sol**2)
hin_sol = np.sqrt(hinx_sol**2 + hiny_sol**2 + hinz_sol**2)
ain_sol = hin_sol**2 / G /(m1 + m2) / (1 - ein_sol**2)
I_sol = np.arccos(hinz_sol/hin_sol)*180.0/np.pi

#plt.plot(t,ain_sol)
#plt.axis([0,t.max(),0.0,8.0])
#plt.plot(t,I_sol)
#plt.axis([0,t.max(),0.0,180.0])
#plt.plot(t,ein_sol)
#plt.axis([0,t.max(),0.0,1.0])
#plt.plot(t,2*np.pi/Omega2z)
plt.plot(t,houtz_sol)
plt.show()
