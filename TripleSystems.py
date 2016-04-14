import numpy as np
import scipy.integrate as integ
import mercury

k2 = 39.4751488

def quadrupole_ode(y,t,m0,m1,m2,a1,a2,Gtot):

    ####################################################
    # twelve time-dependent variables

    G1 = y[0]
    G2 = y[1]
    H1 = y[2]
    H2 = y[3]

    g1 = y[4]
    g2 = y[5]
    h1 = y[6]
    h2 = y[7]
    e1 = y[8]
    e2 = y[9]

    #CosI1 = y[10]
    #CosI2 = y[11]
    CosI1 = (Gtot**2 + G1**2 - G2**2)/(2 * Gtot * G1)
    CosI2 = (Gtot**2 + G2**2 - G1**2)/(2 * Gtot * G2)
    #CosI1 = H1/G1
    #CosI2 = H2/G2
    ######################################################
    # some combination of variables for convenience
 
    CosItot = (Gtot**2 - G1**2 - G2**2)/(2 * G1 * G2)
    CosItot2 = CosItot**2
    SinItot2 = 1.0 - CosItot2
    SinItot = np.sqrt(SinItot2)
    Sin2Itot = 2 * SinItot * CosItot

    SinI1 = np.sqrt(1.0 - CosI1**2)
    SinI2 = np.sqrt(max(1.0 - CosI2**2,0.0))

    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

    C_2 = k2**2 /16 * (m0 + m1)**7 * m2**7 * L1**4 / (m0 + m1 + m2)**3 / m0**3 / m1**3 /L2**3 / G2**3

    ###################################################
    # twelve coupled differential equations

    dG1_dt = -30 * C_2 * e1**2 * SinItot2 * np.sin(2*g1)

    dG2_dt = 0

    #dH1_dt = -30 * C_2 * e1**2 * SinI2 * SinItot * np.sin(2*g1)
    
    dH1_dt = SinI2 / SinItot * dG1_dt - SinI1 / SinItot * dG2_dt 
 
    dH2_dt = -dH1_dt

    dg1_dt = 6 *C_2 *(1/G1 * (4 * CosItot2 + (5 * np.cos(2*g1) - 1)
                                * (1 - e1**2 - CosItot2))
                      + CosItot / G2 * (2 +e1**2*(3 - 5* np.cos(2*g1))))

    dg2_dt = 3 *C_2 *(2 * CosItot/G1 * (2 + e1**2*(3 - 5 * np.cos(2*g1)))
                      + 1/G2 * (4 + 6 * e1**2 + (5* CosItot2 - 3) *
                                  (2 + e1**2 * (3 - 5* np.cos(2*g1)))))
    
    dh1_dt =  -3 * C_2 / G1 / SinI1 * Sin2Itot * (2 + 3 * e1**2 - 5 *e1**2 *np.cos(2*g1))

    dh2_dt =  -3 * C_2 / G2 / SinI2 * Sin2Itot * (2 + 3 * e1**2 - 5 *e1**2 *np.cos(2*g1))

    de1_dt = C_2 * (1 - e1**2)/G1 * 30 * e1 * SinItot2 * np.sin(2*g1)

    de2_dt = 0

    #dCosI1_dt = dH1_dt / G1 - dG1_dt/G1 * CosI1

    #dCosI2_dt = dH2_dt / G2 - dG2_dt/G2 * CosI2


    return [dG1_dt, dG2_dt,  dH1_dt,  dH2_dt, dg1_dt, dg2_dt, dh1_dt, dh2_dt, de1_dt, de2_dt]#, dCosI1_dt,dCosI2_dt]


######################################################################################################

def octupole_ode(y,t,m0,m1,m2,a1,a2,Gtot):

    ####################################################
    # ten time-dependent variables

    G1 = y[0]
    G2 = y[1]
    H1 = y[2]
    H2 = y[3]

    g1 = y[4]
    g2 = y[5]
    h1 = y[6]
    h2 = y[7]
    e1 = y[8]
    e2 = y[9]

    #CosI1 = y[10]
    #CosI2 = y[11]

    G1sq, G2sq = G1**2, G2**2
    G1_1, G2_1 = 1.0/G1, 1.0/G2

    CosI1 = (Gtot**2 + G1sq - G2sq)/(2 * Gtot) * G1_1
    CosI2 = (Gtot**2 + G2sq - G1sq)/(2 * Gtot) * G2_1
    CosItot = (Gtot**2 - G1sq - G2sq)/2 * G1_1 * G2_1

    #CosI1 = H1/G1
    #CosI2 = H2/G2


    ######################################################
    # some combination of variables for convenience
 
    if (np.abs(CosItot - 1.0) < 1.0e-6):
        CosItot = 1.0
        SinItot = 0.0
        SinItot2 = 0.0
        Sin2Itot = 0.0
    else:
        CosItot2 = CosItot**2
        SinItot2 = 1.0 - CosItot2
        SinItot = np.sqrt(SinItot2)
        Sin2Itot = 2 * SinItot * CosItot

    if (np.abs(CosI1 - 1.0) < 1.0e-5):
        CosI1 = 1.0
        SinI1 = 0.0
    else:
        SinI1 = np.sqrt(1.0 - CosI1**2)

    if (np.abs(CosI2 - 1.0) < 1.0e-5):
        CosI2 = 1.0
        SinI2 = 0.0
    else:
        SinI2 = np.sqrt(1.0 - CosI2**2)


    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

    C_2 = k2**2 /16 * (m0 + m1)**7 * m2**7 * L1**4 \
        / (m0 + m1 + m2)**3 / m0**3 / m1**3 /L2**3 / G2**3

    C_3 = -15.0/16 * k2**2 /4 * (m0 + m1)**9 * m2**9 * (m0 - m1) * L1**6 \
        / (m0 + m1 + m2)**4 / m0**5 / m1**5 /L2**3 / G2**5 

    B = 2 + 5 * e1**2 - 7 * e1**2 * np.cos(2*g1)

    A = 4 + 3 * e1**2 - 2.5 * B * SinItot2

    Cosphi = -np.cos(g1) * np.cos(g2) - CosItot * np.sin(g1) * np.sin(g2)

    ###################################################
    # ten coupled differential equations

    dG1_dt = -30.0 * C_2 * e1**2 * SinItot2 * np.sin(2*g1) + \
        C_3 * e1 * e2 * (-35 * e1**2 * SinItot2 * np.sin(2*g1) * Cosphi
                          + A * (np.sin(g1) * np.cos(g2) - CosItot * np.cos(g1) * np.sin(g2))
                          + 10 * CosItot * SinItot2 * (1 - e1**2) * np.cos(g1) * np.sin(g2))

    dG2_dt = C_3 * e1 * e2 * (A * (np.cos(g1) *  np.sin(g2) - CosItot * np.sin(g1) * np.cos(g2)) \
                                  + 10 * CosItot * SinItot2 * (1 - e1**2) * np.sin(g1) * np.cos(g2))

    dH1_dt = SinI2 / SinItot * dG1_dt - SinI1 / SinItot * dG2_dt 

    dH2_dt = -dH1_dt

    dg1_dt = 6 *C_2 *(1.0/G1 * (4* CosItot2 + (5.0 * np.cos(2*g1) - 1.0) * (1.0 - e1**2 - CosItot2)) \
                      + CosItot / G2 * (2.0 +e1**2*(3.0 - 5.0* np.cos(2*g1)))) \
                      -C_3 * e2 *(e1 * (1.0/G2 + CosItot/G1) * (np.sin(g1) * np.sin(g2) *( \
                10 *(3*CosItot2 - 1)*(1.0 - e1**2) + A) - 5 * B * CosItot * Cosphi) \
                                      - (1.0 - e1**2) / e1 / G1 * (np.sin(g1) * np.sin(g2) * 10 * \
                                                                     CosItot * SinItot2 * (1.0 - 3 * e1**2)\
                                                                       + Cosphi * (3 * A - 10 * CosItot2 + 2)))

    
    dg2_dt = 3 *C_2 *(2* CosItot/G1 * (2.0 + e1**2*(3.0 - 5.0 * np.cos(2*g1))) \
                          + 1.0/G2 * (4.0 + 6.0 * e1**2 + (5.0* CosItot2 - 3.0) * \
                                          (2.0 + e1**2 * (3.0 - 5.0* np.cos(2*g1)))) \
                          ) \
                          + C_3 * e1 * (np.sin(g1) * np.sin(g2) *((4*e2**2 + 1) / e2 / G2 * 10 * CosItot * \
                                                                      SinItot2 * (1.0 - e1**2) \
                                                                      - e2 * (1.0 / G1 +  CosItot / G2) * \
                                                                      (A + 10 *(3 * CosItot2 -1)*(1.0 - e1**2))) \
                                            + Cosphi * (5 * B * CosItot * e2 * (1.0 / G1 + CosItot / G2) + \
                                                            (4 * e2**2 + 1) / e2 / G2 * A))
    
    dh1_dt =  -3.0 * C_2 / G1 / SinI1 * Sin2Itot * (2.0 + 3 * e1**2 - 5 *e1**2 *np.cos(2*g1)) \
        - C_3 * e1 * e2 * (5 * B * CosItot * Cosphi - A * np.sin(g1) * np.sin(g2) + \
                               10 * (1.0 - 3 * CosItot2) * (1.0 - e1**2) * np.sin(g1) * np.sin(g2)) * \
                               SinItot / G1 / SinI1
    
    dh2_dt =  dh1_dt
    
    de1_dt = C_2 * (1.0 - e1**2)/G1 * 30 * e1 * SinItot2 * np.sin(2*g1) +\
        C_3 * e2 * (1.0 - e1**2) / G1 * (35 * Cosphi * SinItot2 * e1**2 * np.sin(2*g1)\
                                           -10 *CosItot * SinItot2 * np.cos(g1) * np.sin(g2) * (1.0 - e1**2) \
                                           -A *(np.sin(g1) * np.cos(g2) - CosItot * np.cos(g1) * np.sin(g2)))
        
    de2_dt = -C_3 * e1 *(1.0 - e2**2) / G2 * (10 * CosItot * SinItot2 * (1.0 - e1**2) * np.sin(g1)* np.cos(g2)\
                                                + A * (np.cos(g1) * np.sin(g2) - \
                                                           CosItot * np.sin(g1) * np.cos(g2)))

    #dCosI1_dt = dH1_dt / G1 - dG1_dt/G1 * CosI1

    #dCosI2_dt = dH2_dt / G2 - dG2_dt/G2 * CosI2


    return [dG1_dt, dG2_dt,  dH1_dt,  dH2_dt, dg1_dt, dg2_dt, dh1_dt, dh2_dt, de1_dt, de2_dt]#, dCosI1_dt,dCosI2_dt]


def octupole_ode_alternative(y,t,m0,m1,m2,a1,a2,Gtot):

    ####################################################
    # eight time-dependent variables
    # only orbital elements in this case, not canonical variables
    
    e1 = y[0]
    e2 = y[1]
    I1 = y[2]
    I2 = y[3]
    g1 = y[4]
    g2 = y[5]
    h1 = y[6]
    h2 = y[7]

    

    CosI1 = np.cos(I1)
    CosI2 = np.cos(I2)
    CosItot = np.cos(I1 + I2)


    ######################################################
    #Constants in the system
    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)
    
    # some combination of variables for convenience
    CosItot2 = CosItot**2
    SinItot2 = 1.0 - CosItot2
    SinItot = np.sqrt(SinItot2)
    Sin2Itot = 2 * SinItot * CosItot
    
    SinI1 = np.sqrt(1.0 - CosI1**2)
    SinI2 = np.sqrt(max(1.0 - CosI2**2,0.0))
    
    G1 = L1 * np.sqrt(1.0 - e1**2)
    G2 = L2 * np.sqrt(1.0 - e2**2)
    
    C_2 = k2**2 /16 * (m0 + m1)**7 * m2**7 * L1**4 \
        / (m0 + m1 + m2)**3 / m0**3 / m1**3 /L2**3 / G2**3
    
    C_3 = -15.0/16 * k2**2 /4 * (m0 + m1)**9 * m2**9 * (m0 - m1) * L1**6 \
        / (m0 + m1 + m2)**4 / m0**5 / m1**5 /L2**3 / G2**5 
    
    B = 2 + 5 * e1**2 - 7 * e1**2 * np.cos(2*g1)
    
    A = 4 + 3 * e1**2 - 2.5 * B * SinItot2
    
    Cosphi = -np.cos(g1) * np.cos(g2) - CosItot * np.sin(g1) * np.sin(g2)
    
    # in addition to some derivatives that will not be solved directly...
    dG1_dt = -30.0 * C_2 * e1**2 * SinItot2 * np.sin(2*g1) + \
        C_3 * e1 * e2 * (-35 * e1**2 * SinItot2 * np.sin(2*g1) * Cosphi
                          + A * (np.sin(g1) * np.cos(g2) - CosItot * np.cos(g1) * np.sin(g2))
                          + 10 * CosItot * SinItot2 * (1 - e1**2) * np.cos(g1) * np.sin(g2))
        
    dG2_dt = C_3 * e1 * e2 * (A * (np.cos(g1) *  np.sin(g2) - CosItot * np.sin(g1) * np.cos(g2)) \
                                  + 10 * CosItot * SinItot2 * (1 - e1**2) * np.sin(g1) * np.cos(g2))
    
    dH1_dt = SinI2 / SinItot * dG1_dt - SinI1 / SinItot * dG2_dt 
    
    dH2_dt = -dH1_dt

    ####################################################################################################
    # eight coupled differential equations
    #eccentricities
    de1_dt = C_2 * (1.0 - e1**2)/G1 * 30 * e1 * SinItot2 * np.sin(2*g1) +\
                          C_3 * e2 * (1.0 - e1**2) / G1 * (35 * Cosphi * SinItot2 * e1**2 * np.sin(2*g1)\
                                             -10 *CosItot * SinItot2 * np.cos(g1) * np.sin(g2) * (1.0 - e1**2) \
                                           -A *(np.sin(g1) * np.cos(g2) - CosItot * np.cos(g1) * np.sin(g2)))
        
    de2_dt = -C_3 * e1 *(1.0 - e2**2) / G2 * (10 * CosItot * SinItot2 * (1.0 - e1**2) * np.sin(g1)* np.cos(g2)\
                                                + A * (np.cos(g1) * np.sin(g2) - \
                                                           CosItot * np.sin(g1) * np.cos(g2)))

    #inclinations
    dI1_dt = -(dH1_dt / G1 - dG1_dt/G1 * CosI1)  / SinI1

    dI2_dt = -(dH2_dt / G2 - dG2_dt/G2 * CosI2) / SinI2
    

    #dCosI1_dt = dH1_dt / G1 - dG1_dt/G1 * CosI1
    
    #dCosI2_dt = dH2_dt / G2 - dG2_dt/G2 * CosI2

    # argument of pericenter
    dg1_dt = 6 *C_2 *(1.0/G1 * (4* CosItot2 + (5.0 * np.cos(2*g1) - 1.0) * (1.0 - e1**2 - CosItot2)) \
                      + CosItot / G2 * (2.0 +e1**2*(3.0 - 5.0* np.cos(2*g1)))) \
                      -C_3 * e2 *(e1 * (1.0/G2 + CosItot/G1) * (np.sin(g1) * np.sin(g2) *( \
                10 *(3*CosItot2 - 1)*(1.0 - e1**2) + A) - 5 * B * CosItot * Cosphi) \
                                      - (1.0 - e1**2) / e1 / G1 * (np.sin(g1) * np.sin(g2) * 10 * \
                                                                     CosItot * SinItot2 * (1.0 - 3 * e1**2)\
                                                                       + Cosphi * (3 * A - 10 * CosItot2 + 2)))

    
    dg2_dt = 3 *C_2 *(2* CosItot/G1 * (2.0 + e1**2*(3.0 - 5.0 * np.cos(2*g1))) \
                          + 1.0/G2 * (4.0 + 6.0 * e1**2 + (5.0* CosItot2 - 3.0) * \
                                          (2.0 + e1**2 * (3.0 - 5.0* np.cos(2*g1)))) \
                          ) \
                          + C_3 * e1 * (np.sin(g1) * np.sin(g2) *((4*e2**2 + 1) / e2 / G2 * 10 * CosItot * \
                                                                      SinItot2 * (1.0 - e1**2) \
                                                                      - e2 * (1.0 / G1 +  CosItot / G2) * \
                                                                      (A + 10 *(3 * CosItot2 -1)*(1.0 - e1**2))) \
                                            + Cosphi * (5 * B * CosItot * e2 * (1.0 / G1 + CosItot / G2) + \
                                                            (4 * e2**2 + 1) / e2 / G2 * A))
    #longitude of nodes

    dh1_dt =  -3.0 * C_2 / G1 / SinI1 * Sin2Itot * (2.0 + 3 * e1**2 - 5 *e1**2 *np.cos(2*g1)) \
        - C_3 * e1 * e2 * (5 * B * CosItot * Cosphi - A * np.sin(g1) * np.sin(g2) + \
                               10 * (1.0 - 3 * CosItot2) * (1.0 - e1**2) * np.sin(g1) * np.sin(g2)) * \
                               SinItot / G1 / SinI1
    
    dh2_dt =  dh1_dt


    return [de1_dt, de2_dt,  dI1_dt,  dI2_dt, dg1_dt, dg2_dt, dh1_dt, dh2_dt]

######################################################################################################

def octupole_ode_bin(y,t,m0,m1,m2,a1,a2,Gtot):

    ####################################################
    # seven time-dependent variables

    e1 = y[0]
    e2 = y[1]
    I1 = y[2]
    I2 = y[3]
    g1 = y[4]
    g2 = y[5]
    h1 = y[6]


    ######################################################
    # some combination of variables for convenience
    CosI1 = np.cos(I1)
    CosI2 = np.cos(I2)
    CosItot = np.cos(I1 + I2)

 
    if (np.abs(CosItot - 1.0) < 1.0e-6):
        CosItot = 1.0
        SinItot = 0.0
        SinItot2 = 0.0
        Sin2Itot = 0.0
        Cos2Itot = 1.0
    else:
        CosItot2 = CosItot**2
        SinItot2 = 1.0 - CosItot2
        SinItot = np.sqrt(SinItot2)
        Sin2Itot = 2 * SinItot * CosItot
        Cos2Itot = CosItot2 - SinItot2

    if (np.abs(CosI1 - 1.0) < 1.0e-5):
        CosI1 = 1.0
        SinI1 = 0.0
    else:
        SinI1 = np.sqrt(1.0 - CosI1**2)

    if (np.abs(CosI2 - 1.0) < 1.0e-5):
        CosI2 = 1.0
        SinI2 = 0.0
    else:
        SinI2 = np.sqrt(1.0 - CosI2**2)


    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

    A = k2 * m0 * m1 * m2 * a1**2 \
        / 8 / a2**3 / (m0 + m1)

    B = 15.0 / 16 * k2 * m0 * m1 * m2 * (m0 - m1) * a1**3 \
        / 4 / a2**4 / (m0 + m1)**2


    ###################################################
    # seven coupled differential equations

    de1_dt = np.sqrt(1-e1**2) / 8 / (1-e2**2)**2.5 / L1 * \
        (-120 * A * e1 * (-1 + e2**2) * SinItot2 * np.sin(2*g1) + B * e2 * np.cos(g2) * \
              ((4 + 3 * e1**2)*(3 + 5* Cos2Itot) * np.sin(g1) + 210 * e1**2 * SinItot2 * np.sin(3*g1)) - \
              2 * B * e2 * CosItot * np.cos(g1) * (15* (2 + 5 * e1**2)* Cos2Itot + \
                                                       7* (-2 -9 * e1**2 + 30* e1**2 * np.cos(2*g1)\
                                                                * SinItot2 ))* np.sin(g2))
    
    de2_dt = B * e1 /8 / (-1 + e2**2)**2 /L2 * \
        (-2 * CosItot * np.cos(g2) * (5* (6 + e1**2) * Cos2Itot + 7 \
                                          * (-2 + e1**2 + 10 * e1**2 *np.cos(2*g1) * SinItot2))\
              * np.sin(g1) + 2 * np.cos(g1) *(6 - 13* e1**2 + 5*(2 + 5* e1**2)* Cos2Itot \
                                                  + 70 * e1**2 * np.cos(2*g1) * SinItot2 )* np.sin(g2) ) 
             
    dI1_dt = e1 / 8 / np.sqrt(1-e1**2)/ (1-e2**2)**2.5 / L1 * ( 20 * Sin2Itot * (-B * e2\
             * (2 + 5 * e1**2 + 7 * e1**2 * np.cos(2*g1)) *  np.cos(g2) * np.sin(g1) + 3 * A * e1 * (-1 +e2**2)\
             * np.sin(2 * g1)) + 2 * B * e2 * np.cos(g1)* (-26 - 37* e1**2 + 35* e1**2 * np.cos(2*g1) + 15\
             * Cos2Itot * (-2 - 5 * e1**2 + 7 * e1**2 *np.cos(2*g1))) * SinItot * np.sin(g2))

    dI2_dt = e1 /8 / (-1 + e2**2)**2 /L2 * (2 * B * e2 * (26 + 107 * e1**2 + 5* (6 + e1**2)* Cos2Itot\
             -35* e1**2 *(-5 + Cos2Itot)* np.cos(2*g1)) * np.cos(g2) * SinItot * np.sin(g1)\
             + 20* (-6 * A * e1 * (-1 + e2**2) * SinItot * np.sin(2*g1) + B * e2 * np.cos(g1) * (2\
             + 5* e1**2 - 7 * e1**2 * np.cos(2*g1)) * Sin2Itot * np.sin(g2))) 

    dg1_dt = 1.0 / 32 / e1 / np.sqrt(1-e1**2)/ (1-e2**2)**2.5 / L1 / SinI1 *(24 * A * e1 * (-1 + e2**2)\
             * (2 *(-1 + e2**2) * SinI1 + 10 * ((-1 + 2* e1**2) * CosI2 + np.cos(2*I1 + I2))\
             * np.cos(2*g1) * SinItot + np.sin(I1 + 2*I2) - 6*e1**2*np.sin(I1 + 2*I2) \
             - 5 * np.sin(3*I1 + 2*I2)) + 4 * B * e2 * np.cos(g1) * np.cos(g2) * (2 * (6 - 45*e1**2\
             + 39* e1**4 + 10* Cos2Itot) * SinI1 + 5 * e1**2 * (14*((3-5*e1**2)* CosI2\
             + (-3 + e1**2) * np.cos(2*I1 + I2))* np.cos(2*g1) * SinItot + (-9 + 25*e1**2)\
             * np.sin(I1 +2*I2) + (17 - 5 * e1**2)* np.sin(3*I1 +2*I2))) + 2 * B * e2 *(60 *np.cos(3*(I1+I2))\
             * SinI1 + (-2 - 53 * e1**2 + 76* e1**4)* SinI2 + 2*(-70 * e1**4 * (CosI1\
             + 3* np.cos(I1 + 2*I2))*np.cos(2*g1)* SinItot + np.sin(2*I1 + I2))+ e1**2\
             * (840* CosItot * np.cos(2*g1) * SinI1 * SinItot2 + (57 - 38* e1**2)\
             * np.sin(2*I1 + I2) + 15 * (7 + 2*e1**2) * np.sin(2*I1 + 3*I2) + 75* np.sin(4*I1 + 3*I2)))\
             * np.sin(g1) * np.sin(g2))  
    
    dg2_dt = 1.0 / 32 / e2 / (-1 + e2**2)**3 / L2 / SinI2 * (8* B * e1 * np.cos(g1) * np.cos(g2)\
             * (5 * CosI2 *(-1 - 4*e2**2 + (1 + 6 * e2**2)* np.cos(2*I2)) * (-2 - 5* e1**2 + 7 \
             * e1**2 * np.cos(2*g1)) * np.sin(2*I1) + (5 * np.cos(2*I1)*(2*e2**2 + (1 + 6* e2**2)\
             * np.cos(2*I2))*(-2 -5*e1**2 + 7* e1**2 * np.cos(2*g1))- (1+ 4*e2**2)*(6-13*e1**2 + 35*e1**2\
             * np.cos(2*g1)))*SinI2) + 12* A * e2 *(-1 + e2**2)*((4 + 6* e1**2 + 30* e1**2 * np.cos(2*g1))\
             * SinI2 + (-2 -3*e1**2 + 5* e1**2 * np.cos(2*g1))*(np.sin(2*I1 + I2)- 5 * np.sin(2*I1 + 3*I2)))\
             + 2*B*e1*((2 + 19*e1**2)*(1 + 3*e2**2)*SinI1 + 70*e1**2*(2* e2**2*CosI2-(1+e2**2)\
             * np.cos(2*I1 + I2)+ (1+7*e2**2)*np.cos(2*I1 + 3*I2))* np.cos(2*g1) * SinItot\
             - 2*np.sin(I1 + 2*I2) - 19 * e1**2*np.sin(I1 + 2*I2) - 10*e2**2 * np.sin(I1 + 2*I2)\
             -95* e1**2 * e2**2 *np.sin(I1 + 2*I2) + 30* np.sin(3*I1 + 2*I2) + 5*e1**2 * np.sin(3*I1 + 2*I2)\
             + 30 * e2**2 * np.sin(3*I1 + 2*I2) + 5* e1**2 * e2**2 * np.sin(3*I1 + 2*I2) - 30*np.sin(3*I1 + 4*I2)\
             - 5*e1**2 * np.sin(3*I1 + 4*I2) - 210 * e2**2 * np.sin(3*I1 + 4*I2) - 35* e1**2 * e2**2 \
             * np.sin(3*I1 + 4*I2))* np.sin(g1) * np.sin(g2))
    
    dh1_dt = 1.0 / 8 / np.sqrt(1-e1**2)/ (1-e2**2)**(2.5) / L1 / SinI1 * (-4 * (3 * A * (-1 +e2**2)\
             * (-2 - 3 * e1**2 + 5 * e1**2 * np.cos(2*g1)) + 5 * B * e1 * e2 * np.cos(g1)\
             * (2 + 5 * e1**2 - 7 * e1**2 * np.cos(2*g1))* np.cos(g2)) * Sin2Itot\
             + 2 * B * e1 * e2 * (-46 - 17* e1**2 - 15*(6+e1**2)* Cos2Itot + 35\
             * e1**2 * (1 + 3* Cos2Itot) * np.cos(2*g1)) * SinItot\
             * np.sin(g1) * np.sin(g2)) 


    return [de1_dt, de2_dt,  dI1_dt,  dI2_dt, dg1_dt, dg2_dt, dh1_dt]


############################################################################################################
############################################################################################################

def octupole_evolution(m0,m1,m2,a1,a2,e1,e2,Itot,omega1,omega2,Omega1,Omega2,time):
    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

  
    #find angles respect to invariable plane
    G1_0 = L1 * np.sqrt(1.0 - e1**2)
    G2_0 = L2 * np.sqrt(1.0 - e2**2)
    Gtot = np.sqrt(np.cos(Itot) * 2 * G1_0 * G2_0 + G1_0**2 + G2_0**2)
    CosI2 = (Gtot**2 + G2_0**2 - G1_0**2)/2/Gtot/G2_0
    CosI1 = (Gtot**2 + G1_0**2 - G2_0**2)/2/Gtot/G1_0
    I2 = np.arccos(CosI2)
    I1 = Itot - I2
    #I1 = np.arccos(CosI1)
    #print "Individual inclinations: i1=",I1/np.pi*180.0,"i2=",I2/np.pi*180.0

    H1_0 = G1_0 * CosI1
    H2_0 = G2_0 * CosI2
    g1_0 = omega1
    g2_0 = omega2
    h1_0 = Omega1
    h2_0 = Omega2
    e1_0 = e1
    e2_0 = e2

    ics = [G1_0, G2_0, H1_0, H2_0, g1_0, g2_0, h1_0, h2_0, e1_0, e2_0]
    
    #Gtot = H1_0 + H2_0

    sol = integ.odeint(octupole_ode,ics,time,args=(m0,m1,m2,a1,a2,Gtot,),rtol=10e-20)
    
    return sol

def octupole_evolution_alternative(m0,m1,m2,a1,a2,e1,e2,Itot,omega1,omega2,Omega1,Omega2,time):
    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

   
    # only the relative inclination is a free parameter, individual inclinations depend on the
    # angular momenta:
    #find angles respect to invariable plane
    G1_0 = L1 * np.sqrt(1.0 - e1**2)
    G2_0 = L2 * np.sqrt(1.0 - e2**2)
    Gtot = np.sqrt(np.cos(Itot) * 2 * G1_0 * G2_0 + G1_0**2 + G2_0**2)
    CosI2 = (Gtot**2 + G2_0**2 - G1_0**2)/2/Gtot/G2_0
    CosI1 = (Gtot**2 + G1_0**2 - G2_0**2)/2/Gtot/G1_0
    I2 = np.arccos(CosI2)
    I1 = Itot - I2
    #I1 = np.arccos(CosI1)
    #print "Individual inclinations: i1=",I1/np.pi*180.0,"i2=",I2/np.pi*180.0

    e1_0 = e1
    e2_0 = e2
    I1_0 = I1
    I2_0 = I2
    g1_0 = omega1
    g2_0 = omega2
    h1_0 = Omega1
    h2_0 = Omega2


    ics = [e1_0, e2_0, I1_0, I2_0, g1_0, g2_0, h1_0, h2_0]
    
    sol = integ.odeint(octupole_ode_alternative,ics,time,args=(m0,m1,m2,a1,a2,Gtot,),rtol=10e-20)
    
    return sol

def octupole_evolution_binliu(m0,m1,m2,a1,a2,e1,e2,Itot,omega1,omega2,Omega1,Omega2,time):
    L1 = m0 * m1 / (m0 + m1)* np.sqrt(k2*(m0 + m1) * a1)
    L2 = (m0 + m1)* m2  / (m0 + m1 + m2)* np.sqrt(k2*(m0 + m1 + m2) * a2)

   
    # only the relative inclination is a free parameter, individual inclinations depend on the
    # angular momenta:
    #find angles respect to invariable plane
    G1_0 = L1 * np.sqrt(1.0 - e1**2)
    G2_0 = L2 * np.sqrt(1.0 - e2**2)
    Gtot = np.sqrt(np.cos(Itot) * 2 * G1_0 * G2_0 + G1_0**2 + G2_0**2)
    CosI2 = (Gtot**2 + G2_0**2 - G1_0**2)/2/Gtot/G2_0
    CosI1 = (Gtot**2 + G1_0**2 - G2_0**2)/2/Gtot/G1_0
    I2 = np.arccos(CosI2)
    I1 = Itot - I2
    #I1 = np.arccos(CosI1)
    #print "Individual inclinations: i1=",I1/np.pi*180.0,"i2=",I2/np.pi*180.0

    e1_0 = e1
    e2_0 = e2
    I1_0 = I1
    I2_0 = I2
    g1_0 = omega1
    g2_0 = omega2
    h1_0 = Omega1


    ics = [e1_0, e2_0, I1_0, I2_0, g1_0, g2_0, h1_0]
    
    sol = integ.odeint(octupole_ode_bin,ics,time,args=(m0,m1,m2,a1,a2,Gtot,),rtol=10e-20)
    
    return sol

