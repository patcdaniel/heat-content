import numpy as np
#import numba
#@numba.jit(nopython=True)
def stretching(Vstretching,theta_s,theta_b,hc,N,kgrid, report):
    '''
#function [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, ...
#                          report)
#
# STRETCHING:  Compute ROMS vertical coordinate stretching function
#
# [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)
#
# Given vertical terrain-following vertical stretching parameters, this
# this routine computes the vertical stretching function used used in
# ROMS vertical coordinate transformation. Check the following link for
# details:
#
#    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
#
# On Input:
#
#    Vstretching   Vertical stretching function:
#                    Vstretching = 1,  original (Song and Haidvogel, 1994)
#                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
#                    Vstretching = 3,  R. Geyer BBL refinement
#                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
#    theta_s       S-coordinate surface control parameter (scalar)
#    theta_b       S-coordinate bottom control parameter (scalar)
#    hc            Width (m) of surface or bottom boundary layer in which
#                    higher vertical resolution is required during
#                    stretching (scalar)
#    N             Number of vertical levels (scalar)
#    kgrid         Depth grid type logical switch:
#                    kgrid = 0,        function at vertical RHO-points
#                    kgrid = 1,        function at vertical W-points
#    report        Flag to report detailed information (OPTIONAL):
#                    report = false,   do not report
#                    report = true,    report information
#
# On Output:
#
#    s             S-coordinate independent variable, [-1 <= s <= 0] at
#                    vertical RHO- or W-points (vector)
#    C             Nondimensional, monotonic, vertical stretching function,
#                    C(s), 1D array, [-1 <= C(s) <= 0]
#

# svn $Id: stretching.m 595 2012-01-17 19:44:46Z arango $
#=========================================================================#
#  Copyright (c) 2002-2012 The ROMS/TOMS Group                            #
#    Licensed under a MIT/X style license                                 #
#    See License_ROMS.txt                           Hernan G. Arango      #
#=========================================================================#
    '''
    if (Vstretching < 1)or(Vstretching > 4):
        print('*** Error Stretching-Illegal parameter Vstretching='+str(Vstretching))
        return
    # this code is unnecessary
    #s=np.array([])
    #C=np.array([])
    # need code for checking number of arguments
    Np=N+1
    if(Vstretching==1):
        ds=1.0/N
        if (kgrid==1):
            Nlev=Np
            lev=np.arange(0,N+1)
            s=(lev-N)*ds
        else:
            Nlev=N
            lev=(np.arange(1,N+1))-0.5
            s=(lev-N)*ds
        if (theta_s > 0):
            Ptheta=np.sinh(theta_s*s)/np.sinh(theta_s)
            Rtheta=np.tanh(theta_s*(s+0.5))/(2*np.tanh(0.5*theta_s))-0.5
            C=(1.0-theta_b)*Ptheta+theta_b*Rtheta
        else:
            C=s
    elif (Vstretching==2):
        alfa=1.0
        beta=1.0
        ds=1.0/N
        if(kgrid==1):
            Nlev=Np
            lev=np.arange(0,N+1)
            s=(lev-N)*ds
        else:
            Nlev=N
            lev=np.arange(1,N+1)-0.5
            s=(lev-N)*ds
        if (theta_s > 0):
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if (theta_b > 0):
                Cbot=-1.0+np.sinh(theta_b*(s+1.0))/np.sinh(theta_b)
                weigth=np.power((s+1.0),alfa)*(1+(alfa/bet)*(1.0-np.power((s+1.0),beta)))
                C=weigth*Csur+(1-weigth)*Cbot
            else:
                C=Csur
        else:
            C=s
    # R. Geyer BBL vertical stretching function
    elif(Vstretching==3):
        ds=1.0/N
        if(kgrid==1):
            Nlev=Np
            lev=np.arange(0,N+1)
            s=(lev-N)*ds
        else:
            Nlev=N
            lev=np.arange(1,N+1)-0.5
            s=(lev-N)*ds
        if(theta_s > 0):
            exp_s=theta_s
            exp_b=theta_b
            alpha=3
            Cbot=np.log(np.cosh(alpha*np.power((s+1),exp_b)))/np.log(np.cosh(alpha))-1
            Csur=-np.log(np.cosh(alpha*np.power(np.abs(s),exp_s)))/np.log(np.cosh(alpha))
            weight=(1-np.tanh(alpha*(s+0.5)))/2
            C=weight*Cbot+(1-weight)*Csur
        else:
            C=s
    # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
    # with bottom refinement
    elif(Vstretching==4):
        ds=1.0/N
        if(kgrid==1):
            Nlev=Np
            lev=np.arange(0,N+1)
            s=(lev-N)*ds
        else:
            Nlev=N
            lev=np.arange(1,N+1)-0.5
            s=(lev-N)*ds
        if (theta_s > 0):
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=np.power(-s,2)
        if (theta_b > 0):
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
    if (report):
        print(' ')
        if(Vstretching==1):
            print('Vstretching='+str(Vstretching)+' Song and Haidvogel (1994)')
        elif(Vstretching==2):
            print('Vstretching='+str(Vstretching)+' Shchepetking (2005)')
        elif(Vstretching==3):
            print('Vstretching='+str(Vstretching)+' Geyer (2009), BBL')
        elif(Vstretching==4):
            print('Vstretching='+str(Vstretching)+' Shchepetkin (2010)')
        if(kgrid==1):
            print('kgrid = '+str(kgrid)+' at vertical W-points')
        else:
            print('kgrid = '+str(kgrid)+' at vertical RHO-points')
        print('theta_s= '+str(theta_s))
        print('theta_b= '+str(theta_b))
        print('hc= '+str(hc))
        print(' S-coordinate curves: k, s(k), C(k)')
        if (kgrid==1):
            for k in np.arange(Nlev,1,-1):
                print(' ')
                print(k-1,s[k],C[k])
        else:
            for k in np.arange(Nlev,1,-1):
                print(' ')
                print(k,s[k],C[k])
    return [s,C]
