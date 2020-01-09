import numpy as np
from . import stretching
def set_depth(Vtransform, Vstretching,theta_s, theta_b, hc, N,igrid, h, zeta, report):
    '''
	SET_DEPTH:  Compute ROMS grid depth from vertical stretched variables
	
	 [z]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
	               igrid, h, zeta);
	               
	 Given a batymetry (h), free-surface (zeta) and terrain-following parameters,
	 this function computes the 3D depths for the request C-grid location. If the
	 free-surface is not provided, a zero value is assumef resulting in unperturb
	 depths.  This function can be used when generating initial conditions or
	 climatology data for an application. Check the following link for details:
	 
	    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
	    
	 On Input:
	 
	    Vtransform    Vertical transformation equation:
	                    Vtransform = 1,   original transformation
	                    
	                      z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
	                      
	                      Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
	                      
	                    Vtransform = 2,   new transformation
	                    
	                      z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
	                      
	                       Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
	    Vstretching   Vertical stretching function:
	                    Vstretching = 1,  original (Song and Haidvogel, 1994)
	                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
	                    Vstretching = 3,  R. Geyer BBL refinement
	                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
	 theta_s       S-coordinate surface control parameter (scalar)
	 theta_b       S-coordinate bottom control parameter (scalar)
	 hc            Width (m) of surface or bottom boundary layer in which
	                    higher vertical resolution is required during
	                    stretching (scalar)
	  N             Number of vertical levels (scalar)
	  igrid         Staggered grid C-type (integer):
	                    igrid=1  => density points
	                    igrid=2  => streamfunction points
	                    igrid=3  => u-velocity points
	                    igrid=4  => v-velocity points
	                    igrid=5  => w-velocity points
	  h             Bottom depth, 2D array at RHO-points (m, positive),
	                    h(1:Lp+1,1:Mp+1).
	  zeta          Free-surface, 2D array at RHO-points (m), OPTIONAL,
                    zeta(1:Lp+1,1:Mp+1).
     report        Flag to report detailed information (OPTIONAL):
                     report = 0,       do not report
                     report = 1,       report information
 
  On Output:
 
     z             Depths (m, negative), 3D array.
 

  svn $Id: set_depth.m 586 2012-01-03 20:19:25Z arango $
 ===========================================================================#
   Copyright (c) 2002-2012 The ROMS/TOMS Group                              #
     Licensed under a MIT/X style license                                   #
     See License_ROMS.txt                           Hernan G. Arango        #
 ===========================================================================#
    '''
    # the line below is in the original matlab code but is useless for python
    #z=np.array([])
    # Note zeta needs to be only 2-D (squeeze and remove the time coordinate)
    # the array will be too big otherwise and cause problems.
    # lots of parameter checking code missing
    if(Vtransform < 1)or(Vtransform > 2):
            print('*** Error: Set_depth- Illegal parameter Vtransform='+str(Vtransform))
            return
    if(Vstretching < 1)or(Vstretching > 4):
            print('*** Error: Set_depth- Illegal parameter Vstretching='+str(Vstretching))
            return
    if(hc > np.min(h))and(Vtransform==1):
            print('*** Error: Set_depth- critical depth exceeds minimum bathymetry value.')
            print('Vtransform= '+str(Vtransform))
            print('hc= '+str(hc))
            print('hmax= '+str(np.min(hc)))
            return
    Np=N+1
    # Find the shape of the bathymetry data
    [Lp,Mp]=h.shape
    L=Lp-1
    M=Mp-1
    # get the max and min seafloor depths
    hmin=np.min(h)
    hmax=np.max(h)
    # set up a grid parameter 
    if (igrid==5):
        kgrid=1
    else:
        kgrid=0
    # compute stretching terms
    [s,C]=stretching.stretching(Vstretching,theta_s,theta_b,hc,N,kgrid,report)
    if (igrid==1):
        hr=h
        zetar=zeta
    elif(igrid==2):
        hp=0.25*(h[0:L,0:M]+h[1:Lp,0:M]+h[0:L,1:Mp]+h[1:Lp,1:Mp])
        zetap=0.25*(zeta[0:L,0:M]+zeta[1:Lp,0:M]+zeta[0:L,1:Mp]+zeta[1:Lp,1:Mp])
    elif(igrid==3):
        # I'm thinking that L and M are opposite of what they are in matlab
        #hu=0.5*(h[0:L,0:Mp]+h[1:Lp,0:Mp])
        #zetau=0.5*(zeta[0:L,0:Mp]+zeta[1:Lp,0:Mp])
        hu=0.5*(h[0:Lp,0:M]+h[0:Lp,1:Mp])
        zetau=0.5*(zeta[0:Lp,0:M]+zeta[0:Lp,1:Mp])
    elif(igrid==4):
        hv=0.5*(h[0:L,0:Mp]+h[1:Lp,0:Mp])
        zetav=0.5*(zeta[0:L,0:Mp]+zeta[1:Lp,0:Mp])
        #hv=0.5*(h[0:Lp,0:M]+h[0:Lp,1:Mp])
        #zetav=0.5*(zeta[0:Lp,0:M]+zeta[0:Lp,1:Mp])
    elif(igrid==5):
        hr=h
        zetar=zeta
    # Now actually compute the grid of depths.  This may not be efficient.
    if(Vtransform==1):
        if(igrid==1):
            for k in np.arange(0,N):
                z0=(s[k]-C[k])*hc+C[k]*hr
                zlev=z0+zetar*(1.0+z0/hr)
                if k==0:
                    z=zlev
                else:
                    z=np.dstack((z,zlev))
        elif(igrid==2):
            for k in np.arange(0,N):
                z0=(s[k]-C[k])*hc+C[k]*hp
                zlev=z0+zetap*(1.0+z0/hp)
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
        elif(igrid==3):
            for k in np.arange(0,N):
                z0=(s[k]-C[k])*hc+C[k]*hu
                zlev=z0+zetau*(1.0+z0/hu)
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
        elif(igrid==4):
            for k in np.arange(0,N):
                z0=(s[k]-C[k])*hc+C[k]*hv
                zlev=z0+zetav*(1.0+z0/hv)
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
        elif(igrid==5):
            z=-hr
            for k in np.arange(1,Np):
                z0=(s[k]-C[k])*hc+C[k]*hr
                zlev=z0+zetar*(1.0+z0/hr)
                z=np.dstack((z,zlev))
    elif(Vtransform==2):
        if(igrid==1):
            for k in np.arange(0,N):
                z0=(hc*s[k]+C[k]*hr)/(hc+hr)
                zlev=zetar+(zeta+hr)*z0
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
                #z[:,:,k]=zetar+(zeta+hr)*z0
        elif(igrid==2):
            for k in np.arange(0,N):
                z0=(hc*s[k]+C[k]*hp)/(hc+hp)
                zlev=zetap+(zetap+hp)*z0
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
                #z[:,:,k]=zetap+(zetap+hp)*z0
        elif(igrid==3):
            for k in np.arange(0,N):
                z0=(hc*s[k]+C[k]*hu)/(hc+hu)
                zlev=zetau+(zetau+hu)*z0
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
                #z[:,:,k]=zetau+(zetau+hu)*z0
        elif(igrid==4):
            for k in np.arange(0,N):
                z0=(hc*s[k]+C[k]*hv)/(hc+hv)
                zlev=zetav+(zetav+hv)*z0
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
                #z[:,:,k]=zetav+(zetav+hv)*z0
        elif(igrid==5):
            for k in np.arange(0,Np):
                z0=(hc*s[k]+C[k]*hr)/(hc+hr)
                zlev=zetar+(zetar+hr)*z0
                if k==0:
                	z=zlev
                else:
                	z=np.dstack((z,zlev))
                #z[:,:,k]=zetar+(zetar+hr)*z0
    return z
