# -*- coding: utf-8 -*-
# lesson 11 - Source-vortex panel method
# started 8 Apr 2014
# Ian Carr

# importing libraries
import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# reading geometry from file
coords = np.loadtxt(fname='C:/Users/Ian/Documents/GitHub/AeroHydro/resources/naca0012.dat')
xp,yp = coords[:,0],coords[:,1]

# --------- creating and applying panels ----------

# defining class Panel with info on one panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya         # defining 1st point
        self.xb,self.yb = xb,yb         # defining 2nd point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2   # defining center point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2) # length of the panel
        
        # orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados' 
        
        self.sigma = 0          # source strength
        self.vt = 0             # tangential velocity
        self.Cp = 0             # pressure coeff
        
# the function below discretized the geometry into panels
# the discretization is done with a circle to distribute points
def definePanels(N,xp,yp):
    
    R = (max(xp)-min(xp))/2             # radius of the circle
    xCenter = (max(xp)-min(xp))/2       # x-coord of the center
    xCircle = xCenter + R*np.cos(np.linspace(0,2*pi,N+1)) # x of the circle pts
    
    x = np.copy(xCircle)    # projection of the x-coord on the airfoil
    y = np.empty_like(x)    # initialization of the y coord array
    
    xp,yp = np.append(xp,xp[0]),np.append(yp,yp[0]) # entending arrays
    
    I = 0
    for i in range(N):
        while(I<len(xp)-1):
            if (xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b = yp[I+1]-a*xp[I+1]
        y[i] = a*x[i]+b
    y[N] = y[0]
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
        
    return panel
    
N = 100                      #  number of panels <----------------
panel = definePanels(N,xp,yp)   # discretization of the geometry into panels

# ------------ defining freestream conditions -----------

# class Freestream containg the freestream conditions
class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf = Uinf                # bulk flow velo
        self.alpha = alpha*pi/180       # angle of attack

# defining parameters for above class
Uinf = 100.0                              # freestream velocity
alpha = 0.0                             # angle of attack
freestream = Freestream(Uinf,alpha)     # instant of object freestream

# ---------- building a linear system -----------------

# this system differs from the system in lesson 10 because it incorperates 
# vortices onto the panel to impose circulation to calculate lift which 
# requires the use of the kutta condition and tangential velocity/Cp

# function to evaluate the integral Iij(zi)
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
                +(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
               /((xci-(pj.xa-sin(pj.beta)*s))**2\
               + (yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    
# building the matrix for the source portion
def sourceMatrix(p):
    N = len(p)
    A = np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],+cos(p[i].beta),\
                                    +sin(p[i].beta))
    return A
    # where A is the matrix which is multiplied by the strengths on
    # one side of the linear system of eqn
    
# funciton to build the vortex portion of the system
        
def vortexArray(p):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
        for j in range(N):
            if (j!=i):
                B[i] -= 0.5/pi*I(p[i].xc,p[i].yc,p[j],+sin(p[i].beta),\
                                -cos(p[i].beta))
    return B
    
# and the kutta condition contribution to the array

def kuttaArray(p):
    N = len(p)
    B = np.zeros(N+1,dtype=float)
    for j in range(N):
        if (j==0):
            B[j] = 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],-sin(p[N-1].beta),\
                            +cos(p[N-1].beta))
        elif (j==N-1):
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],-sin(p[0].beta),\
                            +cos(p[0].beta))
        else:
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],-sin(p[0].beta),+cos(p[0].beta))\
                + 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],-sin(p[N-1].beta),+cos(p[N-1].beta))
            B[N] -= 0.5/pi*I(p[0].xc,p[0].yc,p[j],+cos(p[0].beta),+sin(p[0].beta))\
                + 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],+cos(p[N-1].beta),+sin(p[N-1].beta))
    return B
    
# these are the three parts of the linear system
# we now assemble all these into the final matrix, called the global matrix

# function to assemble the global matrix
def buildMatrix(panel):
    N = len(panel)
    A = np.empty((N+1,N+1),dtype=float)
    AS = sourceMatrix(panel)
    BV = vortexArray(panel)
    BK = kuttaArray(panel)
    A[0:N,0:N],A[0:N,N],A[N,:] = AS[:,:],BV[:],BK[:]
    return A
    
# ----------- building linear system -------------    
    
# funciton to build the RHS of the linear system
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N+1,dtype=float)
    for i in range(N):
        B[i] = -fs.Uinf*cos(fs.alpha-p[i].beta)
    B[N] = -fs.Uinf*(sin(fs.alpha-p[0].beta)+sin(fs.alpha-p[N-1].beta))
    return B
    
A = buildMatrix(panel)              # calculating the singularity matrix
B = buildRHS(panel,freestream)      # calculating the freestream RHS

# solving the linear system
var = np.linalg.solve(A,B)
for i in range(len(panel)):
    panel[i].sigma = var[i]
gamma = var[-1]


# ----------- imposing tangential velocity ------------------

# now getting the tangential velocity to impose the kutta condition
def getTangentialVelocity(p,fs,gamma):
    N = len(p)
    A = np.zeros((N,N+1),dtype=float)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),\
                                +cos(p[i].beta))
                A[i,N] -= 0.5/pi*I(p[i].xc,p[i].yc,p[j],+cos(p[i].beta),\
                                +sin(p[i].beta))

    B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
    var = np.empty(N+1,dtype=float)
    var = np.append([pp.sigma for pp in p],gamma)
    vt = np.dot(A,var)+B
    for i in range(N):
        p[i].vt = vt[i]
        
getTangentialVelocity(panel,freestream,gamma) # getting tangential velocity
 
# ------------ BEGINNING OF NEW CODE ----------

# find the control point on the geometry just after the stagnation point
I = 0
for i in range(len(panel)):
    if (panel[i].vt/panel[0].vt<0.):
        I = i
        break

# control point before stagnation point on upper side
V1,x1,y1 = panel[I-1].vt,panel[I-1].xc,panel[I-1].yc
# control point before stagnation point on lower side
V2,x2,y2 = panel[I].vt,panel[I].xc,panel[I].yc

# interpolation to find the location of the stagnation point
xStagn,yStagn = x1-V1*(x2-x1)/(V2-V1) , y1-V1*(y2-y1)/(V2-V1)

sUpper,sLower = np.zeros(I+1),np.zeros(N-I+1)
VeUpper,VeLower = np.zeros(I+1),np.zeros(N-I+1)

sUpper[1] = sqrt((xStagn-panel[I-1].xc)**2+(yStagn-panel[I-1].yc)**2)
VeUpper[1] = -panel[I-1].vt
sLower[1] = sqrt((xStagn-panel[I].xc)**2+(yStagn-panel[I].yc)**2)
VeLower[1] = panel[I].vt

for i in range(1,I):
    sUpper[i+1] = sUpper[i] + panel[I-1-i].length/2
    VeUpper[i+1] = -panel[I-1-i].vt
for i in range(1,N-I):
    sLower[i+1] = sLower[i] + panel[I+i].length/2
    VeLower[i+1] = panel[I+i].vt
    
rho = 1.2                   # density of air kg/m**3
mu = 1.9*10**-5             # dynamic viscosity of air kg/ms
nu = mu/rho                 # kinematic viscosity
    
thetaUpper = np.zeros(len(sUpper),dtype=float)
thetaLower = np.zeros(len(sLower),dtype=float)

# computing the first value at the stagnation point
thetaUpper[0] = sqrt(0.75*nu/abs((VeUpper[1]-VeUpper[0])/(sUpper[1]-sUpper[0])))
thetaLower[0] = sqrt(0.75*nu/abs((VeLower[1]-VeLower[0])/(sLower[1]-sLower[0])))

# integration using the trapezoidal rule: Numpy function np.trapz
for i in range(1,len(thetaUpper)):
    thetaUpper[i] = sqrt(0.45*nu/VeUpper[i]**6*np.trapz(VeUpper[0:i+1]**5,sUpper[0:i+1]))
for i in range(1,len(thetaLower)):
    thetaLower[i] = sqrt(0.45*nu/VeLower[i]**6*np.trapz(VeLower[0:i+1]**5,sLower[0:i+1]))
 
# code credit here and up Olivier Mesnard   
    
dVedsUpper = np.zeros(len(VeUpper),dtype=float)
dVedsLower = np.zeros(len(VeLower),dtype=float)

# defining fucntion to calculate derivative
def gradient(y,x):
    return np.gradient(y,np.append(x[1]-x[0],x[1:len(x)]-x[0:len(x)-1]))

# calculating the gradient of velocity
dVedsUpper = gradient(VeUpper,sUpper)
dVedsLower = gradient(VeLower,sUpper)

# calculating the pressure gradient parameter           
lambdaUpper = np.zeros(len(sUpper),dtype=float)
lambdaLower = np.zeros(len(sLower),dtype=float)

for i in range(len(lambdaUpper)):
    lambdaUpper[i] = ((rho*thetaUpper[i]**2)/mu)*(dVedsUpper[i])
    lambdaUpper[i] = ((rho*thetaUpper[i]**2)/mu)*(dVedsLower[i])
    
# calculating the shape factor from lambda
HUpper = np.zeros(len(sUpper),dtype=float)
HLower = np.zeros(len(sLower),dtype=float)

for i in range(len(HUpper)):
    if lambdaUpper[i]>0 and lambdaUpper[i]<0.1:
        HUpper[i] = 2.61-3.75*lambdaUpper[i]+5.24*lambdaUpper[i]**2
    if lambdaUpper[i]>-0.1 and lambdaUpper[i]<=0:
        HUpper[i] = 2.088+(0.0731/(lambdaUpper[i]+0.14))
        
# calculating the parameter l
lUpper = np.zeros(len(sUpper),dtype=float)
lLower = np.zeros(len(sLower),dtype=float)

for i in range(len(sUpper)):
    if lambdaUpper[i]>0 and lambdaUpper[i]<0.1:
        lUpper[i] = 0.22+1.402*lambdaUpper[i] + (0.018*lambdaUpper[i])/\
                    (lambdaUpper[i]+0.107)
    if lambdaUpper[i]>-0.1 and lambdaUpper[i]<=0:
        lUpper[i] = 2.088+(0.0731)/(lambdaUpper[i]+0.14)
        
# calculating the coefficient of friction
cfUpper = np.zeros(len(sUpper),dtype=float)
cfLower = np.zeros(len(sLower),dtype=float)

for i in range(len(sUpper)):
    cfUpper[i] = 2*lUpper[i]/(VeUpper[i]*thetaUpper[i]/nu)

# -------- Michaels Tranision Criterion -----------

# the criterion is entirely based on the reynolds numbers computed below
ReTheta = np.zeros_like(VeUpper)
ReS = np.zeros_like(ReTheta)
mc = np.zeros_like(ReS)                         # transition criterion   

for i in range(1,len(mc)-1):
    ReTheta[i] = (rho*VeUpper[i]*thetaUpper[i])/mu     # Re based on momentum thickness
    ReS[i] = (rho*VeUpper[i]*sUpper[i])/mu             # Re based on position
    mc[i] = 1.174*(1+(22400/ReS[i]))*ReS[i]**0.46      # transition criterion
    if mc[i]<ReTheta[i]:
        print 'Transition point at: ', sUpper[i]
        iTrans = i
        sTrans = sUpper[i]                      # used to overwrite in head's
        thetaTrans = thetaUpper[i]
        HTrans = HUpper[i-1]
        break
        
# ------------- Head's Method --------------

# defining the four equations to solve for the four unknowns

# skin friction
def cf(theta,H,Ve,nu):
    return 0.246*10.**(-0.678*H)*(theta*Ve/nu)**(-0.268)
    
# entrainment shape factor
def fH1(H):
    if (H<=1.6): 
        return 3.3 + 0.8234*(H-1.1)**(-1.287)
    else: return 3.3 + 1.5501*(H-0.6678)**(-3.064)

# shape factor, the inverse function
def fH(H1):
    H = 1.1 + ((H1-3.3)/0.8234)**(-(1./1.287))
    if (H<=1.6):
        return H
    else: return 0.6678 + ((H1-3.3)/1.5501)**(-(1./3.064))
    
# von Karman momentum integral equation    
def F(H1):
    return 0.0306*(H1-3.)**(-0.6169)
    
def RHS1(theta,H,Ve,nu,dVeds):
    return 0.5*cf(theta,H,Ve,nu) - theta/Ve*(2+H)*dVeds

def RHS2(H1,theta,H,Ve,nu,dVeds):
    return -(H1/theta)*RHS1(theta,H,Ve,nu,dVeds) - (H1/Ve)*dVeds + F(H1)/theta
            
H1Upper = np.zeros(len(sUpper),dtype=float)

# initial condition (at the transition point)
H1Upper[iTrans] = fH1(HTrans)
HUpper[iTrans] = HTrans
thetaUpper[iTrans] = thetaTrans

# advancing to the next point
for i in range(iTrans,len(sUpper)-1):
    h = sUpper[i+1]-sUpper[i]
    thetaUpper[i+1] = thetaUpper[i] + h*RHS1(thetaUpper[i],HUpper[i],\
                                            VeUpper[i],nu,dVedsUpper[i])
    H1Upper[i+1] = H1Upper[i] + h*RHS2(H1Upper[i],thetaUpper[i],HUpper[i],\
                                        VeUpper[i],nu,dVedsUpper[i])
    HUpper[i+1] = fH(H1Upper[i+1])
    cfUpper[i] = cf(thetaUpper[i],HUpper[i],VeUpper[i],nu)
    
    
deltaUpper = np.zeros(len(sUpper),dtype=float)
deltaLower = np.zeros(len(sLower),dtype=float)

# calculating displacement thickness
deltaUpper = HUpper*thetaUpper
 
plt.figure(figsize=(10,6))
plt.plot(sUpper,deltaUpper)

# --------------- Plotting --------------

# creating separate position arrays for plotting
xcLower = [p.xc for p in panel if p.loc=='intrados']
ycLower = [p.yc for p in panel if p.loc=='intrados']

xcUpper = [p.xc for p in panel if p.loc=='extrados']
ycUpper = [p.yc for p in panel if p.loc=='extrados']

betaUpper = [p.beta for p in panel if p.loc=='extrados']
betaLower = [p.beta for p in panel if p.loc=='intrados']

# adding the height of the displacement thickness to the airfoil
xDeltaUpper = np.zeros(len(sUpper),dtype=float)
yDeltaUpper = np.zeros(len(sUpper),dtype=float)

# flipping the order such that it can be mapped
deltaUpper = deltaUpper[::-1]

# the boundary layer is very thin so for visualization purposes it's 
# magnitude is multiplied by 10

for i in range(len(xcUpper)):
    xDeltaUpper[i] = (deltaUpper[i]*np.cos(betaUpper[i]))*20+xcUpper[i]
    yDeltaUpper[i] = (deltaUpper[i]*np.sin(betaUpper[i]))*20+ycUpper[i]

# plotting the discretized geometry
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valY*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 20
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.title('Displacement Thickness')
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(xDeltaUpper,yDeltaUpper,linewidth=2)
plt.show()  

# calculating the cord length
cord = max([p.xa for p in panel])-min([p.xa for p in panel])

# calculating the wall shear stress from the skin friction
tauWall = np.zeros(len(sUpper),dtype=float)
tauWallTotal = np.zeros(len(sUpper),dtype=float)

# integrating the value of wall shear stress over the airfoil surface
for i in range(1,len(sUpper)-1):
    tauWall[i] = cfUpper[i]*0.5*rho*VeUpper[i]
    tauWallTotal[i] = (sUpper[i+1]-sUpper[i])*tauWall[i]
print 'Wall Shear Stress = ', sum(tauWallTotal)

# doulbing the wall shear stress for both sides.
totalWallShear = 2*sum(tauWallTotal)

# calculating drag coefficient
Cd = 2*(totalWallShear/(.5*Uinf**2*cord))

print 'Coefficient of drag = ', Cd

# calculating airfoil Reynold number
Re = Uinf*(cord)/nu
print 'Reynolds Number = ', Re

