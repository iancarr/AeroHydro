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
    
N = 20                      #  number of panels <----------------
panel = definePanels(N,xp,yp)   # discretization of the geometry into panels

# ------------ defining freestream conditions -----------

# class Freestream containg the freestream conditions
class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf = Uinf                # bulk flow velo
        self.alpha = alpha*pi/180       # angle of attack

# defining parameters for above class
Uinf = 1.0                              # freestream velocity
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

# plotting tangential velocity
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
vtmin,vtmax = min([p.vt for p in panel]),max([p.vt for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = vtmin-valY*(vtmax-vtmin),xmax+valY*(vtmax-vtmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$V_t$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
         [-p.vt for p in panel if p.loc=='extrados'],'-bo')
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
         [p.vt for p in panel if p.loc=='intrados'],'-ro')
plt.show()
 
# ------------ BEGINNING OF NEW CODE ----------

print testGrad
# defining function to get velocity gradient
N = len(panel)/2
dvdxInt = np.zeros(N)
dvdxExt = np.zeros_like(dvdxInt)
for i in range(N-2):
    if panel[i].loc=='intrados':
        print i
        a = (panel[i+2].xc-panel[i].xc)/(panel[i+1].xc-panel[i].xc)
        dvdxInt[i] = (1/(panel[i+2].xc-panel[i+1].xc))*(a*(panel[i+1].vt-panel[i].vt)\
                -((panel[i+2].vt-panel[i].vt)/a))

for i in range(N-2):
    if panel[i].loc=='extrados':
        print i
        a = (panel[i+2].xc-panel[i].xc)/(panel[i+1].xc-panel[i].xc)
        dvdxExt[i] = (1/(panel[i+2].xc-panel[i+1].xc))*(a*(panel[i+1].vt-panel[i].vt)\
                        -((panel[i+2].vt-panel[i].vt)/a))

                
# combining the two velocity gradient parts 
dvdxInt = np.append(dvdxInt,np.zeros(2))

# plotting
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('Airfoil Surface',fontsize=16)
plt.ylabel('Velocity Gradient',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='intrados'], dvdxInt,'-ro')
plt.plot([p.xc for p in panel if p.loc =='extrados'], dvdxExt,'-bo')
plt.show()

# calculating Reynolds number based on freestream velocity
# the material properties are based on air
rho = 1.2                   # density of air kg/m**3
mu = 1.9*10**-5             # dynamic viscosity of air
nu = mu/rho                 # kinematic viscosity
L = max(xp)-min(xp)
Re = rho*Uinf*L/mu

# calculating integral in the momentum thickness equation
intVe = np.zeros_like(dvdx,dtype=float) # integral in theta calculation
theta = np.zeros_like(intVe,dtype=float)

# intermediate calcualtion of integral
for i in range(N-1):
    if panel[i].loc=='intrados':
        panel[i].vt = -panel[i].vt   
    intVe[i] = (panel[i].xc-panel[0].xc)*((panel[i].vt+panel[0].vt)/2)**5

# calcating momentum thickness
for i in range(N-1):
    if panel[i].loc=='intrados':
        panel[i].vt = -panel[i].vt
    theta[i] = np.sqrt((0.45/panel[i].vt**6)*intVe[i])
    
# assigning value of theta at stagnation
theta[0] = np.sqrt((0.075*nu)/dvdx[0])

# calculating the pressure gradient parameter
lam = np.empty_like(theta)
for i in range(N):
    lam[i] = (rho*theta[i]**2/mu)*(dvdx[i]) 
    
# calculating shape factor H from lambda (lam)
