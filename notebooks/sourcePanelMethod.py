# lesson 10 - Source Panel Method
# Started 25 Mar 2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# ------------ defining airfoil shape --------------

# data from NACA is in resources file
# reading shape data from naca file
coords = np.loadtxt(fname='/Users/ian/GitHub/AeroHydro/resources/naca0012.dat')
xp,yp = coords[0:,0],coords[:,1]            # separating data into two arrays

# plotting the geometry
# in notebook %matplotlib inline - command to show figures embedded in notebook

valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.show()

# ----------- creating panels on airfoil ---------------

# defining a class to define a panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                 # first point
        self.xb,self.yb = xb,yb                 # second point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2   # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)# length of panel
        
        # orientation of panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel (on outer or inner curve)
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.sigma = 0.                         # source strength
        self.vt = 0.                            # tangential velocity
        self.Cp = 0.                            # pressure coeff


# defining function to perform discreization
def definePanels(N,xp,yp):
    
    R = (max(xp)-min(xp))/2                         # radius of the circle
    xc,yc = (max(xp)+min(xp))/2,(max(yp)+min(yp))/2 # center point of circle
    xCircle = xc + R*np.cos(np.linspace(0,2*pi,N+1))
    yCircle = yc + R*np.sin(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xCircle[0:-1])
    y = np.empty_like(x)
    
    I = 0
    for i in range(N):
        while(I<len(xp)-1):
            if (xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I+= 1
        a = (yp[(I+1)%len(yp)]-yp[I])/(xp[(I+1)%len(yp)]-xp[I])
        b = yp[(I+1)%len(yp)]-a*xp[(I+1)%len(xp)]
        y[i] = a*x[i]+b
        
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[(i+1)%N],y[(i+1)%N])
    
    return panel
    
N = 20                      # number of panels
panel = definePanels(N,xp,yp) # discretization of the geometry into panels

# plotting the geometry with the panels
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
xStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
    np.append([p.ya for p in panel],panel[0].ya),\
    linestyle='-',linewidth=1,\
    marker='o',markersize=6,color='r')
plt.show()

# -------------- establishing freestream conditions -------------

# class for freestream
class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf = Uinf            # freestream velocity mag
        self.alpha = alpha*pi/180   # angle of attack conversion
        
# freestream parameters
Uinf = 1.0                          # velocity mag
alpha = 5.0                         # angle of attack
freestream = Freestream(Uinf,alpha)

# ------------- imposing flow tangency condition --------------

# function to evaluate the integral Iij(zi)
def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
                 +(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
                /((xci-(pj.xa-sin(pj.beta)*s))**2\
                + (yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]

# ------------ building the linear system of equ -------------

# function to build the source matrix
def buildMatrix(p):
    N = len(p)
    A = np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],cos(p[i].beta),\
                sin(p[i].beta))
    return A

# function to make RHS of linear system
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
        B[i] = -fs.Uinf*cos(fs.alpha-p[i].beta)
    return B

A = buildMatrix(panel)          # calculating the singularity matrix
B = buildRHS(panel,freestream)  # calculating freestream RHS

# solving the linear system
var = np.linalg.solve(A,B)
for i in range(len(panel)):
    panel[i].sigma = var[i]
    
# ----------- calculating the surface pressure coeff -----------

# function to calculate the tangential velocity to the control points
def getTangentialVelocity(p,fs,gamma):
    N = len(p)
    A = np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if (i!=j):
                A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),\
                        cos(p[i].beta))
    B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
    var = np.array([pp.sigma for pp in p])
    vt = np.dot(A,var)+B
    for i in range(N):
        p[i].vt = vt[i]

getTangentialVelocity(panel,freestream,gamma)

# defining function to get pressure coefficient
def getPressureCoefficient(p,fs):
    for i in range(len(p)):
        p[i].Cp = 1-(p[i].vt/fs.Uinf)**2
        
getPressureCoefficient(panel,freestream)

# plotting pressure coefficient
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]),max([p.Cp for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
        [p.Cp for p in panel if p.loc=='extrados'],\
        'ro-',linewidth=2)
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
        [p.Cp for p in panel if p.loc=='intrados'],
        'bo-',linewidth=2)
plt.legend(['extrados','intrados'],'best',prop={'size':14})
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.gca().invert_yaxis()
plt.title('Number of panels : %d'%len(panel))
plt.show()

# ------- producing/plotting stream line plot --------

# function to calculate velocity field given a mesh
def getVelocityField(panel,freestream,gamma,X,Y):
    Nx,Ny = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.Uinf*cos(freestream.alpha)\
                + 0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,1,0) for p in panel])
            v[i,j] = freestream.Uinf*sin(freestream.alpha)\
                + 0.5/pi*sum([p.sigma*I(X[i,j],Y[i,j],p,0,1) for p in panel])
    return u,v
    
# defining the mesh
Nx,Ny = 20,20
valX,valY = 1.0,2.0
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
X,Y = np.meshgrid(np.linspace(xStart,xEnd,Nx),np.linspace(yStart,yEnd,Ny))

# get the velocity field on the mesh grid
u,v = getVelocityField(panel,freestream,gamma,X,Y)

# plotting the velocity field
size=12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=1,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of the velocity field')

# ---------- calculating/plotting pressure contour --------

Cp = 1.0-(u**2+v**2)/freestream.Uinf**2

# plotting the pressure field
size=12
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.fill([p.xc for p in panel],[p.yc for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of the pressure field')
plt.show()
