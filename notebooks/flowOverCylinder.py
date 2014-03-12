# notebook 9 -- Flow over a cylinder source sheet method 
# started 8mar2014
# Ian Carr

# importing familiar libraries
import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# defining a cylinder
R = 1.0
theta = np.linspace(0,2*pi,100)
xCylinder,yCylinder = R*np.cos(theta),R*np.sin(theta)

# plotting the cylinder
size = 4
plt.figure(num = None,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.plot(xCylinder,yCylinder,c='b',ls='-',lw=2)
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1)

# creating class for source panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya # defining first point
        self.xb,self.yb = xb,yb # second point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2 # center control point
        self.length = np.sqrt((xb-xa)**2+(yb-ya)**2)
    
    # orinetation of the panel
    if (xb-xa<=0): self.beta = acos((yb-ya)/self.length)
    elif (xb-xa>0): self.beta = pi+acos(-(yb-ya)/self.length)
    
    self.sigma = 0 # source strength
    self.Vt = 0 # tangential velocity
    self.Cp = 0 # pressure coefficient
    
Np =10 # number of panels

# defining end points of panels
xb = R*np.cos(np.linspace(0,2*pi,Np+1))
yb = R*np.sin(np.linspace(0,2*pi,Np+1))

# defining the points on panels between end points
panel = np.empty(Np,dtype=object)
for i in range(Np):
    panel[i] = Panel(xb[i],yb[i],xb[i+1],yb[i+1])
    
# plotting
size = 6
plt.figure(num=None,figsize=(size,size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.plot(xCylinder,yCylinder,c='b',ls='-',lw=1)
plt.plot(xb,yb,c='r',ls='-',lw=1)
plt.scatter([p.xa for p in panel],[p.ya for p in panel],c='r',s=30)
plt.scatter([p.xc for p in panel],[p.ya for p in panel],c='r',s=30)
plt.legend(['cylinder','panels','end point','center point'],\
            loc='best',prop={'size':16})
plt.xlim(-1.1,1.1)
plt.ylim(-1.1,1.1)