# notebook 8 -- Source Sheet
# started 6mar2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 100
xStart,xEnd = -1.0,1.0
yStart,yEnd = -1.5,1.5
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

# defining the velocity and stream function of a freestream
Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros_like(uFreestream)

psiFreestream = Uinf*Y

# defining a class which describes a source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    # get source velocity

