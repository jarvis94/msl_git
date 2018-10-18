#!/usr/bin/env python

#Solve Poisson's equation
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg

# Parameters for solving problem in the interval 0 < x < L
N = 120                    # No. of coordinate points
L = 6
x,dx = np.linspace(0,L,N,retstep=True);       # Coordinate vector
rho=np.zeros(N)

#Three-point finite-difference representation of Laplacian
Lap = (-2*np.diag([1]*N,k=0)+np.diag([1]*(N-1),k=1)+np.diag([1]*(N-1),k=-1))#/(dx**2)
print(Lap)
# Next modify Lap so that it is consistent with f(0) = f(L) = 0
Lap[0,0] = -2
#Lap[0,2] = 1
#Lap[1,0] = 1

#Lap[N-1,N-2] = 1; 
#Lap[N-2,N-1] = 1; 
Lap[N-1,N-1] = -1;

#Lap=np.delete(Lap,(0),axis=0) #Delete first row
#Lap=np.delete(Lap,(N-2),axis=0) #Delete last row


myLap=np.matrix(Lap)
print(myLap)

#Potential
##
# rho[0]=1.0
# rho[1:int(N/2)]=1.0
# rho[N-1]=-1.0/2
# rho[int(N/2):N-1]=-1.0
###
# rho[0:int(N/4)]=0.0
# rho[int(N/4):int(N/2)]=1.0
# rho[int(N/2):int(3*N/4)]=-1.0
# rho[N-1]=0.0
###
rho[0:int(N/4)]=0.0
rho[int(N/4):int(N/2)]=1.0
rho[int(N/2):int(3*N/4)]=-1.0
rho[N-1]=0.0



#rho=np.delete(rho,(0),axis=0) #Delete first row
#rho=np.delete(rho,(N-2),axis=0) #Delete last row
myrho=np.matrix(rho)
print(myrho)

plt.plot(rho,'g*-')
plt.grid()
plt.show()
print(myrho.T)

V = -myLap.I*myrho.T
print(V)
plt.plot((V),'r+-')
plt.grid()
plt.show()
