# Author: Pieters Michaël
# Date: 30/11/2020

import matplotlib.pyplot as plt
import numpy as np
from math import *

#########################
## INPUT PARAMETERS #####
#########################

N_x = 200                       # number of x-cells
N_time = 100                    # number of time steps
time = 0.0
E_y = np.zeros(N_x)             # starting values for y-comp of electric field
H_z = np.zeros(N_x)             # starting values for z-comp of magnetic field

t_zero = 40.0                   # center of incident pulse
width = 12                      # width of incident pulse


#########################
## SOLUTION #############
#########################

for step in range(N_time):
    time += 1
    pulse = exp(-0.5*((t_zero - time)/width)**2)    # gaussian pulse
    E_y[120] = pulse
    # vectorized FTDT scheme in free space
    E_y[1:N_x] = E_y[1:N_x] - 0.5*(H_z[1:N_x] - H_z[0:N_x-1])
    H_z[0:N_x-1] = H_z[0:N_x-1] - 0.5*(E_y[1:N_x] - E_y[0:N_x-1])
    
        
x = np.arange(N_x)

fig = plt.figure(figsize = [8.0, 10.0])
ax1 = plt.subplot(2, 1, 1)
ax1.plot(x, E_y, linewidth=1.5)
ax1.set_xlim(0, 200)
ax1.set_ylim(-1, 1)
plt.xlabel('FTDT cells')
plt.ylabel('y-component of E-field')
plt.title('Propagation of a Gaussian Pulse')
ax2 = plt.subplot(2, 1, 2)
ax2.plot(x, H_z, linewidth=1.5)
ax2.set_xlim(0, 200)
ax2.set_ylim(-1, 1)
plt.xlabel('FTDT cells')
plt.ylabel('z-component of H-field')
plt.show()
plt.close(fig)
