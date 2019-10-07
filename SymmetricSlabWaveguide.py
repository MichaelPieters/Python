# Author: Pieters Michaël
# Date: 7/10/2019

import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy.optimize import ridder
from scipy.integrate import quad

#########################
## INPUT PARAMETERS #####
#########################

nf = 3.6                        # breaking index of film
n0 = 3.55                       # breaking index of substrate
a = 2.5                         # halfwidth thickness (in micrometer)
wavelength = 1.0                # vacuum wavenumber (in micrometer)

#########################
## CONSTANTS ############
#########################

Z = 377                         # impedance of free space (in ohms)


#########################
## FUNCTIONS ############
#########################

k = 2*pi/wavelength
R = a*k*sqrt(nf**2 - n0**2)

def odd_modes(y):
    return y/tan(y) + sqrt(R**2 - y**2)

def even_modes(y):
    return y*tan(y) - sqrt(R**2 - y**2)

def kappa(y):
    return y/a

def gamma(y):
    return np.sqrt(R**2 - y**2)/a

def TE_ODD(x, kappa, gamma):
    if x > a:
        Ey = sin(kappa*a)*exp(-gamma*(x-a))
    elif x < -a:
        Ey = -sin(kappa*a)*exp(gamma*(x+a))
    else:
        Ey = sin(kappa*x)
    return Ey

def TE_EVEN(x, kappa, gamma):
    if x > a:
        Ey = cos(kappa*a)*exp(-gamma*(x-a))
    elif x < -a:
        Ey = cos(kappa*a)*exp(gamma*(x+a))
    else:
        Ey = cos(kappa*x)
    return Ey

def POWER_ODD(x, kappa, gamma):
    beta = sqrt(gamma**2 + n0**2*k**2)
    return (beta/(2*k*Z))*TE_ODD(x, kappa, gamma)**2

def POWER_EVEN(x, kappa, gamma):
    beta = sqrt(gamma**2 + n0**2*k**2)
    return (beta/(2*k*Z))*TE_EVEN(x, kappa, gamma)**2

TE_ODD_VEC = np.vectorize(TE_ODD)
TE_EVEN_VEC = np.vectorize(TE_EVEN)

#########################
## SOLUTION #############
#########################

n = 0

X = np.linspace(-10, 10, 500)

while n*pi/2 < R:
    yi = n*pi/2
    yf = min((n+1)*pi/2 - 0.001, R)
    if n%2 == 0:
        root = ridder(even_modes, yi, yf)
        A = quad(POWER_EVEN, -10, 10, args=(kappa(root), gamma(root)))[0]
        Y = TE_EVEN_VEC(X, kappa(root), gamma(root))/sqrt(A)
    elif n%2 == 1:
        root = ridder(odd_modes, yi, yf)
        A = quad(POWER_ODD, -10, 10, args=(kappa(root), gamma(root)))[0]
        Y = TE_ODD_VEC(X, kappa(root), gamma(root))/sqrt(A)
    plt.plot(X, Y, linewidth=1.5, label=str(n))
    n += 1
        

plt.legend(loc=4)
plt.title('Symmetric Slab Waveguide Modes')
plt.xlabel('$x$ in $\mu m$')
plt.show()
plt.close()

