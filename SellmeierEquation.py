# Author: Pieters Michaël
# Date: 28/11/2020

import matplotlib.pyplot as plt
import numpy as np
from math import *

#########################
## INPUT PARAMETERS #####
#########################

# Data for SiO2

G_1 = 0.696749                  # Sellmeier coefficient
G_2 = 0.408218                  # Sellmeier coefficient
G_3 = 0.890815                  # Sellmeier coefficient

lambda_1 = 0.0690606            # Sellmeier coefficient (in micrometer)
lambda_2 = 0.115652             # Sellmeier coefficient (in micrometer)
lambda_3 = 9.900559             # Sellmeier coefficient (in micrometer)


#########################
## SOLUTION #############
#########################

lambda_p = np.linspace(0.5, 1.8, 101)

term1 = (G_1*lambda_p**2)/(lambda_p**2 - lambda_1**2)
term2 = (G_2*lambda_p**2)/(lambda_p**2 - lambda_2**2)
term3 = (G_3*lambda_p**2)/(lambda_p**2 - lambda_3**2)

ref_index_sq = 1.0 +term1 + term2 + term3
ref_index = np.sqrt(ref_index_sq)

plt.plot(lambda_p, ref_index, linewidth=1.5)
plt.title('Sellmeier Equation')
plt.xlabel('$\lambda$ in $\mu m$')
plt.ylabel('Refractive index for Si02')
plt.show()
plt.close()
