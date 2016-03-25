from astropy.analytic_functions import blackbody_lambda
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np


# temperature of black body
temperature= 6000 * u.K

#     Plots B as a function of lambda
#                        2 * h * c^2                   1
#     B ( lambda , T ) = ---------------  *  -----------------------------
#                           lambda^5          e^(h*c / (lambda*k_b*T) - 1


# wein's constant
b_wein= (2.8977729* 10**-3)*u.m*u.K

# finding the maximum wavelength using wein's formula
wavelength_max= (b_wein / temperature).to(u.AA)

# setting up the x axis upto the maximum wavelength
wavelengths=np.logspace(0
                       ,np.log10(wavelength_max.value + 10 * wavelength_max.value)
                       ,num=1000)* u.AA

#using the function blackbody_lambda from the astropy library
flux = blackbody_lambda(wavelengths, temperature)

# plot flux vs the wavelength
plt.plot(wavelengths.value, flux.value)

plt.xlabel('wavelength ( Angstroms )')

plt.ylabel('Flux ( ' + str(flux.unit) +' )' )

#show the plot
plt.show()