from astropy.analytic_functions import blackbody_lambda
from astropy import units as u
import numpy as np
from scipy import integrate

# Radius of star
radius= 1 # solar radii

# distance of the star from earth
distance = 10 # parsecs

# temperature of black body
temperature = 6000 * u.K # Kelvin

'''
   Constants
'''
# U band
# units: Angstrom, Angstrom, ergs/( cm**2 * s * sr)
u_start, u_end, u_F_0 = 3320, 3980, 1.387 * 10**-6

# B band
# units: Angstrom, Angstrom, ergs/( cm**2 * s * sr)
b_start, b_end, b_F_0= 3980 , 4920, 1.570 * 10**-6

# V band
# units: Angstrom, Angstrom, ergs/( cm**2 * s * sr)
v_start, v_end, v_F_0 = 5070, 5950, 1.000 * 10**-6

# R band
# units: Angstrom, Angstrom, ergs/( cm**2 * s * sr )
r_start, r_end, r_F_0 = 5890, 7270, 1.054 * 10**-6


# light year 'LY'
LY = 9.4607 * 10**15 # m

# parsec 'parsec'
parsec = 3.26 * LY # m

# solar radius 'SR'
SR = 6.95700 * 10**8 # m

'''
#     flux is a function of lambda and Temperature
#                        2 * h * c^2                   1                       R^2
#     F ( lambda , T ) = ---------------  *  ----------------------------- * --------
#                           lambda^5          e^(h*c / (lambda*k_b*T) - 1      r^2
#     where
#       h is the Plank's Constant
#       c is the speed of light
#       k_b is the Boltzmann's Constant
#       R is the radius of the black body
#       r is the distance from the black body
'''
flux = lambda x: (
   float( blackbody_lambda(x,temperature).value ) *
   float( (float(SR * radius)/float(parsec * distance) )**2 ) )

# using the formula for magnitude
# m = m_0 - 2.5 * logbase10 ( F / F_0 )
u_magnitude= -2.5*np.log10(integrate.quad(flux,u_start,u_end)[0]/u_F_0)
print "Approximate magnitude when U filter is applied: " + str(u_magnitude)

# using the formula for magnitude
# m = m_0 - 2.5 * logbase10 ( F / F_0 )
b_magnitude= -2.5*np.log10(integrate.quad(flux,b_start,b_end)[0]/b_F_0)
print "Approximate magnitude when B filter is applied: " + str(b_magnitude)

# using the formula for magnitude
# m = m_0 - 2.5 * logbase10 ( F / F_0 )
v_magnitude= -2.5*np.log10(integrate.quad(flux,v_start,v_end)[0]/v_F_0)
print "Approximate magnitude when V filter is applied: " + str(v_magnitude)

# using the formula for magnitude
# m = m_0 - 2.5 * logbase10 ( F / F_0 )
r_magnitude= -2.5*np.log10(integrate.quad(flux,r_start,r_end)[0]/r_F_0)
print "Approximate magnitude when R filter is applied: " + str(r_magnitude)


