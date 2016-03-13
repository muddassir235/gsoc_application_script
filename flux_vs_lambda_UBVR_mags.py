from astropy.analytic_functions import blackbody_lambda
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

# Radius of black body
radius=1

# temprature of black body
temprature= 1000 * u.K

'''
    Constants
'''
# U band
# units: Angstrom, Angstrom, ergs/(Angstrom * cm**2 * s)
u_start, u_end, u_F_0 = 3320, 3980, 1.81* 10**-6

# B band
# units: Angstrom, Angstrom, ergs/(Angstrom * cm**2 * s)
b_start, b_end, b_F_0= 3980, 4920, 4260* 10**-9

# V band
# units: Angstrom, Angstrom, ergs/(Angstrom * cm**2 * s)
v_start, v_end, v_F_0 = 5070, 5950, 3640* 10**-9

# R band
# units: Angstrom, Angstrom, ergs/(Angstrom * cm**2 * s)
r_start, r_end, r_F_0 = 5890, 7270, 3080* 10**-9

# speed of light 'c'
# units: m / s
c = 299792458

# light year 'LY'
# units: m
LY = 9.4607 * 10**15

# parsec 'parsec'
# units: m
parsec = 3.26 * LY

# solar radius 'SR'
# units: m
SR = 6.95700 * 10**8

'''
#     flux is a function of lambda and Temprature
#                        2 * h * c^2                   1                       R^2
#     F ( lambda , T ) = ---------------  *  ----------------------------- * --------
#                           lambda^5          e^(h*c / (lambda*k_b*T) - 1      r^2
#     where
#       h is the Plank's Constant
#       c is the speed of light
#       k_b is the Boltzman's Constant
#       R is the radius of the black body
#       r is the distance from the black body
'''
flux = lambda x: ( float(blackbody_lambda(x,temprature).value) * float( (float(SR*radius)/float(parsec*10))**2 ) )

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

#     Plots B as a function of lambda
#                        2 * h * c^2                   1
#     B ( lambda , T ) = ---------------  *  -----------------------------
#                           lambda^5          e^(h*c / (lambda*k_b*T) - 1
b_wein= (2.8977729* 10**-3)*u.m*u.K
wavelength_max= (b_wein / temprature).to(u.AA)
wavelengths=np.logspace(0,np.log10(wavelength_max.value + 10* wavelength_max.value),num=1000)* u.AA
flux = blackbody_lambda(wavelengths, temprature)
plt.plot(wavelengths.value, flux.value)
plt.xlabel('wavelength ( Angstroms )')
plt.ylabel('Flux ( ' + str(flux.unit) +' )' )
plt.show()




























