from astropy.analytic_functions import blackbody_lambda
import astropy.units as u
from scipy import integrate

'''
  Constants
'''
# light year 'LY'
# units: m
LY = 9.4607 * 10**15

# parsec 'parsec'
# units: m
parsec = 3.26 * LY

# solar radius 'SR'
# units: m
SR = 6.95700 * 10**8

# U band
# Angstroms, Angstroms
U_start, U_end = 3320, 3980

# B band
# Angstroms, Angstroms
B_start, B_end= 3980, 4920

# V band
# Angstroms, Angstroms
V_start, V_end = 5070, 5950

# R band
# Angstroms, Angstroms
R_start, R_end = 5890, 7270


# temperature of Vega
temperature = 9602 * u.K #Kelvin

# radius of Vega
radius = 2.818 # solar radii

# distance of Vega from the earth
distance_from_earth = 25.05 # light years

flux_Vega = lambda x: (
   float(blackbody_lambda(x,temperature).value) *
   float( (float(SR*radius)/float(LY*distance_from_earth))**2 )
)

print 'Flux in U band of Vega: ' + \
     str( integrate.quad(flux_Vega,U_start,U_end)[0] )
print 'Flux in B band of Vega: ' + \
     str( integrate.quad(flux_Vega,B_start,B_end)[0] )
print 'Flux in V band of Vega: ' + \
     str( integrate.quad(flux_Vega,V_start,V_end)[0] )
print 'Flux in R band of Vega: ' + \
     str( integrate.quad(flux_Vega,R_start,R_end)[0] )