#Integration der d'Alembertgleichung

# import

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.integrate import odeint


# const SI wenn nicht anders deklariert

c = 299792458 # m/s
pi = np.pi
h = 6.62607015e-34
hbar = h/(2*pi)
k_B = 1.380649e-23
G = 6.674e-11
yr_SI = 60*60*24*365.25
pc = 1.029e8*c


############################
my = 10e12 # Grössenordnung 1/Compton 
n = 1#1e38# Grössenordnung? extrem gross
H0 = 67.7 # H(a=1) unit: km/sMpc
H0_SI = H0*1e-3/pc
L_H = c/H0_SI #Hubble Länge


print('H0 = %.4e s^-1' %H0_SI)
print('Hubble-Länge = %.4e' %L_H)

###############################
R = L_H
#ny2 = (c**2 * n*(n+2) * R**4 + c**2 * my**2 * R**6)

ny2 = (H0_SI*R/c)**4 * n*(n+2)  #+ (H0/c)**4 * my**2 * R**6 #passt etwa von der grössenordnung mit schrödingers näerung

#ny2 = H0**2 * R**4 * c**(-1) * n*(n+2)  + H0**2 * c**(-1) * my**2 * R**6
#ny2 = (H0*R/c)**4 * n*(n+2)  * my**(-2)+ (H0/c)**4  * R**6

# ny1 = R**3 *c *(2*pi)**(-1) * (n*(n+2)/R**2 + my)**(-0.5)# aus skript Gl.11

print(ny2**-0.5) # kehrwert von ny
print(R)
print(n*(n+2))

#z = c**3 * H0**(-2) * R**(-3) # z ist tau-schlange
    
# DGL (tau_schlange)

def fkt(y, z):
    f, u = y
    dfdz = u
    dudz = -ny2 *f
    
    dydz = [dfdz, dudz]
    return dydz



z_lst = np.linspace(0,10,1000)
y0 = [0, 1]

#Lsg der DGL
sol = odeint(fkt, y0, z_lst)

plt.plot(z_lst, sol[:, 0], label = 'f_r') 
#plt.plot(z_lst, sol[:, 0].imag, label = 'f_i')
#plt.plot(z_lst, sol[:, 1], label = 'u')
plt.legend()
plt.xlabel('tau_schlange')
plt.ylabel('f')
plt.show()

    
    
    
    
    