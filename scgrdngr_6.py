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
n =  1#1e38# Grössenordnung? extrem gross
H0 = 67.7 # H(a=1) unit: km/sMpc
H0_SI = H0*1e-3/pc
L_H = c/H0_SI #Hubble Länge


#print(pc)
print('H0 = %.4e s^-1' %H0_SI)
print('Hubble-Länge = %.4e' %L_H)

###############################
#R = L_H*(1+z)
b = H0_SI
k2 = c**2 * n*(n+2)/ b**2 -1
k = k2**0.5
print(k2, k)

def fkt(y, x):
    w, u = y
    R = L_H*(1+x)
    #z = 1+b*x
    z = my*R*c/b
    dwdz = u
    dudz = -1/z*u - (1+k**2/z**2)*w

    
    dydz = [dwdz, dudz]
    return dydz


#z_lst = np.linspace(0,1e61,10)
z_lst = np.linspace(0,10,1000) #nicht bei 0 beginnen --> division by zero
y0 = [0, 1] 

sol = odeint(fkt, y0, z_lst)
#print(sol)

plt.plot(z_lst, sol[:, 0], label = 'w_r') #subplts benutzen
plt.plot(z_lst, sol[:, 0].imag, label = 'w_i')
plt.plot(z_lst, sol[:, 1], label = 'u')
plt.legend()
plt.xlabel('z~')#bzw z-schlange
plt.ylabel('f')
plt.show()         

    
    
    
    
    