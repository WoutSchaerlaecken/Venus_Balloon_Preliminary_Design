import numpy as np



n = 460
T = 300
p = 300 * 10**5
SF = 2


"-----Tank Specifications-----"
rho = 2200
tensile_strength = 3 * 10**9

R = 8.31446261815324

V = n*R*T/p
r = ((V * 3)/(4*np.pi))**(1/3)



t = p * r /(2*tensile_strength) * SF
print("Thickness: ", t)
mass = rho * 4 * np.pi * r**2 * t
print("Mass: ", mass)

