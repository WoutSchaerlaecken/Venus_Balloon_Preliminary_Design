from functions import *
from inputs import *
import numpy as np
from constants import *


gondola_mass = payload_mass / probe.payload_percentage + 3
print("gondola mass: ", gondola_mass)

total_mass = gondola_vs_total_mass(gondola_mass)
print("total mass: ", total_mass)


for i in range(20):
    n_zp_2 = Zero_pressure_moles(total_mass, T2 , V_SP, p2 ,R_Venus , R)
    #print("Zero pressure moles: ", n_zp_2)

    V_zp = ZP_balloon_volume(n_zp_2, R, T2, p2)
    #print("Zero pressure volume: ", V_zp)
    #print("diameter: ", (V_zp*6/np.pi)**(1/3))

    balloon_mass = Balloon_mass(V_zp, balloon_material_density)
    print("Balloon mass: ", balloon_mass)

    total_mass = gondola_mass + balloon_mass 
    print("total mass: ", total_mass)
print("Final total mass: ", total_mass)


import matplotlib.pyplot as plt

# Initialize lists to store iteration data
iteration = []
balloon_masses = []
total_masses = []

# Re-run the loop to collect data for plotting
total_mass = gondola_vs_total_mass(gondola_mass)
for i in range(20):
    n_zp_2 = Zero_pressure_moles(total_mass, T2, V_SP, p2, R_Venus, R)
    V_zp = ZP_balloon_volume(n_zp_2, R, T2, p2)
    balloon_mass = Balloon_mass(V_zp, balloon_material_density)
    total_mass = gondola_mass + balloon_mass 

    # Append data to lists
    iteration.append(i)
    balloon_masses.append(balloon_mass)
    total_masses.append(total_mass)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(iteration, balloon_masses, label="Balloon Mass (kg)", marker='o')
plt.plot(iteration, total_masses, label="Total Mass (kg)", marker='x')
plt.xlabel("Iteration")
plt.ylabel("Mass (kg)")
plt.title("Balloon and Total Mass vs Iteration")
plt.legend()
plt.grid()
plt.show()

""" n_zp_1 = Zero_pressure_moles(total_mass, T1 , V_SP, p1 ,R_Venus , R)
print("Zero pressure moles: ", n_zp_1)

print(n_zp_2-n_zp_1)
print((p1/T1 - p2/T2)*V_SP/R)

V_zp_min = ZP_balloon_volume(n_zp_1, R, T2, p2)
print("Minimum Zero pressure volume: ", V_zp_min)

n_sp_2 = (1.5*p2*V_SP) / (R * T2) 
print("Super pressure moles: ", n_sp_2)

super_p_1 = (n_zp_2-n_zp_1+n_sp_2 - (p1*V_SP)/(R*T1)) * ((R*T1) / (V_SP))
print("Super pressure 1: ", super_p_1)

n_sp_1 = (p1+super_p_1*V_SP)/(R * T1)
print("Super pressure moles: ", n_sp_1)



plot_balloons(V_SP, V_zp, V_zp_min) """


