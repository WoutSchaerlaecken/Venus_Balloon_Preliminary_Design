from functions import *
from inputs import *
import numpy as np
from constants import *
import matplotlib.pyplot as plt
from Venus_atmosphere import venus_atmosphere

H1 = 45
H2 = 55


T1, p1,_ = venus_atmosphere(H1)
T2, p2,_ = venus_atmosphere(H2)

balloon_ratio = 1/3


inflation_mass = 2
sample_collection_mass = 1

gondola_mass = payload_mass / probe.payload_percentage + sample_collection_mass + inflation_mass

total_mass = 2*gondola_mass
DP1 = 3000
DP2 = T2/T1 * DP1

iterations = [0]
ZP_balloon_masses = [0]
SP_balloon_masses = [0]
balloon_masses = [0]

total_masses = [total_mass]


for i in range(20):
    n_sp1 = ((p1 + DP1)*V_SP)/(R * T1)
    #print("n_sp: ", n_sp1)

    n_sp2 = ((p2 + DP2)*V_SP)/(R * T2)
    #print("n_sp: ", n_sp2)

    delta_n = n_sp1 - n_sp2
    #print("delta_n: ", delta_n)

    n_zp1 = (R_Venus/R)* (total_mass - p1/(R_Venus*T1)*V_SP)
    #print("n_zp_1: ", n_zp1)
    n_zp2 = (R_Venus/R)* (total_mass - p2/(R_Venus*T2)*V_SP)
    #print("n_zp_2: ", n_zp2)

    delta_n = n_zp2 - n_zp1
    #print("delta_n: ", delta_n)
    V_zp1 = n_zp1 * R * T1 / p1 + V_SP
    V_zp = n_zp2 * R * T2 / p2 + V_SP
    r_zp = (3 * V_zp / (4 * np.pi)) ** (1/3)
    
    r_sp = r_zp * balloon_ratio

    V_SP = 4/3 * np.pi * r_sp**3
    #print("V_zp: ", V_zp)

    ZP_balloon_mass = Balloon_mass(V_zp, balloon_material_density)
    #print("balloon mass: ", ZP_balloon_mass)

    SP_Balloon_mass = Balloon_mass(V_SP, SP_balloon_material_density)
    #print("balloon mass: ", SP_Balloon_mass)

    He_mass = (n_sp1 + n_sp1)*Mhe

    total_mass = gondola_mass + ZP_balloon_mass+ SP_Balloon_mass + He_mass
    #print("total mass: ", total_mass)

    iterations.append(i+1)
    ZP_balloon_masses.append(ZP_balloon_mass)
    SP_balloon_masses.append(SP_Balloon_mass)
    balloon_masses.append(ZP_balloon_mass + SP_Balloon_mass)
    total_masses.append(total_mass)

#print("Final total mass: ", total_mass)
gondola_volume = gondola_mass / gondola_density
gondola_size = gondola_volume **(1/3)

max_radius = (V_zp * 3 / (4 * np.pi)) ** (1/3)
max_diameter = max_radius * 2


print("\nVarying altitude balloon sizing for a paylad mass of",payload_mass,"kg and an altitude ranging from",H1,"km to",H2,"km:")
print("assuming a gondola mass of", gondola_mass, "kg, and an inflation system mass of", inflation_mass, "kg, a super pressure balloon volume of", V_SP, "m^3 and a maximum super pressure of", DP1, "Pa:")
print("\ngondola mass: ", gondola_mass)
print("\nsuper pressure balloon mass: ", round(SP_Balloon_mass,2),"kg")
print("zero pressure balloon mass: ", round(ZP_balloon_mass,2),"kg")
print("super pressure balloon diameter:", round(r_sp*2,2),"m")
print("zero pressure balloon diameter:", round(r_zp*2,2),"m")
print("total helium moles: ", round(n_sp1 + n_zp1,2),"mol")
print("helium mass: ", round(He_mass,2),"kg")
print("Super pressure at ", H1, "km: ", DP1,"Pa")
print("Super pressure at ", H2, "km: ", round(DP2, 2),"Pa")
print("moles of gas in super pressure balloon at ", H1, "km: ", round(n_sp1, 2),"mol")
print("moles of gas in super pressure balloon at ", H2, "km: ", round(n_sp2, 2),"mol")
print("zero pressure balloon volume at",H1,"km : ", round(V_zp1, 2),"m^3")
print("zero pressure balloon volume at",H2,"km : ", round(V_zp, 2),"m^3")
print("\nFinal total mass: ", total_mass, "kg")
print("\n")

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(iterations, ZP_balloon_masses, label="Zero Pressure Balloon Mass (kg)")
plt.plot(iterations, SP_balloon_masses, label="Super Pressure Balloon Mass (kg)")
plt.plot(iterations, balloon_masses, label="Total Balloon Mass (kg)")
plt.plot(iterations, total_masses, label="Total Mass (kg)")
plt.xlabel("Iteration")
plt.ylabel("Mass (kg)")
plt.title("Balloon and Total Mass vs Iteration")
plt.legend()
plt.grid()
plt.show()




