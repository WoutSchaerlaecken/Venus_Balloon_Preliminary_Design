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

He_permeability = 700 #cm^3/m^2/day for wrinkled at 50C

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

#calculate the duration based on Helium leaking


zp_2_surface = r_zp**2 * 4 * np.pi
n_lost = p1*zp_2_surface * He_permeability * 1e-6 / (R * T1)
a = zp_2_surface * He_permeability * 1e-6 

B = (n_zp2-n_lost*6)*R/R_Venus + p1/(R_Venus*T1)*V_SP

days = 1/n_lost * (n_zp2 + p1/(R*T1) * V_SP - total_mass * R_Venus/R)






print("\nVarying altitude balloon sizing for a payload mass of", payload_mass, "kg and an altitude ranging from", H1, "km to", H2, "km:")
print("Assuming a gondola mass of", gondola_mass, "kg, and an inflation system mass of", inflation_mass, "kg,")
print("a super pressure balloon volume of", round(V_SP, 2), "m^3 and a maximum super pressure of", DP1, "Pa:\n")

print(f"{'Gondola mass:':<40} {gondola_mass:.2f} kg")
print(f"{'Super pressure balloon mass:':<40} {SP_Balloon_mass:.2f} kg")
print(f"{'Zero pressure balloon mass:':<40} {ZP_balloon_mass:.2f} kg\n")
print(f"{'Super pressure balloon diameter:':<40} {r_sp*2:.2f} m")
print(f"{'Zero pressure balloon diameter:':<40} {r_zp*2:.2f} m\n")
print(f"{'Total helium moles:':<40} {n_sp1 + n_zp1:.2f} mol")
print(f"{'Helium mass:':<40} {He_mass:.2f} kg\n")
print(f"{'Super pressure at ' + str(H1) + ' km:':<40} {DP1:.2f} Pa")
print(f"{'Super pressure at ' + str(H2) + ' km:':<40} {DP2:.2f} Pa\n")
print(f"{'Moles in SP balloon at ' + str(H1) + ' km:':<40} {n_sp1:.2f} mol")
print(f"{'Moles in SP balloon at ' + str(H2) + ' km:':<40} {n_sp2:.2f} mol\n")
print(f"{'ZP balloon volume at ' + str(H1) + ' km:':<40} {V_zp1:.2f} m^3")
print(f"{'ZP balloon volume at ' + str(H2) + ' km:':<40} {V_zp:.2f} m^3\n")
print(f"{'Final total mass:':<40} {total_mass:.2f} kg")
print(f"{'The balloon will stay within the altitude range of ' + str(H1) + ' km to ' + str(H2) + ' km for:':<40} {days:.2f} days")


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




