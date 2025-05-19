from functions import *
from inputs import *
import numpy as np
from constants import *
import matplotlib.pyplot as plt
from Venus_atmosphere import venus_atmosphere
from Power_Balloon import Pd

H1 = 40
H2 = 57


T1, p1,rho = venus_atmosphere(H1)
T2, p2,_ = venus_atmosphere(H2)


print(T2)

balloon_ratio = 1/3

He_permeability = 700 #cm^3/m^2/day for wrinkled at 50C

inflation_mass = 2
sample_collection_mass = 1

power_mass = 2


probe_mass = payload_mass / probe.payload_percentage 
print(probe_mass*probe.power_percentage)
print("probe mass: ", probe_mass)
structures_mass = probe_mass * probe.structures_percentage
ttc_mass = probe_mass * probe.ttc_percentage
thermal_mass = probe_mass * probe.thermal_percentage
processing_mass = probe_mass * probe.processing_percentage
other_mass = probe_mass * probe.other_percentage



# Print subsystem masses
print("Subsystem Masses:")
print(f"{'Payload Mass:':<30} {payload_mass:.2f} kg")
print(f"{'Structures Mass:':<30} {structures_mass:.2f} kg")
print(f"{'Thermal Mass:':<30} {thermal_mass:.2f} kg")
print(f"{'Power Mass:':<30} {power_mass:.2f} kg")
print(f"{'TTC Mass:':<30} {ttc_mass:.2f} kg")
print(f"{'Processing Mass:':<30} {processing_mass:.2f} kg")
print(f"{'Other Mass:':<30} {other_mass:.2f} kg")
print(f"{'Sample Collection Mass:':<30} {sample_collection_mass:.2f} kg")
print(f"{'Inflation Mass:':<30} {inflation_mass:.2f} kg")



print(f"{'Payload Power:':<30} {probe.payload_power * Pd:.2f} W")
print(f"{'Power Power:':<30} {probe.power_power * Pd:.2f} W")
print(f"{'TTC Power:':<30} {probe.ttc_power * Pd:.2f} W")
print(f"{'Processing Power:':<30} {probe.processing_power * Pd:.2f} W")
#probe_mass = probe_mass + power_mass - probe_mass * probe.power_percentage


gondola_mass = payload_mass + structures_mass + ttc_mass + thermal_mass + processing_mass + other_mass + sample_collection_mass + inflation_mass + power_mass

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

    He_mass = (n_sp2 + n_zp2)*Mhe

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

B = (n_zp2-n_lost*6)*R/R_Venus + p1/(R_Venus*T1)*V_SP

days = 1/n_lost * (n_zp2 + p1/(R*T1) * V_SP - total_mass * R_Venus/R)



"""---------------Calculating Energy Required for Altitude Control------------------"""
m1 = n_sp2 * Mhe 

m2 = n_sp1 * Mhe

a = (p2-p1)/(m2-m1)

b = p1 - a*m1

a1 = (DP2-DP1)/(m2-m1)
b1 = DP1 - a1*m1

import numpy as np
from scipy.integrate import quad

# Define the function to integrate
def integrand(x):
    numerator = a * x + b - a1 * x + b1
    denominator = a * x + b
    power_term = (numerator / denominator) ** (gamma_He-1 / gamma_He)
    return Cp_He * T1 * power_term

# Set integration bounds
lower_bound = m1  # adjust as needed
upper_bound = m2  # adjust as needed

# Perform the numerical integration
E, error = quad(integrand, lower_bound, upper_bound)

#print(f"Estimated error: {error:.2e}")




print("Energy required for altitude control: ", E)




# Create a breakdown of the masses for a stacked bar chart
import matplotlib.pyplot as plt

# Mass breakdowns
gondola_labels = [
    "Payload", "Structures", "Thermal", "Power", "Communications",
    "Command & Data Handling", "Sampling System"
]
gondola_values = [
    payload_mass, structures_mass, thermal_mass, power_mass, ttc_mass,
    processing_mass, sample_collection_mass
]

balloon_labels = [
    "Pressure system", "Super Pressure Balloon", "Zero Pressure Balloon", "Helium"
]
balloon_values = [
    inflation_mass, SP_Balloon_mass, ZP_balloon_mass, He_mass
]

# System-level aggregation
system_labels = ["Gondola", "Balloon"]
system_values = [sum(gondola_values), sum(balloon_values)]

# Plot
fig, ax = plt.subplots(figsize=(10, 8))

# X positions for bars
bar_positions = [0, 1, 2]  # Gondola, Balloon, System
bar_width = 0.6

# Plot Gondola breakdown
bottom = 0
for label, value in zip(gondola_labels, gondola_values):
    ax.bar(bar_positions[0], value, bottom=bottom, label=label if bar_positions[0] == 0 else "")
    bottom += value

# Plot Balloon breakdown
bottom = 0
for label, value in zip(balloon_labels, balloon_values):
    ax.bar(bar_positions[1], value, bottom=bottom, label=label if bar_positions[1] == 1 else "")
    bottom += value

# Plot System breakdown (just two bars stacked)
bottom = 0
for label, value in zip(system_labels, system_values):
    ax.bar(bar_positions[2], value, bottom=bottom, label=label if bar_positions[2] == 2 else "")
    bottom += value

# X-axis labels and formatting
ax.set_xticks(bar_positions)
ax.set_xticklabels(["Gondola", "Balloon", "System"])
ax.set_ylabel("Mass (kg)")
ax.set_title("Mass Breakdown of Gondola, Balloon and System for a Payload Mass of " + str(payload_mass) + " kg and an Altitude Range of " + str(H1) + " km to " + str(H2) + " km")
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1))
plt.tight_layout()




print("\nVarying altitude balloon sizing for a payload mass of", payload_mass, "kg and an altitude ranging from", H1, "km to", H2, "km:")
print("Assuming a gondola mass of", gondola_mass, "kg, and an inflation system mass of", inflation_mass, "kg,")
print("a super pressure balloon volume of", round(V_SP, 2), "m^3 and a maximum super pressure of", DP1, "Pa:\n")

print(f"{'Gondola mass:':<40} {gondola_mass:.2f} kg")
print(f"{'Gondola size:':<40} {gondola_size:.2f} m\n")
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
#plt.show()


balloon_center = np.array([0, 0, r_zp * 1.5])  # Balloon center at (0, 0, balloon_radius * 1.5)
gondola_center = np.array([0, 0, -gondola_size / 2])  # Gondola center at (0, 0, -gondola_size / 2)

# Create the figure and 3D axis
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the balloon as a sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = r_zp * np.outer(np.cos(u), np.sin(v)) + balloon_center[0]
y = r_zp * np.outer(np.sin(u), np.sin(v)) + balloon_center[1]
z = r_zp * np.outer(np.ones(np.size(u)), np.cos(v)) + balloon_center[2]
ax.plot_surface(x, y, z, color='blue', alpha=0.5, label="Balloon")

# Plot the gondola as a cube
gondola_size_half = gondola_size / 2
gondola_vertices = np.array([
    [-gondola_size_half, -gondola_size_half, -gondola_size_half],
    [gondola_size_half, -gondola_size_half, -gondola_size_half],
    [gondola_size_half, gondola_size_half, -gondola_size_half],
    [-gondola_size_half, gondola_size_half, -gondola_size_half],
    [-gondola_size_half, -gondola_size_half, gondola_size_half],
    [gondola_size_half, -gondola_size_half, gondola_size_half],
    [gondola_size_half, gondola_size_half, gondola_size_half],
    [-gondola_size_half, gondola_size_half, gondola_size_half],
]) + gondola_center

# Define the edges of the cube
gondola_edges = [
    [0, 1], [1, 2], [2, 3], [3, 0],  # Bottom face
    [4, 5], [5, 6], [6, 7], [7, 4],  # Top face
    [0, 4], [1, 5], [2, 6], [3, 7]   # Vertical edges
]

# Plot the gondola edges
for edge in gondola_edges:
    ax.plot3D(*zip(*gondola_vertices[edge]), color='red', alpha=0.7)

# Draw a line connecting the balloon and gondola
ax.plot(
    [balloon_center[0], gondola_center[0]],
    [balloon_center[1], gondola_center[1]],
    [balloon_center[2] - r_zp, gondola_center[2] + gondola_size_half],
    color='black', linestyle='--'
)

# Set the aspect ratio and limits
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
ax.set_xlim(-r_zp * 1.5, r_zp * 1.5)
ax.set_ylim(-r_zp * 1.5, r_zp * 1.5)
ax.set_zlim(-gondola_size * 2, r_zp * 3.5)

# Add labels and legend
ax.set_title("3D Balloon and Gondola Visualization")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")



# 2D Visualization (side view)
fig2, ax2 = plt.subplots(figsize=(8, 8))

# Draw balloon as a circle
balloon_circle = plt.Circle((0, r_zp * 1.5), r_zp, color='blue', alpha=0.5, label='Balloon')
ax2.add_patch(balloon_circle)

# Draw gondola as a square
gondola_bottom = -gondola_size / 2
gondola_left = -gondola_size / 2
gondola_square = plt.Rectangle(
    (gondola_left, gondola_bottom), gondola_size, gondola_size,
    color='red', alpha=0.7, label='Gondola'
)
ax2.add_patch(gondola_square)

# Draw line connecting balloon and gondola
ax2.plot(
    [0, 0],
    [r_zp * 1.5 - r_zp, gondola_bottom + gondola_size],
    color='black', linestyle='--'
)

# Set limits and labels
ax2.set_xlim(-r_zp * 1.5, r_zp * 1.5)
ax2.set_ylim(-gondola_size * 2, r_zp * 3.5)
ax2.set_aspect('equal')
ax2.set_xlabel("X (m)")
ax2.set_ylabel("Z (m)")
ax2.set_title("2D Side View: Balloon and Gondola")
ax2.legend(loc="upper right", fontsize=25)
plt.tight_layout()
plt.show()
