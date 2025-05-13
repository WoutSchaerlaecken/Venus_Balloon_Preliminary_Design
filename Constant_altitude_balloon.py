from functions import *
from inputs import *
import numpy as np
from constants import *
from Power_Balloon import total_power_mass
from mpl_toolkits.mplot3d import Axes3D



gondola_mass = payload_mass / probe.payload_percentage

structures_mass = gondola_mass * probe.structures_percentage
thermal_mass = gondola_mass * probe.thermal_percentage
power_mass = gondola_mass * probe.power_percentage
ttc_mass = gondola_mass * probe.ttc_percentage
processing_mass = gondola_mass * probe.processing_percentage
adcs_mass = gondola_mass * probe.adcs_percentage
other_mass = gondola_mass * probe.other_percentage


gondola_mass = gondola_mass + 1 
print("gondola mass: ", gondola_mass)
total_mass = gondola_vs_total_mass(gondola_mass) 
print("total mass: ", total_mass)



Balloon_volume = simple_balloon_volume(p3, Ma, R, T3, Mhe, total_mass)


balloon_mass = Balloon_mass(Balloon_volume, balloon_material_density)



iteration = [0]
balloon_masses = [balloon_mass]
total_masses = [total_mass]

for i in range(20):
    Balloon_volume = simple_balloon_volume(p3, Ma, R, T3, Mhe, total_mass)

    balloon_mass = Balloon_mass(Balloon_volume, balloon_material_density)

    total_mass = gondola_mass + balloon_mass

    iteration.append(i+1)
    balloon_masses.append(balloon_mass)
    total_masses.append(total_mass)

gondola_volume = gondola_mass / gondola_density
gondola_size = gondola_volume **(1/3)

balloon_radius = (Balloon_volume * 3 / (4 * np.pi)) ** (1/3)
print("Balloon volume: ", Balloon_volume)
print("Balloon radius: ", balloon_radius)

print("Gondola size: ", gondola_size)


print('\n For a payload mass of ', payload_mass, 'kg:')
print("Final total mass: ", total_mass)

print("\n Mass Budget Breakdown:")
print("Final balloon mass: ", balloon_mass)
print("Final gondola mass: ", gondola_mass)
print("Final structures mass: ", structures_mass)
print("Final thermal mass: ", thermal_mass)
print("Final power mass: ", power_mass)
print("Final ttc mass: ", ttc_mass)
print("Final processing mass: ", processing_mass)
print("Final adcs mass: ", adcs_mass)
print("Final other mass: ", other_mass)


plt.figure(figsize=(10, 6))
plt.plot(iteration, balloon_masses, label="Balloon Mass (kg)", marker='o')
plt.plot(iteration, total_masses, label="Total Mass (kg)", marker='x')
plt.xlabel("Iteration")
plt.ylabel("Mass (kg)")
plt.title("Balloon and Total Mass vs Iteration")
plt.legend()
plt.grid()
plt.show()



import matplotlib.pyplot as plt

# Define the balloon and gondola dimensions
balloon_center = np.array([0, 0, balloon_radius * 1.5])  # Balloon center at (0, 0, balloon_radius * 1.5)
gondola_center = np.array([0, 0, -gondola_size / 2])  # Gondola center at (0, 0, -gondola_size / 2)

# Create the figure and 3D axis
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the balloon as a sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = balloon_radius * np.outer(np.cos(u), np.sin(v)) + balloon_center[0]
y = balloon_radius * np.outer(np.sin(u), np.sin(v)) + balloon_center[1]
z = balloon_radius * np.outer(np.ones(np.size(u)), np.cos(v)) + balloon_center[2]
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
    [balloon_center[2] - balloon_radius, gondola_center[2] + gondola_size_half],
    color='black', linestyle='--'
)

# Set the aspect ratio and limits
ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
ax.set_xlim(-balloon_radius * 1.5, balloon_radius * 1.5)
ax.set_ylim(-balloon_radius * 1.5, balloon_radius * 1.5)
ax.set_zlim(-gondola_size * 2, balloon_radius * 3.5)

# Add labels and legend
ax.set_title("3D Balloon and Gondola Visualization")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
plt.show()