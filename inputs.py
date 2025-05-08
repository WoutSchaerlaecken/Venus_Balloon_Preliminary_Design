H1 = 40000
T1 = 416.15
p1 = 354738

H2 = 62000
T2 = 254.15
p2 = 17000

V_SP = 4


total_mass = 25
balloon_material_density = 0.25  # kg/m^3

probe_fraction = 0.4
satellite_fraction = 1-probe_fraction
satellite_mass = total_mass * satellite_fraction
probe_mass = total_mass * probe_fraction

Gondola_mass = 0.5 * total_mass  # kg
