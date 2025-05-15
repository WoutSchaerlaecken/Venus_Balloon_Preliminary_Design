from functions import *
H1 = 40000
T1 = 416.15
p1 = 354738

H2 = 62000
T2 = 254.15
p2 = 17000

H3 = 53000
T3 = 323.0
p3 = 71162




V_SP = 5

Mhe = 4.0026 /1000 # kg/mol
Ma = 43.449175/1000 # kg/mol


payload_mass = 1

#total_mass = 25
balloon_material_density = 0.2  # kg/m^2
SP_balloon_material_density = 0.25  # kg/m^2
gondola_density = 700  # kg/m^3



probe = SpacecraftBudgets(
    payload_percentage=23.8,
    structures_percentage=17.4,
    thermal_percentage=6.5,
    power_percentage=28.2,
    ttc_percentage=13.5,
    processing_percentage=10.5,
    adcs_percentage = 0,
    other_percentage= 0,
    payload_power=50,
    power_power=5,
    ttc_power=25,
    processing_power=20,
)


