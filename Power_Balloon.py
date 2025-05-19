from functions import *
from inputs import *
import numpy as np
from constants import *

"""just realised that this is actually for earth, so should be better on Venus"""
t_round = 6 #days
Pd = 25 #W
Pe = 1 #W
Psp = 75 #specific power in W/kg
P_delta = 5/1000 #area per power in m^2/W 
E_sp = 130 #specific power W/kg


td = 24*3600*t_round * 0.5
te = 24*3600*t_round * 0.5

eta_d = 0.8
eta_e = 0.6
eta_bat = 0.9
DOD = 0.7

SC = 1000 #W/m^2 conservative estimate for around 50km
eta_1 = 0.18 #conversion efficiency of solar cells
eta_2 = 0.8 #packing density of solar cells
rho = 2  #kg/m^2 solar array mass density

gondola_mass = payload_mass / probe.payload_percentage
power_mass = gondola_mass * probe.power_percentage

P_sa = 1/td *(Pd*td/eta_d + Pe*te/eta_e)
print("P_sa: ", P_sa)

array_area = P_sa/(eta_1 * eta_2 * SC)
#array_area = P_sa/300
array_mass = rho * array_area

#array_area = P_delta * P_sa
#array_mass = P_sa/Psp

E_bat = (Pe*te/3600)/(eta_bat*DOD)

M_bat = E_bat/E_sp


print("\nOriginal power mass: ", power_mass)

print("\nSolar array properties:")
print("Solar array area: ", array_area)
print("Solar array mass: ", array_mass)
print("\nBattery mass: ", M_bat)

total_power_mass = (M_bat + array_mass)*3/2
M_pmd = 1/3 * total_power_mass

print("\nPower mass distribution: ", M_pmd)
print("\nUpdated total power mass: ", total_power_mass)

if __name__ == "__main__":
    # Execute the main functionality of the script
    pass


