from functions import *
from inputs import *
import numpy as np
from constants import *
from Power_Balloon import total_power_mass
from mpl_toolkits.mplot3d import Axes3D
from Power_Balloon import Pd

power = Pd * probe.ttc_power
print("Available Communications Power: ", power)

gondola_mass = payload_mass / probe.payload_percentage
power_mass = gondola_mass * probe.power_percentage

print("Power mass: ", power_mass)