from functions import *
from inputs import *
import numpy as np
from constants import *


n_zp_2 = Zero_pressure_moles(total_mass, T2 , V_SP, p2 ,R_Venus , R)
print("Zero pressure moles: ", n_zp_2)

V_zp = ZP_balloon_volume(n_zp_2, R, T2, p2)
print("Zero pressure volume: ", V_zp)

balloon_mass = Balloon_mass(V_zp, balloon_material_density)
print("Balloon mass: ", balloon_mass)


n_zp_1 = N_zp_min(n_zp_2, p1, T1, p2, T2, R, V_SP)
print("Minimum Zero pressure moles: ", n_zp_1)

V_zp_min = ZP_balloon_volume(n_zp_1, R, T2, p2)
print("Minimum Zero pressure volume: ", V_zp_min)

plot_balloons(V_SP, V_zp)
plot_balloons(V_zp, V_zp_min)

x = (p1/T1-p2/T2) * (V_SP / R)
print("X: ", x)