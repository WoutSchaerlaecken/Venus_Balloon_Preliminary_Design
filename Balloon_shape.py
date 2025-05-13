import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# Constants
g = 9.80665  # m/s^2
w = 0.1  # kg/m^2, areal density of the balloon film
Rz = 10.0  # m, outer envelope radius
Lg = 20.0  # m, gore length
payload_mass = 50.0  # kg
sigma_c = 100.0  # N/m (assumed)
T = 273.15 + 20  # K
gas_constant = 8.314  # J/(mol·K)
M_gas = 4.0026e-3  # kg/mol (Helium)
M_air = 28.97e-3  # kg/mol (Air)

# Sample atmospheric and gas densities
rho_atm = 1.2  # kg/m^3
rho_gas = 0.1785  # kg/m^3 (He at sea level)

# Initial conditions
def initial_conditions():
    rg0 = 1.0  # initial radius guess
    theta0 = np.pi / 4  # 45 degrees
    sigma_m0 = 50.0  # N/m (assumed)
    A0 = 0.0
    V0 = 0.0
    z0 = 0.0
    return [rg0, theta0, sigma_m0, A0, V0, z0]

# System of ODEs
def shape_odes(s, y, rho_atm, rho_gas, z0):
    rg, theta, sigma_m, A, V, z = y
    b = g * (rho_atm - rho_gas)
    
    dtheta_ds = (1 / (rg * sigma_m)) * (
        sigma_c * np.cos(theta) - rg * w * np.sin(theta) - rg * b * (z - z0)
    )
    dsigma_m_ds = (1 / rg) * (
        sigma_c * np.sin(theta) - rg * w * np.cos(theta) - sigma_m * np.sin(theta)
    )
    drg_ds = np.sin(theta)
    dz_ds = np.cos(theta)
    dA_ds = 2 * np.pi * rg
    dV_ds = np.pi * rg**2 * np.cos(theta)

    return [drg_ds, dtheta_ds, dsigma_m_ds, dA_ds, dV_ds, dz_ds]

# Shooting method objective
def shooting_objective(z0_guess, rho_atm, rho_gas):
    y0 = initial_conditions()
    y0[-1] = z0_guess
    sol = solve_ivp(
        shape_odes,
        [0, Lg],
        y0,
        args=(rho_atm, rho_gas, z0_guess),
        dense_output=False
    )
    final_z = sol.y[5, -1]  # final z
    z_target = 15.0  # desired end height (example)
    return final_z - z_target

# Solve for correct z0
def solve_shape():
    res = root_scalar(
        shooting_objective,
        args=(rho_atm, rho_gas),
        bracket=[-10, 10],
        method='brentq'
    )
    z0_final = res.root
    y0 = initial_conditions()
    y0[-1] = z0_final
    sol = solve_ivp(
        shape_odes,
        [0, Lg],
        y0,
        args=(rho_atm, rho_gas, z0_final),
        dense_output=False
    )
    return sol

# Run and plot
solution = solve_shape()
rg_vals, z_vals = solution.y[0], solution.y[5]

plt.plot(rg_vals, z_vals)
plt.xlabel("Radius rg (m)")
plt.ylabel("Height z (m)")
plt.title("Balloon Midsection Shape")
plt.grid(True)
plt.axis('equal')
plt.show()


import matplotlib.pyplot as plt
import numpy as np

def plot_multiple_balloon_shapes(rg_solutions, z_solutions, R_sp, R_zp, z_shift):
    """
    rg_solutions: list of np.arrays for radius profiles
    z_solutions: list of np.arrays for z profiles (same length as rg_solutions)
    R_sp, R_zp: radii of the bottom and top spherical caps
    z_shift: vertical distance from SP center to base of midsection
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot balloon shapes
    for rg, z in zip(rg_solutions, z_solutions):
        ax.plot(rg, z, linewidth=1)
        ax.plot(-rg, z, linewidth=1)

    # Add SP bottom cap (circle)
    theta_sp = np.linspace(np.pi, 2*np.pi, 200)
    r_sp = R_sp * np.cos(theta_sp)
    z_sp = R_sp * np.sin(theta_sp) + z_shift
    ax.plot(r_sp, z_sp, 'k--', label="Super-Pressure Balloon")
    ax.plot(-r_sp, z_sp, 'k--')

    # Add ZP top cap (circle)
    theta_zp = np.linspace(0, np.pi, 200)
    r_zp = R_zp * np.cos(theta_zp)
    z_zp = R_zp * np.sin(theta_zp) + max([z[-1] for z in z_solutions])
    ax.plot(r_zp, z_zp, 'k-.', label="Zero-Pressure Balloon")
    ax.plot(-r_zp, z_zp, 'k-.')

    # Optional red dot (e.g., payload location)
    ax.plot(rg_solutions[0][-1], z_solutions[0][-1], 'ro')

    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.axis("equal")
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()

# Example usage
# You would have to generate rg_solutions and z_solutions for multiple volumes
# Generate example radius and height profiles for multiple balloon shapes
rg_solutions = [np.linspace(0, 5 + i, 100) for i in range(100)]
z_solutions = [np.linspace(0, 10 + 2 * i, 100) for i in range(100)]

# Then call:
plot_multiple_balloon_shapes(rg_solutions, z_solutions, R_sp=1.0, R_zp=1.0, z_shift=1.0)

