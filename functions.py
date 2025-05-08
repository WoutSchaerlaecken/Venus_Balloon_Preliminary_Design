import matplotlib.pyplot as plt

def Balloon_mass(Volume, balloon_mat_desity):
    radius = ((3 * Volume) / (4 * 3.14)) ** (1/3)
    surface_area = 4 * 3.14 * radius ** 2
    mass = balloon_mat_desity * surface_area
    return mass


def Zero_pressure_moles(mass, temp, Vsp, pressure, Ra, R):
    n_zp = (Ra/R)*(mass-(pressure/(Ra*temp)*Vsp))
    return n_zp

def ZP_balloon_volume(n, R, T, P):
    V = (n * R * T) / P
    return V

def N_zp_min(n_zp_max, P1, T1, P2, T2, R, V_sp):
    n_zp_min = n_zp_max - ((P1/T1 - P2/T2) * (V_sp / R))
    return n_zp_min




def plot_balloons(super_pressure_volume, zero_pressure_volume):
    def sphere_radius(volume):
        return ((3 * volume) / (4 * 3.14)) ** (1/3)
    
    super_pressure_radius = sphere_radius(super_pressure_volume)
    zero_pressure_radius = sphere_radius(zero_pressure_volume)
    
    fig, ax = plt.subplots()
    
    # Draw the circles
    super_pressure_circle = plt.Circle((0, 0), super_pressure_radius, color='blue', alpha=0.5, label='Super Pressure Balloon')
    zero_pressure_circle = plt.Circle((0, zero_pressure_radius-super_pressure_radius), zero_pressure_radius, color='orange', alpha=0.5, label='Zero Pressure Balloon')
    
    ax.add_artist(super_pressure_circle)
    ax.add_artist(zero_pressure_circle)
    
    # Set limits to fit the balloons
    ax.set_xlim(-zero_pressure_radius * 1.3,zero_pressure_radius * 1.3)
    ax.set_ylim(-super_pressure_radius*1.3, zero_pressure_radius * 2)
    
    ax.set_aspect('equal')
    plt.legend()
    plt.title('Visual Representation of Balloons')
    plt.show()