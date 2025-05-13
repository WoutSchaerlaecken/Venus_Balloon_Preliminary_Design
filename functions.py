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


def SP_balloon_pressure(n, R, T, V):
    P = (n * R * T) / V
    return P


def plot_balloons(super_pressure_volume, zero_pressure_volume, zero_pressure_1_volume):
    def sphere_radius(volume):
        return ((3 * volume) / (4 * 3.14)) ** (1/3)
    
    super_pressure_radius = sphere_radius(super_pressure_volume)
    zero_pressure_radius = sphere_radius(zero_pressure_volume)
    zero_pressure_1_radius = sphere_radius(zero_pressure_1_volume)
    
    fig, ax = plt.subplots()
    
    # Draw the circles
    super_pressure_circle = plt.Circle((0, 0), super_pressure_radius, color='blue', alpha=0.5, label='Super Pressure Balloon')
    zero_pressure_circle = plt.Circle((0, zero_pressure_radius-super_pressure_radius), zero_pressure_radius, color='orange', alpha=0.5, label='Zero Pressure Balloon')
    zero_pressure_1_circle = plt.Circle((0, zero_pressure_1_radius-super_pressure_radius), zero_pressure_1_radius, color='red', alpha=0.5, label='Zero Pressure Balloon 1')

    ax.add_artist(super_pressure_circle)
    ax.add_artist(zero_pressure_circle)
    ax.add_artist(zero_pressure_1_circle)
    
    # Set limits to fit the balloons
    ax.set_xlim(-zero_pressure_radius * 1.3,zero_pressure_radius * 1.3)
    ax.set_ylim(-super_pressure_radius*1.3, zero_pressure_radius * 2)
    
    ax.set_aspect('equal')
    plt.legend()
    plt.title('Visual Representation of Balloons')
    plt.show()

class SpacecraftBudgets:
    def __init__(self, payload_percentage, structures_percentage, thermal_percentage, power_percentage, ttc_percentage, processing_percentage, adcs_percentage, other_percentage, payload_power, power_power, ttc_power, processing_power):
        self.payload_percentage = payload_percentage / 100
        self.structures_percentage = structures_percentage / 100
        self.thermal_percentage = thermal_percentage / 100
        self.power_percentage = power_percentage / 100
        self.ttc_percentage = ttc_percentage / 100
        self.processing_percentage = processing_percentage / 100
        self.adcs_percentage = adcs_percentage / 100
        self.other_percentage = other_percentage / 100
        self.payload_power = payload_power / 100
        self.power_power = power_power / 100
        self.ttc_power = ttc_power / 100
        self.processing_power = processing_power / 100



    def summary(self):
        print(f"Payload Percentage       : {self.payload_percentage:.1f}%")
        print(f"Structures Percentage    : {self.structures_percentage:.1f}%")
        print(f"Thermal Percentage       : {self.thermal_percentage:.1f}%")
        print(f"Power Percentage         : {self.power_percentage:.1f}%")
        print(f"TTC Percentage           : {self.ttc_percentage:.1f}%")
        print(f"Processing Percentage    : {self.processing_percentage:.1f}%")
        print(f"ADCS Percentage          : {self.adcs_percentage:.1f}%")
        print(f"Other Percentage         : {self.other_percentage:.1f}%")
        print(f"Payload Power       : {self.payload_power:.1f}%")
        print(f"Structures Power    : {self.structures_power:.1f}%")
        print(f"TTC Power           : {self.ttc_power:.1f}%")
        print(f"Processing Power    : {self.processing_power:.1f}%")



def gondola_vs_total_mass(gondola_mass):
    tot_mass = gondola_mass * 1.4497 + 15.147
    return tot_mass

def simple_balloon_volume(p, Ma, R, T, Mg, mass):
    V = (mass * R * T) / (p*(Ma-Mg))
    return V


