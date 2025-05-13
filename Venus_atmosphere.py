import os 
import numpy as np 

#---------------Venus atmosphere--------------
#From 0 - 250km
#---------------------------------------------

# Add altitude in km 

def venus_atmosphere(alt):
    current_dir = os.path.dirname(os.path.abspath(__file__))

    file_dens_low = os.path.join(current_dir, "VIRALow.txt")
    file_dens_mid = os.path.join(current_dir, "VIRAMid.txt")
    file_dens_high = os.path.join(current_dir, "VIRAHi.txt")

    altitudes = []
    densities = []
    Temp = []
    pressure = [] #Nm^2 

    with open(file_dens_low, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and parts[0].replace('.', '', 1).isdigit():
                try:
                    z_km = float(parts[0])
                    rho = float(parts[2])
                    pres = float(parts[3])
                    temp = float(parts[4])
                    altitudes.append(z_km)
                    densities.append(rho)
                    Temp.append(temp)
                    pressure.append(pres)
                except ValueError:
                    continue  

    altitudes.pop()
    densities.pop()
    Temp.pop()
    pressure.pop()

    with open(file_dens_mid, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and parts[0].replace('.', '', 1).isdigit():
                try:
                    z_km = float(parts[0])
                    rho = float(parts[8])
                    pres = float(parts[9])
                    temp = float(parts[10])
                    altitudes.append(z_km)
                    densities.append(rho)
                    Temp.append(temp)
                    pressure.append(pres)
                except ValueError:
                    continue  

    altitudes.pop()
    densities.pop()
    Temp.pop()
    pressure.pop()

    with open(file_dens_high, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and parts[0].replace('.', '', 1).isdigit():
                try:
                    z_km = float(parts[0])
                    rho = float(parts[9])
                    pres = float(parts[10])
                    temp = float(parts[11])
                    altitudes.append(z_km)
                    densities.append(rho)
                    Temp.append(temp)
                    pressure.append(pres)
                except ValueError:
                    continue 

    temp_alt = float(np.interp(alt, altitudes, Temp))
    pres_alt = float(np.interp(alt, altitudes, pressure))
    dens_alt = float(np.interp(alt, altitudes, densities))

    return temp_alt, pres_alt, dens_alt

