def Balloon_mass(Volume, balloon_mat_desity):
    radius = ((3 * Volume) / (4 * 3.14)) ** (1/3)
    surface_area = 4 * 3.14 * radius ** 2
    mass = balloon_mat_desity * surface_area
    return mass