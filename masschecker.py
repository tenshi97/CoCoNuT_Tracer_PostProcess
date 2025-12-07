import h5reader
import numpy as np
h = h5reader.hydrof(model="20m_new",subdirectory="models")
r_len = len(h.xzn())
theta_len = len(h.yzn())
M_ALL = 0
c = 299792458e2
M_S = 1.989e33
h.set_index(h.ngroups-1)
for x in range(r_len):
    for y in range(theta_len):
        r_1 = h.xzl()[x]
        r_2 = h.xzr()[x]
        r_mid = h.xzn()[x]
        theta_1 = h.yzl()[y]
        theta_2 = h.yzr()[y]
        theta_mid = h.yzn()[y]
        DV = 2 * np.pi * (np.cos(theta_1) - np.cos(theta_2)) * (r_2 ** 3 - r_1 ** 3) / 3
        # DV = 2*np.pi*np.sin(theta_mid)*r_mid * 0.5 *(theta_2-theta_1)*(r_2**2-r_1**2)
        LF = 1.0 / np.sqrt(1 - (h.vex()[x][y][0] ** 2 + h.vey()[x][y][0] ** 2) / (c ** 2))  # Lorentz Factor
        DM = h.den()[x][y][0] * DV * (h.phi()[x][y][0]**6) * LF
        M_ALL += DM
print(M_ALL/M_S)