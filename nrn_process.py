
from pathlib import Path
import os

import numpy
import numpy as np
import matplotlib.pyplot as plt
import phycon_and_nuc_table
import re
folder_pathname = "F:/compastro/Tracers16000/nucelosynthesis/CCSN_14m_old/"
folder_path = Path(folder_pathname)

subfolders = [f.name for f in folder_path.iterdir() if f.is_dir()]
count_bad = 0
mass_assign_file = open("20m_new_mass_assign.txt","r")
mass_assign = []
M_S = 1.989e33
total_mass = 0
inspect = phycon_and_nuc_table.pc_nuc_inspect_all
time_points = np.linspace(0,1.5,100)
mass_inspect_time = [[] for i in range(len(inspect))]
for i in range(len(mass_inspect_time)):
    mass_inspect_time[i] = [0 for x in range(len(time_points))]
for line in mass_assign_file:
    if(": " not in line):
        continue
    print(line)
    mass = float(line.split(": ")[1])
    mass_assign.append(mass)
total_weighed_result = [0 for i in range(195)]
MASSES = [np.zeros(len(time_points)) for i in range(len(inspect))]
INIT = []
for subfolder in subfolders:
    id = int(subfolder.split("_")[3])
    mass = mass_assign[int(id)-1]
    total_mass += mass
    winnet_tracer = folder_pathname + subfolder + "/"
    snapshot_folder = folder_pathname + subfolder +"/snaps/"
    if(not os.path.exists(winnet_tracer+"finabsum.dat")):
        count_bad += 1
        continue
    nrn_result = np.loadtxt(winnet_tracer+"finabsum.dat",skiprows=1)
    snapdir = Path(snapshot_folder)
    Times = []
    XNUS = [[] for i in range(len(inspect))]
    INIT_MASS = [0 for i in range(200)]
    for file in snapdir.iterdir():
        if(file.is_file()):
            with file.open("r") as snap:
                fs = snap.readline()
                fs = snap.readline()
                time = float(fs.split()[0])
                xnus = np.loadtxt(file,skiprows=3)
                Times.append(time)
                if(time==0):
                    for i in range(len(inspect)):
                        for ele in xnus:
                            neutron = int(ele[0])
                            proton = int(ele[1])
                            xnu = float(ele[3])
                            if(neutron+proton<=100 and xnu>1e-50):
                                INIT_MASS[neutron+proton] += mass * xnu / M_S
                for i in range(len(inspect)):
                    founded = False
                    for ele in xnus:
                        neutron = int(ele[0])
                        proton = int(ele[1])
                        xnu = float(ele[3])
                        inspect_ele = inspect[i]
                        if(neutron+proton == int(inspect_ele[1]) and proton == int(inspect_ele[2])):
                            XNUS[i].append(xnu)
                            founded = True
                            break
                    if(not founded):
                        XNUS[i].append(0)
    for i in range(len(inspect)):
        MASSES[i] += numpy.interp(time_points,Times,XNUS[i]) * mass



    for i in range(nrn_result.shape[0]):
        total_weighed_result[i] += mass*nrn_result[i,2]/M_S

plt.plot([i for i in range(1,100)],total_weighed_result[1:100])
plt.plot([i for i in range(1,100)],INIT_MASS[1:100])
plt.yscale("log")
plt.xlabel("Mass Number A")
plt.ylabel(r"$Mass \qquad Fraction/M_{\odot}$")
plt.title(f"Nuclei Mass Fraction(Ejecta Mass:{total_mass/M_S:.2f}"+r"$M_{\odot}$)")
plt.show()
for i in range(len(inspect)):
    plt.figure()
    plt.plot(time_points,MASSES[i],label=inspect[i][0])
    plt.yscale("log")
    plt.legend()
    plt.show()



