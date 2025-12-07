import numpy as np
import matplotlib.pyplot as plt
import h5reader
import astropy
import diag_2d
a = h5reader.hydrof(model="20m_old",subdirectory="models")
ebind = diag_2d.get_binding_energy(a)
f = open("20m_old_mass_assign.txt")
ids = []
plt.contour(a.xzn(),a.yzn(),ebind[:,:,0].transpose(),levels=np.linspace(0,np.max(ebind[:,:,0]),3))
for line in f:
    if("Total" not in line):
        id = int(line.split(":")[0])
        ene = float(line.split(":")[1].strip())
        if(ene<1):
            continue
        ids.append(id)
for i in range(len(a.tracer_id())):
    print(i,a.tracer_id()[i])
    if(int(a.tracer_id()[i]) in ids):
        plt.scatter(a.tracer_x()[i],a.tracer_y()[i],c="red")
plt.xscale("log")
plt.show()
