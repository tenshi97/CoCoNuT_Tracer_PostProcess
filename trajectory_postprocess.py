import os

import astropy.constants.astropyconst40
import h5py
from pathlib import Path
import numpy as np
a = h5py.File("F:/compastro/Tracers16000/12m_new_cmb.h5")
print(list(a['vertex']))
print(a['vertex']['n_zones'][()])
print(list(a['vertex']['rri'])[0:3])
print(list(a['vertex']['rhoi']))
os.system("pause")
cellnum = a['vertex']['n_zones'][()]
xzls = a["vertex"]["rri"][:-1]
xzrs = a["vertex"]["rri"][1:]
m_enclosed = []
lr = 0
m = 0
rri = list(a['vertex']['rri'])
rhoi = list(a['vertex']['rhoi'])
print(rhoi)
for i in range(cellnum):
    r = rri[i]
    rho = rhoi[i]
    dr = r - lr
    if(i==0):
        dm = 4*np.pi/3 * (r**3) * rho
    else:
        dm = 4*np.pi*(r**2) * dr * rho
    m += dm
    m_enclosed.append(m)
    lr = r
print(m_enclosed)
os.system("pause")
print(xzls[0])
print(xzrs[0])
folder_path = Path("F:/compastro/Tracers16000/trajectory/ejected")
for file_path in folder_path.rglob('*.dat'):
    if("traj" in file_path.name):
        continue
    f = open(file_path,"r")
    line1 = f.readline()
    ri = float(line1.split(" ")[2])
    em = 0
    for i in range(cellnum):
        if(ri>xzls[i] and ri<=xzrs[i]):
            em = m_enclosed[i]
            print(f"this tracer belongs to Cell {i}, enclosed mass is {em}")
    original_text = f.read()
    headers = (f"# version:  10001\n"
               f"# model: s15-2007\n"
               f"# mass:   {em}\n"
               f"# \n"
               "#                   time                       dt                   radius                  density              temperature                      Y_e                  L_nu_e-                  L_nu_e+           L_nu_{tau+tau}                  E_nu_e-                  E_nu_e+            E_nu_{mu/tau}\n"
               "#                      s                        s                       cm                    g/cm3                        K                                             erg/s                    erg/s                    erg/s                      MeV                      MeV                      MeV\n"
               "# \n")
    new_text = headers+original_text
    filename = file_path.name
    filename = filename.replace("s15_16000","12m_traj")
    w = open(str(folder_path)+"/"+"new_"+filename,"w")
    w.write(new_text)
    print(headers)




