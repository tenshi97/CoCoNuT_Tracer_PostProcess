import os

import astropy.constants.astropyconst40
import h5py
from pathlib import Path
import numpy as np
folder_path = Path("F:/compastro/Tracers16000/models")
folder_path_2 = Path("F:/compastro/Tracers16000/trajectory/CCSN_CoCoNuT_Trajectory/")
file_prefix = "14m_old"
columns = [0,4,3,2,5,6,7,8,9,10,11]
firstheader = f"# t[s]      T [GK]       rho[g/cm^3]   r [km]      Ye           Le [erg/s]  Lebar [erg/s]  Ltau [erg/s] Ee [MeV]  Eebar [MeV]   Etau[MeV]\n"
secondheader = "#----------------------------------------------------------------------------------------------------------------------------------------\n"
for files in folder_path.glob(f"{file_prefix}*.dat"):
    print(files)
    a = np.loadtxt(files,skiprows=0)
    b = a[:,columns]
    b[:,1] /= 1e9
    b[:,3] /= 1e5
    filename_out = str(folder_path_2) +"/winnet_"+files.name.replace(".dat","")
    out = open(filename_out,"w")
    formatter = lambda x: f"{x:.4e}".replace("e", "e")
    matrix_output = np.array2string(b, separator=" ",formatter={'all': formatter})
    print(matrix_output)
    np.savetxt(out,b,fmt="%.4e", delimiter=" ")
    out.close()
    with open(filename_out, "r", encoding="utf-8") as f:
        original_content = f.readlines()
    updated_content = [firstheader,secondheader] + original_content
    with open(filename_out, "w", encoding="utf-8") as f:
        f.writelines(updated_content)




