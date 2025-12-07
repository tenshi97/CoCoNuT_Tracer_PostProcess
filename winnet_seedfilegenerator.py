import os

import h5py
import tracer_reader
import phycon_and_nuc_table
from matplotlib import pyplot as plt

from tracer_reader import TracerReader

modelname = "14m_new"
a = tracer_reader.TracerReader(model=modelname, subdirectory="models")
print(r_len,theta_len)

a.set_index(0)
ejected_files = open(f"./trajectory/{modelname}_ejected_tracers_list.txt")
ejected_tracers = []
for lines in ejected_files:
    ejected_tracers.append(int(lines.split("_")[1]))
print(len(ejected_tracers))
for i in range(len(a.tracer_id())):
    tid = int(a.tracer_id()[i])
    if(tid not in ejected_tracers):
        print(tid)
        continue
    print(tid)
    output = open(f"./trajectory/CCSN_Seed_{modelname}/seed_{modelname}_{tid}","w")
    tx = a.getX()
    ty = a.getY()
    (ix,iy) = a.grid_finder(tx,ty)
    Xnu = a.xnu()[ix][iy][0][:]
    print(tid,tx,ty,ix,iy)
    output.write("#    A    Z       X\n")
    for i in range(len(phycon_and_nuc_table.pc_nuc)-3):
        Ax = int(phycon_and_nuc_table.pc_nuc[i][2])
        Zx = int(phycon_and_nuc_table.pc_nuc[i][1])
        Xx = Xnu[i]
        Xstr = f"{Xx:.3e}"
        output.write(f"{str(Ax).rjust(5)}{str(Zx).rjust(5)}{Xstr.rjust(12)}"+"\n")