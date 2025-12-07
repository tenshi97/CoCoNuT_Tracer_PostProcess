import os

import numpy as np

import phycon_and_nuc_table as nt
import h5reader
import shock_tracker
import diag_2d
import matplotlib.pyplot as plt

model_name = '14m_new'
h=h5reader.hydrof(model=model_name, subdirectory="models")
xif = h.xzr()
yif = h.yzr()
r_len = len(h.xzn())
theta_len = len(h.yzn())
pc_nuc = nt.pc_nuc
wc_mb = nt.wc_mb
wc_me = nt.wc_me
pc_meverg = nt.pc_meverg
pc_mb = nt.pc_mb
pc_msol = nt.pc_msol
pc_cl = nt.pc_cl

ebind = diag_2d.get_binding_energy(h)

def grid_finder(x,y):
    st = 0
    ed = r_len - 1
    x_cell = -1
    y_cell = -1
    while(st<ed):
        mid = (st+ed)//2
        if(x<xif[mid]):
            ed = mid
        else:
            st = mid+1
    x_cell = st
    st = 0
    ed = theta_len - 1
    while(st<ed):
        mid = (st+ed)//2
        if(y<yif[mid]):
            ed = mid
        else:
            st = mid+1
    y_cell = st
    return (x_cell,y_cell)

tracer_ids = np.array([int(round(x)) for x in h.tracer_id()])
tracer_num = len(tracer_ids)
tracer_data = [[] for x in range(tracer_num)]
record_ebind = open(f"F:/compastro/Tracers16000/trajectory/{model_name}_tracers_ebind_list.txt","w")
record_ejected= open(f"F:/compastro/Tracers16000/trajectory/{model_name}_ejected_tracers_list.txt","w")
record_not_ejected= open(f"F:/compastro/Tracers16000/trajectory/{model_name}_not_ejected_tracers_list.txt","w")

tracer_ids = np.array([int(round(x)) for x in h.tracer_id()])
last_ts_data = [() for i in range(tracer_num)]
unbound_tracers = []
print(h.index)
count = 0
for i in range(tracer_num):
    gid = tracer_ids[i]
    r = h.tracer_x()[i]
    theta = h.tracer_y()[i]
    (x, y) = grid_finder(r, theta)
    last_ts_data[gid-1]  = (gid,r,theta,x,y)
for i in range(tracer_num):
    gid = last_ts_data[i][0]
    r = last_ts_data[i][1]
    theta = last_ts_data[i][2]
    x = last_ts_data[i][3]
    y = last_ts_data[i][4]
    record_ebind.write(f"tracer_{gid} coordinate r:{r}({x}) theta:{theta}(y) binding energy{ebind[x][y][0]}\n")
    if (ebind[x][y][0] > 0):
        unbound_tracers.append(gid)
        record_ejected.write(f"tracer_{gid}\n")
        count += 1
    else:
        record_not_ejected.write(f"tracer_{gid}\n")
print(count)
record_ebind.close()
record_ejected.close()
record_not_ejected.close()