import numpy as np

import h5reader
import matplotlib.pyplot as plt
import math
import re
import copy
from collections import deque
from openpyxl import Workbook
from openpyxl.styles import PatternFill
from openpyxl.styles import Alignment
import random
h=h5reader.hydrof(model='s15_16000')
h.set_index(h.ngroups-1)
contour = plt.contourf(h.yzn(),h.xzn(),np.log(h.tem()[:,:,0]))
cbar = plt.colorbar(contour)
plt.yscale("log")
plt.show()
tracer_ids = np.array([int(round(x)) for x in h.tracer_id()])
tracer_num = len(tracer_ids)
print(tracer_num)
tracer_data = [[] for i in range(tracer_num)]
for t in range(100,h.ngroups):
    h.set_index(t)
    tracer_ids = np.array([int(round(x)) for x in h.tracer_id()])
    for i in range(len(tracer_ids)):
        gid = tracer_ids[i]
        tracer_data[gid-1].append((h.time(),h.tracer_x()[i],h.tracer_y()[i]))
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
sc = ax.scatter([],[],s=2)
ax.set_rlim(1e5,1e7)
ax.set_yscale("log")
plt.show()
colors = []
steps = []
for t in range(100,h.ngroups):
    data_t = [x[t] for x in tracer_data]
    THETA_t = []
    R_t = []
    for i in range(0,tracer_num):
        if(i%40<=207 and i%40>=196):
            THETA_t.append(data_t[i][2])
            R_t.append(data_t[i][1])
    if(t==0):
        ax.scatter(THETA_t, R_t, s=1)
    else:

        if(t%200==0):
            random_color = (random.random(), random.random(), random.random())
            colors.append(random_color)
            steps.append(t)  # 保存步数
            ax.scatter(THETA_t,R_t,color=random_color,s=2,label=f'Step {t}')

        sc.set_offsets(np.c_[THETA_t,R_t])
    fig.canvas.draw()
    plt.pause(0.03)

    plt.title("Distribution at time point "+str(t))
plt.ioff()
ax.legend(title="Recorded Steps", bbox_to_anchor=(1.1, 1.05))
plt.show()

# data = []
# procs = 48
# tracer_per_proc = 40
# total_tracer = procs*tracer_per_proc
# tracer_data = [[] for x in range(0,total_tracer)]
# time_data = [[] for x in range(0,3700)]
# for i in range(0,48):
#     print(i)
#     timestep_p = 0
#     with open("F:/compastro/hdf5output/debug/filen"+str(i),mode='r') as file:
#         data_proc = []
#         tracer_num_record = []
#         current_timestep = 0
#         data_proc_timestep = {"timestep":0,"tracers": []}
#         local_id = 0
#         global_id = 0
#         pos_r = 0
#         cellindex_r = 0
#         pos_theta = 0
#         cellindex_theta = 0
#         for line in file:
#             if(line.find("Cell")!=-1):
#                 continue
#             if(line.find("IPROC")!=-1):
#                 continue
#             if(line.find("v_r")!=-1 or line.find("v_theta")!=-1):
#                 continue
#             if(line.find("INIT")!=-1):
#                 continue
#             pattern1 = re.compile(r"TRACER NO =\s+(\d+)\s+TIME =\s+(\d+.\d+E?[-+]?\d+)\s+(\d+)")
#             match1 = pattern1.search(line)
#             if(match1):
#                 local_id = int(match1.group(1))
#                 time = float(match1.group(2))
#                 timestep = int(match1.group(3))
#                 if(current_timestep!=timestep):
#                     if(current_timestep!=0):
#                         timestep_p+=1
#                         tracer_num_record.append(len(data_proc_timestep["tracers"]))
#                         data_proc.append(copy.deepcopy(data_proc_timestep))
#                     data_proc_timestep["timestep"] = timestep
#                     data_proc_timestep["time"] = time
#                     data_proc_timestep["tracers"] = []
#                     current_timestep = timestep
#                 continue
#             pattern2 = re.compile(r"GLOBAL ID=\s+(\d+)")
#             match2 = pattern2.search(line)
#             if match2:
#                 global_id = int(match2.group(1))
#                 continue
#             pattern3 = re.compile(r"position r = \s+(\d+.\d+E?-?\d+)\s+(\d+)")
#             match3 = pattern3.search(line)
#             if match3:
#                 pos_r = float(match3.group(1))
#                 cellindex_r = int(match3.group(2))
#                 continue
#             pattern4 = re.compile(r"position theta = \s+(\d+.\d+E?-?\d+)\s+(\d+)")
#             match4 =  pattern4.search(line)
#             if match4:
#                 pos_theta = float(match4.group(1))
#                 cellindex_theta = int(match4.group(2))
#                 data_proc_timestep["tracers"].append((local_id,global_id,pos_r,cellindex_r,pos_theta,cellindex_theta))
#                 tracer_data[global_id-1].append((timestep,time,pos_r,pos_theta))
#                 time_data[timestep_p].append((global_id,pos_r,pos_theta,time))
#                 continue
#         tracer_num_record.append(len(data_proc_timestep["tracers"]))
#         data_proc.append(copy.deepcopy(data_proc_timestep))
# #    for ele in data_proc:
# #        wt.write(str(ele)+"\n")
#     data.append(data_proc)
# test_tracer = 900
#
# tracer_data_i = sorted(tracer_data[test_tracer],key=lambda x:x[0])
# T = [x[1] for x in tracer_data_i]
# R = [x[2] for x in tracer_data_i]
# THETA = [x[3] for x in tracer_data_i]
# plt.plot(T,THETA)
# plt.plot(time_steps,tracer_tpos)
# plt.show()
# plt.plot(T,R)
# plt.plot(time_steps,tracer_rpos)
# plt.show()



