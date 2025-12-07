import numpy as np
import tracer_reader
modelname = "14m_new"
a = tracer_reader.TracerReader(modelname=modelname,subdirectory="models")
f = open(f"trajectory/{modelname}_ejected_tracers_list.txt","r")
ids = []
for lines in f:
    if("_" in lines):
        id = int(lines.split("_")[1])
        print(id)
        ids.append(id)
a.GetTrajectories(ids,f"trajectory/CCSN_Traj_{modelname}")