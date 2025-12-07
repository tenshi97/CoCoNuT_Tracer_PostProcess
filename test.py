import tracer_reader
import matplotlib.pyplot as plt


t = tracer_reader.TracerReader(subdirectory="models",modelname="21m_old")
ebind = t.get_binding_energy()[:,:,0]
t.set_index(t.n)
tracerX = t.get_element("trx")
tracerY = t.get_element("try")
tracerID = t.get_element("trid")
ejected_tracers = []
print(t.getX())
# for i in range(t.tracerNum):
#     x = tracerX[i]
#     y = tracerY[i]
#     id = tracerID[i]
#     ix,iy = t.grid_finder(x,y)
#     if(ebind[ix][iy]>0):
#         ejected_tracers.append(id)
# print(len(ejected_tracers))
# t.GetTrajectories_new(masks=ejected_tracers)