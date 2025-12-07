import tracer_reader
import matplotlib.pyplot as plt


t = tracer_reader.TracerReader(subdirectory="models",modelname="21m_old")
ebind = t.get_binding_energy()[:,:,0]
t.GetWinNetSeedFiles("trajectories")