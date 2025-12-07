import numpy as np
import tracer_reader
import matplotlib.pyplot as plt
trs = tracer_reader.TracerReader(subdirectory="models",modelname="14m_new")
trs.set_index(trs.n-1)
