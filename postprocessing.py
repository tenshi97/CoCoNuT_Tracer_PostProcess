import os

import h5py
import h5reader
import tracer_reader
import numpy as np
import matplotlib.pyplot as plt
a = tracer_reader.TracerReader(modelname="14m_old",subdirectory="models")
#a.ShowTrajectory(1000,False)
#8600 for new and 16194 for old
plt.rcParams['axes.titlesize'] = 20  # 图标题 title
plt.rcParams['axes.labelsize'] = 20 # x/y轴标签 xlabel/ylabel
plt.rcParams['xtick.labelsize'] = 10  # x轴刻度
plt.rcParams['ytick.labelsize'] = 10  # y轴刻度
plt.rcParams['legend.fontsize'] = 14  # 图例 legend
a.ShowInitial()
a.ShowFinal()
a.set_index(a.n)
a.contour_entropy()
a.contour_vex()
a.contour_any("pre")
# a.get_allgridmass()