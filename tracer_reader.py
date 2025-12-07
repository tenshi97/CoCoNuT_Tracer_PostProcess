import h5py
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from phycon_and_nuc_table import pc_cl, wc_mn,pc_meverg,pc_nuc,pc_gc,pc_mb,pc_ergmev,wc_me
from astropy.constants import c,G
from matplotlib.colors import LogNorm,Normalize
from matplotlib.animation import FuncAnimation
import random
import os
import sys
from tqdm import tqdm

import phycon_and_nuc_table



class TracerReader:
    def __init__(self, modelname,subdirectory = ""):
        self.subdir = ""
        self.modelname = modelname
        self.timepoints = []
        self.steps = []
        self.files = []
        self.index_to_file = []
        # if subdirectory is None or subdirectory == "" :
        #     subdir = ""
        # else:
        #     subdir = f"./{subdirectory}/"
        # print(subdir)
        # subdir = Path(subdir)
        if not subdirectory:  # None 或 "" -> 当前目录
            base = Path(".")
        else:
            base = Path(subdirectory)  # 相对路径也 OK

        print("base dir:", base)
        files = [
            f for f in base.iterdir()
            if f.is_file() and f.name.startswith(modelname) and f.suffix.startswith(".o")
        ]
        files = sorted(files, key=lambda f: f.name)
        self.files = [str(f.as_posix()) for f in files]
        print(self.files)
        # self.files = [subdir.name + "/" + f.name for f in subdir.iterdir() if (f.is_file() and f.name.startswith(modelname))]
        self.n = -1
        file_count = 0
        if(len(self.files)==0):
            raise ValueError("No Results have been read")
        for h5file in self.files:
            print(h5file)
            h5obj = h5py.File(h5file)
            for step in h5obj.keys():
                self.steps.append(step)
                self.timepoints.append(h5obj[step]["time"][()])
                self.index_to_file.append(file_count)
                self.n +=1
            file_count +=1
        if(self.n ==-1):
            raise ValueError("No Results have been read")
        self.index = 0
        self.id = 1
        ff = h5py.File(self.files[0])
        self.tracerNum = ff[self.steps[0]]["trid"].shape[0]
        self.gridnx = ff[self.steps[0]]["xzn"].shape[0]
        self.gridny = ff[self.steps[0]]["yzn"].shape[0]
        self.gridxzl = np.array(ff[self.steps[0]]["xzl"]).transpose()
        self.gridyzl = np.array(ff[self.steps[0]]["yzl"]).transpose()
        self.gridxzr = np.array(ff[self.steps[0]]["xzr"]).transpose()
        self.gridyzr = np.array(ff[self.steps[0]]["yzr"]).transpose()
        self.gridxzn = np.array(ff[self.steps[0]]["xzn"]).transpose()
        self.gridyzn = np.array(ff[self.steps[0]]["yzn"]).transpose()
        print(f"{self.tracerNum} tracers have been read")
        print(f"Total Steps{len(self.steps)},{self.n}")
    def set_time(self,time):
        if(time>=self.timepoints[self.n]):
            print(f"Time exceeds the max time, set to index:{self.n}(time:{self.timepoints[self.n]})")
            self.index = self.n
            return
        if(time<0):
            print(f"Time exceeds the min time, set to index:{0}(time:{self.timepoints[0]})")
            self.index = 0
            return
        st, ed = 0,self.n
        while st<ed:
            mid = (st+ed)//2
            if (time < self.timepoints[mid]):
                ed = mid
            else:
                st = mid + 1
        if st == 0:
            self.index = 0
        else:
            if(np.abs(time-self.timepoints[st-1]) < np.abs(time-self.timepoints[st])):
                self.index = st-1
            else:
                self.index = st
        print(f"Current index is {self.index}")
    def set_index(self,index):
        if(index>self.n):
            self.index = self.n
            return
        if(index<0):
            self.index = 0
        self.index = index

    def set_id(self,tid):
        if(tid<0 or tid>self.tracerNum):
            raise ValueError("Tracer Id out of boundary")
        self.id = tid

    def get_element(self,name):
        if name is not None:
            file = self.files[self.index_to_file[self.index]]
            #print(f"open {file}")
            h5obj = h5py.File(file)
            if name in h5obj[self.steps[self.index]].keys():
                #print(f"open {self.steps[self.index]}")
                ret = np.array(h5obj[self.steps[self.index]][name])
                return ret.transpose()
            else:
                if(name=="pre_grad"):
                    pre = self.get_element("pre")
                    ret = np.gradient(pre, axis=0) / pre
                    return ret
                print(f"{name} is not a valid keyword")
                raise ValueError(f"{name} is not a valid keyword")
        else:
            raise ValueError("No Keywords")

    def getX(self):
        tracerID = self.get_element("trid")
        tracerX = self.get_element("trx")
        tracerID = list(tracerID)
        ret = []
        for i in range(len(tracerID)):
            tracerID[i] = int(tracerID[i])
            if(tracerID[i]==self.id):
                return float(tracerX[i])
        raise ValueError("Not Found")
    def getY(self):
        tracerID = self.get_element("trid")
        tracerY = self.get_element("try")
        tracerID = list(tracerID)
        ret = []
        for i in range(len(tracerID)):
            tracerID[i] = int(tracerID[i])
            if(tracerID[i]==self.id):
                return float(tracerY[i])
        raise ValueError("Not Found")

    def XTimeSeriesForOneTracer(self,start=1,end=-1,scatter=False):
        tracerX_plot = []
        time_plot = []
        if(start<0 or start>=self.n):
            start = 0
        if(end<0 or end>=self.n):
            end = self.n-1
        for i in range(start,end+1):
            self.set_index(i)
            if(scatter):
                plt.scatter(self.timepoints[i],self.getX(),c="blue")
                continue
            tracerX_plot.append(self.getX())
            time_plot.append(self.timepoints[i])
        if(not scatter):
            plt.plot(time_plot,tracerX_plot)
        plt.title(f"Radial Coordinates of Tracer {self.id} Time Series")
        plt.show()

    def YTimeSeriesForOneTracer(self,start=1,end=-1,scatter=False):
        tracerY_plot = []
        time_plot = []
        plt.figure()
        if(start<0 or start>=self.n):
            start = 0
        if(end<0 or end>=self.n):
            end = self.n-1
        for i in range(start,end+1):
            self.set_index(i)
            if(scatter):
                plt.scatter(self.timepoints[i],self.getY(),c="blue")
                continue
            tracerY_plot.append(self.getY())
            time_plot.append(self.timepoints[i])
        if(not scatter):
            plt.plot(time_plot,tracerY_plot)
        plt.title(f"Theta Coordinates of Tracer {self.id} Time Series")
        plt.show()

    def contour_vex(self):
        self.set_index(self.n)
        vex = self.get_element("vex")
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        vmax = np.max(np.abs(vex))
        vmin = -vmax
        contour = plt.contourf(xzn,yzn,vex[:,:,0].transpose(),levels=20, cmap='seismic',vmax=vmax,vmin=vmin)

        cbar = plt.colorbar(contour)
        ebind = self.get_binding_energy()
        cs = plt.contour(xzn, yzn, ebind[:,:,0].transpose(), levels=[-1e16,0,1e18,1e19], colors=['blue', 'red','green','yellow'], linewidths=2)
        plt.clabel(cs, fmt={-1e16:"E=-1e16",0: 'E=0', 1e18: 'E=1e18',1e19:'E=1e19'}, colors=['blue', 'red','green','yellow'])
        plt.xscale("log")
        plt.xlabel("Radius/cm")
        plt.ylabel("Theta")
        plt.tight_layout()
        plt.savefig(f"vex_contour_{self.modelname}_vex.png",dpi=640,bbox_inches='tight')
        plt.show()
        plt.close()

    def contour_any(self,name):
        if(name!="ebind"):
            par = self.get_element(name)
        else:
            par = self.get_binding_energy()
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        contour = plt.contourf(xzn, yzn, par[:, :, 0].transpose(), levels=20, cmap='viridis',norm=LogNorm())
        cbar = plt.colorbar(contour)
        ebind = self.get_binding_energy()
        cs = plt.contour(xzn, yzn, ebind[:, :, 0].transpose(), levels=[-1e16, 0, 1e18, 1e19],
                         colors=['blue', 'red', 'green', 'yellow'], linewidths=2)
        plt.clabel(cs, fmt={-1e16: "E=-1e16", 0: 'E=0', 1e18: 'E=1e18', 1e19: 'E=1e19'},
                   colors=['blue', 'red', 'green', 'yellow'])
        plt.xscale("log")
        plt.xlabel("Radius/cm")
        plt.ylabel("Theta")
        plt.tight_layout()
        plt.savefig(f"{name}_contour_{self.modelname}.png",dpi=640,bbox_inches='tight')
        plt.close()

    def get_gridmass(self,x,y):
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        xzl = self.get_element("xzl")
        yzl = self.get_element("yzl")
        xzr = self.get_element("xzr")
        yzr = self.get_element("yzr")
        vex = self.get_element("vex")
        vey = self.get_element("vey")
        den = self.get_element("den")
        phi = self.get_element("phi")
        DV = 2 * np.pi * (np.cos(yzl[y]) - np.cos(yzr[y])) * (xzr[x] ** 3 - xzl[x] ** 3) / 3
        LF = 1.0 / np.sqrt(1 - (vex[x][y][0] ** 2 + vey[x][y][0] ** 2) / (c.cgs.value**2))
        DM = den[x,y,0] * DV * (phi[x,y,0]**6) * LF
        return DM

    def get_allgridmass(self):
        self.set_index(self.n)
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        xzl = self.get_element("xzl")
        yzl = self.get_element("yzl")
        xzr = self.get_element("xzr")
        yzr = self.get_element("yzr")
        vex = self.get_element("vex")
        vey = self.get_element("vey")
        den = self.get_element("den")
        phi = self.get_element("phi")
        xnu = self.get_element("xnu")[:,:,:,:20]
        mnu_cor = np.zeros(20)
        print(xnu.shape)
        xs = xzn.shape[0]
        ys = yzn.shape[0]
        M_all = 0
        M_ejected = 0
        M_ejected_f = 0
        ebind = self.get_binding_energy()
        shock_loc = 1e9
        for x in range(xs):
            for y in range(ys):
                DV = 2 * np.pi * (np.cos(yzl[y]) - np.cos(yzr[y])) * (xzr[x] ** 3 - xzl[x] ** 3) / 3
                LF = 1.0 / np.sqrt(1 - (vex[x][y][0] ** 2 + vey[x][y][0] ** 2) / (c.cgs.value ** 2))
                DM = den[x, y, 0] * DV * (phi[x, y, 0] ** 6) * LF
                M_all += DM
                if(ebind[x,y,0]>0):
                    M_ejected += DM
                if(ebind[x,y,0]<0 and xzn[x]>shock_loc):
                    M_ejected_f += DM
                    mnu_cor += DM* xnu[x,y,0,:]


        print(M_all/1.989e33)
        print(M_ejected/1.989e33)
        print(M_ejected_f/1.989e33)
        np.savetxt(f"{self.modelname}_mnu_cor.txt",mnu_cor, fmt='%.8e', delimiter=' ')
        f = open(f"{self.modelname}_ejectamass.txt","w")
        f.write(f"Innermost Ejecta:{M_ejected/1.989e33:.8f} Solar Mass\n")
        f.write(f"Correct Region Mass:{M_ejected_f/1.989e33:.8f} Solar Mass\n")
        f.write(f"Total Grid Mass:{M_all/1.989e33:.8f} Solar Mass\n")
        f.close()
    def count_ejecta(self):
        self.set_index(self.n)
        shock_loc = 1e9
        ebind = self.get_binding_energy()[:,:,0]
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        #
        # for i in range(ebind.shape[0]):
        #     for j in range(ebind.shape[1]):
        #         if(ebind[i][j]>0 or (ebind[i][j]<0 and xzn[i,j,0]>shock_loc)):


    def contour_entropy(self):
        self.set_index(self.n)
        sto = self.get_element("sto")
        xzn = self.get_element("xzn")
        yzn = self.get_element("yzn")
        contour = plt.contourf(xzn,yzn,sto[:,:,0].transpose(),levels=20, cmap='viridis')
        cbar = plt.colorbar(contour)
        cbar.set_label('Entropy')
        ebind = self.get_binding_energy()
        cs = plt.contour(xzn, yzn, ebind[:,:,0].transpose(), levels=[-1e16,0,1e18,1e19], colors=['blue', 'red','green','yellow'], linewidths=2)
        plt.clabel(cs, fmt={-1e16:"E=-1e16",0: 'E=0', 1e18: 'E=1e18',1e19:'E=1e19'}, colors=['blue', 'red','green','yellow'])
        plt.xscale("log")
        plt.xlabel("Radius/cm")
        plt.ylabel("Theta")
        plt.tight_layout()
        plt.savefig(f"entropy_contour_{self.modelname}.png",dpi=640,bbox_inches='tight')
        plt.close()

    def get_binding_energy(self):

        self.set_index(self.n)

        vex = self.get_element("vex")
        vey = self.get_element("vey")
        vez = self.get_element("vez")
        ene = self.get_element("ene")
        xnu = self.get_element("xnu")
        pre = self.get_element("pre")
        den = self.get_element("den")
        alpha = self.get_element("alpha")
        w = 1 / np.sqrt(1 - (vex ** 2 + vey ** 2 + vez ** 2) / pc_cl ** 2)
        mnuc = np.zeros(21)
        mnuc[0:20] = (pc_nuc[:, 3] + pc_nuc[:, 2] * (pc_mb * pc_cl ** 2) * pc_ergmev) / pc_nuc[:, 2]
        mnuc[20] = wc_me
        e_int = ene * pc_gc / pc_cl ** 2
        e_nuc = np.sum([xnu[:, :, :, s] * mnuc[s] * pc_meverg / pc_mb for s in range(21)], axis=0)
        e_cor = (wc_mn - 8.8) * pc_meverg / pc_mb
        e_int = (e_int + pc_cl ** 2 * (1 - w)) / w + pre * (1 - w ** 2) / (den * w ** 2)
        e_int -= e_nuc * w
        e_int += e_cor * w

        e_bind = alpha * (w * (
                    c.cgs.value ** 2 + e_int + pre / den) - pre / den / w) - c.cgs.value ** 2
        # xzn = self.get_element("xzn")
        # yzn = self.get_element("yzn")
        # contour = plt.contourf(xzn,yzn,e_bind[:,:,0].transpose(),levels=20, cmap='viridis')
        # cbar = plt.colorbar(contour)
        # cbar.set_label('Binding Energy')
        # plt.xscale("log")
        # plt.tight_layout()
        # plt.xlabel("Radius/cm")
        # plt.ylabel("Theta")
        # plt.show()
        # plt.close()
        return e_bind

    def grid_finder(self,x, y):
        r_len = self.gridnx
        theta_len = self.gridny
        xif = list(self.gridxzr)
        yif = list(self.gridyzr)
        st = 0
        ed = r_len - 1
        x_cell = -1
        y_cell = -1
        while (st < ed):
            mid = (st + ed) // 2
            if (x < xif[mid]):
                ed = mid
            else:
                st = mid + 1
        x_cell = st
        st = 0
        ed = theta_len - 1
        while (st < ed):
            mid = (st + ed) // 2
            if (y < yif[mid]):
                ed = mid
            else:
                st = mid + 1
        y_cell = st
        return (x_cell, y_cell)

    def getCellIndexByID(self,id):
        self.id = id
        x = self.getX()
        y = self.getY()
        return self.grid_finder(x,y)

    #def ContourTracers(self,masks):
    def ShowInitial(self):
        self.index = self.n
        eject_mask = []
        tracerx = self.get_element("trx")
        tracery = self.get_element("try")
        tracerids = self.get_element("trid")
        tnum = len(tracerids)
        eject_count = 0
        ebind = self.get_binding_energy()[:,:,0]
        for i in range(tnum):
            x,y = self.grid_finder(tracerx[i],tracery[i])
            id = tracerids[i]
            if(ebind[x,y]>0):
                eject_count += 1
                eject_mask.append(id)
        print(f"Total ejected tracer number is {eject_count}")
        R, Theta = np.meshgrid(self.gridxzn, self.gridyzn)
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
        ax.set_rlim(1e6,1.5e9)
        ax.set_yscale("log")
        ax.set_thetalim(0, np.pi)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi)
        self.set_index(0)
        tracerx = self.get_element("trx")
        tracery = self.get_element("try")
        tracerids = self.get_element("trid")
        ejected = []
        not_ejected = []
        for i in range(len(tracerx)):
            if(tracerids[i] in eject_mask):
                ejected.append([tracerx[i],tracery[i]])
            else:
                not_ejected.append([tracerx[i],tracery[i]])
        ejected = np.array(ejected)
        not_ejected = np.array(not_ejected)
        ax.scatter(ejected[:,1], ejected[:,0], s=1, marker="o",c="#39c5bb",label="Ejected Tracers")
        ax.scatter(not_ejected[:,1],not_ejected[:,0], s=1, marker="o",c="purple",label="Unejected Tracers")
        plt.title(f"Initial Distribution of Tracer Particles for {self.modelname} Model")
        plt.ylabel("Theta")
        plt.xlabel("Radius/cm")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"figures/{self.modelname}_initdist.png",dpi=1080,bbox_inches='tight')
        plt.close(fig)

    def ShowFinal(self):
        self.index = self.n
        eject_mask = []
        tracerx = self.get_element("trx")
        tracery = self.get_element("try")
        tracerids = self.get_element("trid")
        tnum = len(tracerids)
        eject_count = 0
        ebind = self.get_binding_energy()[:,:,0]
        for i in range(tnum):
            x,y = self.grid_finder(tracerx[i],tracery[i])
            id = tracerids[i]
            if(ebind[x,y]>0):
                eject_count += 1
                eject_mask.append(id)
        print(f"Total ejected tracer number is {eject_count}")
        R, Theta = np.meshgrid(self.gridxzn, self.gridyzn)
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
        emax = np.max(ebind)
        contourf = ax.contourf(Theta, R, ebind.transpose(), 50, cmap='seismic',vmax = emax,vmin = -emax)
        cbar = plt.colorbar(contourf, ax=ax, pad=0.1)
        ax.set_rlim(1e5,1e10)
        ax.set_yscale("log")
        ax.set_thetalim(0, np.pi)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi)
        tracerx = self.get_element("trx")
        tracery = self.get_element("try")
        tracerids = self.get_element("trid")
        ejected = []
        not_ejected = []
        plt.tight_layout()
        for i in range(len(tracerx)):
            if(tracerids[i] in eject_mask):
                ejected.append([tracerx[i],tracery[i]])
            else:
                not_ejected.append([tracerx[i],tracery[i]])
        ejected = np.array(ejected)
        not_ejected = np.array(not_ejected)
        ax.scatter(ejected[:,1], ejected[:,0], s=6, marker="o",c="#39c5bb",label="Ejected Tracers")
        ax.scatter(not_ejected[:,1],not_ejected[:,0], s=6, marker="x",c="purple",label="Unejected Tracers")
        plt.title(f"Final Distribution of Tracer Particles for {self.modelname} Model")
        plt.legend()
        plt.grid(True)
        plt.ylabel("Theta")
        plt.xlabel("Radius/cm")
        plt.savefig(f"figures/{self.modelname}_finaldist.png",dpi=1080)
        plt.close(fig)

    def ShowAnimation(self,timestep,snapstep,bg=None,rmin=1e6,rmax=1e10,st=-1,ed=-1):

        self.index = self.n
        eject_mask = []
        snapc =  0
        tracerx = self.get_element("trx")
        tracery = self.get_element("try")
        tracerids = self.get_element("trid")
        tnum = len(tracerids)
        eject_count = 0
        ebind = self.get_binding_energy()[:, :, 0]
        if(st<0):
            st = 0
        if(ed<0):
            ed = self.n

        for i in range(tnum):
            x, y = self.grid_finder(tracerx[i], tracery[i])
            id = tracerids[i]
            if (ebind[x, y] > 0):
                eject_count += 1
                eject_mask.append(id)
        print(f"Total ejected tracer number is {eject_count}")
        R, Theta = np.meshgrid(self.gridxzn, self.gridyzn)
        plt.ion()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
        ax.set_rlim(rmin,rmax)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_thetamin(0)
        ax.set_thetamax(180)
        sc = ax.scatter([], [], s=2)
        x1,y1 = self.grid_finder(rmin,np.pi/2)
        x2,y2 = self.grid_finder(rmax,np.pi/2)
        if(bg!=None):
            field = self.get_element(bg)
            print(field.shape)
            field = self.get_element(bg)[:, :, 0]
            print(field.shape)
            slice_ = field[x1:x2, :]
            vmax = np.max(slice_)
            vmin = np.min(slice_)
            for i in range(tnum+1):
                field = self.get_element(bg)[:, :, 0]
                slice_ = field[x1:x2, :]
                if(np.max(slice_)>vmax):
                    vmax = np.max(slice_)
                if(np.min(slice_)<vmin):
                    vmin = np.min(slice_)
            if(vmin<0):
                norm = Normalize(vmin=vmin, vmax=vmax)
            else:
                norm = LogNorm(vmin=vmin, vmax=vmax)

            cf = ax.contourf(Theta, R, field.transpose(), levels=100, cmap='viridis',norm=norm,alpha=0.6,zorder=1)
            cbar = plt.colorbar(cf, ax=ax, pad=0.1)
        ax.set_rscale("log")
        colors = []
        steps = []

        for t in range(st,ed+1,timestep):
            self.set_index(t)
            tracerx = self.get_element("trx")
            tracery = self.get_element("try")
            tracerids = self.get_element("trid")
            R_t = []
            THETA_t = []
            for i in range(tnum):
                if(tracerids[i] in eject_mask):
                    R_t.append(tracerx[i])
                    THETA_t.append(tracery[i])
            if(t==0 and snapstep>0):
                ax.scatter(THETA_t,R_t,s=1)
            else:
                if(t>=snapc and snapstep>0):
                    random_color = (random.random(), random.random(), random.random())
                    colors.append(random_color)
                    steps.append(t)  # 保存步数
                    ax.scatter(THETA_t, R_t, color=random_color, s=2, label=f'Step {t}',zorder=3)
                    snapc = t+snapstep
            for coll in cf.collections:
                coll.remove()
            sc.set_offsets(np.c_[THETA_t, R_t])
            sc.set_zorder(4)
            if (bg != None):
                field = self.get_element(bg)[:,:,0]
                print(x1,x2)
                cf = ax.contourf(Theta, R, field.transpose(), levels=100, cmap='viridis', norm=norm, alpha=0.6, zorder=1)
                cbar.update_normal(cf)
            plt.title(f"Step {t}")
            fig.canvas.draw()
            plt.pause(0.5)

        plt.ioff()
        ax.legend(title="Recorded Steps", bbox_to_anchor=(1.1, 1.05))
        plt.show()




        #Find All Ejected Tracers, Add X Marker on final animation
    #def ShowTrajectory(self,masks)
    def GetTrajectories(self,masks,dir = "trajectories",isrestart = False,restart_step = -1):
        outdir = dir +f"/CCSN_{self.modelname}_traj"
        os.makedirs(outdir, exist_ok=True)
        restartfile = Path(outdir + "/restart")

        if(len(masks)==0):
            raise ValueError("Empty input for mask")
        if(isrestart):
            start = restart_step
            restartfile = Path(outdir+"/restart")
            if(restartfile.exists()):
                f = open(restartfile,"r")
                try:
                    step = f.readline()
                    start = int(step)
                except Exception as e:
                    print(e)
                    step = 0
        else:
            start = 0
        last_ts = 0
        for index in tqdm(
            range(start, self.n + 1),
            total=self.n + 1,
            desc="Writing trajectories",
            unit="step",
            ncols=80,

        ):
#        for index in range(start,self.n+1,1):
            self.set_index(index)


            den = self.get_element("den")
            tem = self.get_element("tem")
            xnu = self.get_element("xnu")
            alpha = self.get_element("alpha")
            phi = self.get_element("phi")
            fnu = self.get_element("fnu")
            enu = self.get_element("enu")
            dnu = self.get_element("dnu")
            firstheader = f"# t[s]      T [GK]       rho[g/cm^3]   r [km]      Ye           Le [erg/s]  Lebar [erg/s]  Ltau [erg/s] Ee [MeV]  Eebar [MeV]   Etau[MeV]\n"
            secondheader = "#----------------------------------------------------------------------------------------------------------------------------------------\n"
            print(f"Current Index:{index}")
            erg_to_MeV = 624150.913
            for tid in masks:
                self.id = tid
                r_t = self.getX()
                theta_t = self.getY()
                x,y = self.grid_finder(r_t,theta_t)
                time = self.timepoints[self.index]
                dt = time-last_ts
                tempt = tem[x][y][0]
                tempt /= 1e9

                dent = den[x][y][0]
                frac_electron = xnu[x][y][0][20]
                alpha_factor = alpha[x][y][0]
                phi_factor = phi[x][y][0]
                luminosity_nu = 4*np.pi*alpha_factor*(phi_factor**4)*(r_t**2)*fnu[x][y][0]
                enu_average = enu[x][y][0]/dnu[x][y][0]*erg_to_MeV
                enu_average = np.nan_to_num(enu_average, nan=0.0, posinf=0.0, neginf=0.0)
                file_path = f"{outdir}/winnet_{self.modelname}_{tid}"
                trajectory_file = ""
                if(not isrestart and index == start):
                    trajectory_file = open(file_path, "w")
                    trajectory_file.writelines(firstheader)
                    trajectory_file.writelines(secondheader)
                else:
                    trajectory_file = open(file_path, "a")
                trajectory_file.write(f"{time:.4e} {tempt:.4e} {dent:.4e} {r_t/1e5:.4e}  {frac_electron:.4e} {luminosity_nu[0]:.4e} {luminosity_nu[1]:.4e} {luminosity_nu[2]:.4e} {enu_average[0]:.4e} {enu_average[1]:.4e} {enu_average[2]:.4e}\n")
                trajectory_file.close()
            f = open(restartfile,"w")
            f.writelines(str(index))
            f.close()
            last_ts = time

    def GetTrajectories_new(self, masks, outdir="trajectories", isrestart=False, restart_step=-1):
        if (len(masks) == 0):
            raise ValueError("Empty input for mask")
        if (isrestart):
            start = restart_step
        restartfile = Path(outdir + "/restart")
        if (restartfile.exists()):
            f = open(restartfile, "r")
            try:
                step = f.readline()
                start = int(step)
            except Exception as e:
                print(e)
                step = 0
            else:
                start = 0
        last_ts = 0
        buffers = {tid: [] for tid in masks}
        firstheader = f"# t[s] T [GK] rho[g/cm^3] r [km] Ye Le [erg/s] Lebar [erg/s] Ltau [erg/s] Ee [MeV] Eebar [MeV] Etau[MeV]\n"
        secondheader = "#----------------------------------------------------------------------------------------------------------------------------------------\n"
        for index in tqdm(
            range(start, self.n + 1),
            total=self.n + 1,
            desc="Writing trajectories",
            unit="step",
            ncols=80,
        ):
            self.set_index(index)
            print(f"Current Index:{index}")
            erg_to_MeV = 624150.913
            den = self.get_element("den")
            tem = self.get_element("tem")
            xnu = self.get_element("xnu")
            alpha = self.get_element("alpha")
            phi = self.get_element("phi")
            fnu = self.get_element("fnu")
            enu = self.get_element("enu")
            dnu = self.get_element("dnu")
            for tid in masks:
                self.id = tid
                r_t = self.getX()
                theta_t = self.getY()
                x, y = self.grid_finder(r_t, theta_t)
                time = self.timepoints[self.index]
                dt = time - last_ts
                tempt = tem[x][y][0]
                tempt /= 1e9
                dent = den[x][y][0]
                frac_electron = xnu[x][y][0][20]
                alpha_factor = alpha[x][y][0]
                phi_factor = phi[x][y][0]
                luminosity_nu = 4 * np.pi * alpha_factor * (phi_factor ** 4) * (r_t ** 2) * fnu[x][y][0]
                enu_average = enu[x][y][0] / dnu[x][y][0] * erg_to_MeV
                enu_average = np.nan_to_num(enu_average, nan=0.0, posinf=0.0, neginf=0.0)
                line = (f"{time:.4e} {tempt:.4e} {dent:.4e} {r_t / 1e5:.4e} {frac_electron:.4e} {luminosity_nu[0]:.4e} {luminosity_nu[1]:.4e} {luminosity_nu[2]:.4e} {enu_average[0]:.4e} {enu_average[1]:.4e} {enu_average[2]:.4e}\n")
                buffers[tid].append(line)

        for tid in masks:
            with open(f"{outdir}/winnet_{self.modelname}_{tid}", "w") as fp:
                fp.write(firstheader)
                fp.write(secondheader)
                fp.writelines(buffers[tid])

    def ShowTrajectory(self,ID,vexc=False):
        trajectory = []
        for t in range(self.n+1):
            self.set_index(t)
            tracerx = self.get_element("trx")
            tracery = self.get_element("try")
            tracerids = self.get_element("trid")
            for i in range(len(tracerids)):
                id = tracerids[i]
                if(id==ID):
                    trajectory.append([tracerx[i],tracery[i]])
                    break
        trajectory = np.array(trajectory)
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
        ax.set_rlim(1e6,1e9)
        ax.set_yscale("log")
        ax.set_thetalim(0*np.pi/3, 3*np.pi/3)
        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi)
        if(vexc):
            self.set_index(self.n)
            vex = self.get_element("vex")
            xzn = self.get_element("xzn")
            yzn = self.get_element("yzn")
            vmax = np.max(np.abs(vex))
            vmin = -vmax
            R, Theta = np.meshgrid(self.gridxzn, self.gridyzn)
            contour = ax.contourf(Theta,R,vex[:,:,0].transpose(),levels=20, cmap='seismic',vmax=vmax,vmin=vmin)
            cbar = plt.colorbar(contour)


        ax.plot(trajectory[:,1],trajectory[:,0],lw=1,c="green")
        ax.scatter(trajectory[0,1],trajectory[0,0],s=50,marker="o",c="purple")
        ax.scatter(trajectory[-1,1], trajectory[-1, 0], s=50, marker="*",c="darkblue")
        fig.savefig(f"figures/{self.modelname}_{ID}_traj.png",dpi=1080)
        plt.close(fig)

    def ShowTracersAgainstNeutronExcess(self):
        xnu_electron = self.get_element("xnu")[:,:,:,20]
        etas = 1-2*xnu_electron

    def GetWinNetSeedFiles(self,dir):
        outdir = dir +f"/CCSN_{self.modelname}_seed"
        os.makedirs(outdir, exist_ok=True)
        self.set_index(self.n)
        ebind = self.get_binding_energy()[:,:,0]
        tracer_ids = self.get_element("trid")
        tracer_x = self.get_element("trx")
        tracer_y = self.get_element("try")
        tracer_num = len(tracer_ids)
        ejected = []
        self.set_index(0)
        for i in range(tracer_num):
            tx = tracer_x[i]
            ty = tracer_y[i]
            tid = tracer_ids[i]
            ix,iy = self.grid_finder(tx,ty)
            if(ebind[ix][iy]<=0):
                continue
            output = open(f"{outdir}/seed_{self.modelname}_{tid}", "w")
            Xnu = self.get_element("xnu")[ix,iy,0,:]
            print(tid, tx, ty, ix, iy)
            output.write("#    A    Z       X\n")
            for i in range(len(phycon_and_nuc_table.pc_nuc) - 3):
                Ax = int(phycon_and_nuc_table.pc_nuc[i][2])
                Zx = int(phycon_and_nuc_table.pc_nuc[i][1])
                Xx = Xnu[i]
                Xstr = f"{Xx:.3e}"
                output.write(f"{str(Ax).rjust(5)}{str(Zx).rjust(5)}{Xstr.rjust(12)}" + "\n")
            output.close()








#   3D
