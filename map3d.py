#25.18

import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt

# Use this function to map 2D->2D or 3D->3D
def to_cart (R3, Mu, Phi, flipyy = 0):
    R = np.cbrt(3*R3)
    Theta = np.arccos(Mu)
    x = R * np.sin (Theta) * np.cos (Phi)
    y = R * np.sin (Theta) * np.sin (Phi)
    z = R * np.cos (Theta)
    if flipyy == 0:
#        print('flipyy',0)
        return  x,y,z
    else:
#        print('flipyy',1)
        return -x,z,y        

#Use this_function to map 3D->2D
def to_cart_32D (R, Theta, Phi):
    n1, n2, n3 = (1,0,0)
    xi = R * np.sin (Theta) * np.cos (Phi)
    yi = R * np.sin (Theta) * np.sin (Phi)
    zi = R * np.cos (Theta)

    x = np.linalg.norm(np.cross((xi,yi,zi),(n1,n2,n3)))
    y = 0.0
    z = np.dot(a,b)

    return x,y,z

def map3d(model, rl, rr, thl, thr, phil, phir, x0, x1, z0, z1, nx, ny, nz, tol):

    print('Starting Mapping...')

    rho = model.den()
    xnu = model.xnu()

    # print(np.shape(rho))
    # exit()

    thl = np.concatenate((thl, thl))
    thr = np.concatenate((thr, thr))

    nr  = np.size(rl)
    nth = np.size(thl)
    nph = np.size(phil)
    nnuc = np.size(xnu[0,0,0,:])

    yzl = np.concatenate((model.yzl(), model.yzl()))
    yzr = np.concatenate((model.yzr(), model.yzr()))
    yzn = np.concatenate((model.yzn(), model.yzn()))

    yywt = np.concatenate((model.yinyang_weight(), model.yinyang_weight()))

    xigmap = 0.0

# calculate explosion properties
    dv_r = 1./3. * (model.xzr()**3-model.xzl()**3)
    dv_theta = abs(np.cos(yzl)-np.cos(yzr))
    dv_phi = model.zzr() - model.zzl()
    # CHECK YINYANG GRID AND WEIGHT
    dv = dv_r[:,np.newaxis,np.newaxis]*dv_theta[np.newaxis,:,np.newaxis]*dv_phi[np.newaxis,np.newaxis,:]*yywt[np.newaxis,:,:]
    print(np.shape(dv))
    # exit()
    xig=np.sum(xnu[:,:,:,14:nnuc-1],axis=3)
    print('Ejecta mass',np.sum(dv*rho*(1.0-xnu[:,:,:,2])/1.99e33),'M_sun')
    print('Explosion energy',np.sum(dv*rho*model.ene()),'erg')
    print('Nickel/IG mass',np.sum(dv*rho*xig/1.99e33),'M_sum')
    print('He mass',np.sum(dv*rho*xnu[:,:,:,4]/1.99e33),'M_sun')
    print('C mass',np.sum(dv*rho*xnu[:,:,:,5]/1.99e33),'M_sun')
    print('O mass',np.sum(dv*rho*xnu[:,:,:,7]/1.99e33),'M_sun')
    print('Ne mass',np.sum(dv*rho*xnu[:,:,:,8]/1.99e33),'M_sun')
    print('Mg mass',np.sum(dv*rho*xnu[:,:,:,9]/1.99e33),'M_sun')
    print('Si mass',np.sum(dv*rho*xnu[:,:,:,10]/1.99e33),'M_sun')
    print('S mass',np.sum(dv*rho*xnu[:,:,:,11]/1.99e33),'M_sun')
    print('Ar mass',np.sum(dv*rho*xnu[:,:,:,12]/1.99e33),'M_sun')
    print('Ca mass',np.sum(dv*rho*xnu[:,:,:,13]/1.99e33),'M_sun')
    print('Ti mass',np.sum(dv*rho*xnu[:,:,:,14]/1.99e33),'M_sun')

    dm = np.sum (rho[:,:,:]*dv,axis=(1,2))
    xhe = np.sum(rho[:,:,:]*xnu[:,:,:,4]*dv,axis=(1,2)) / dm
    xc  = np.sum(rho[:,:,:]*xnu[:,:,:,5]*dv,axis=(1,2)) / dm
    xo  = np.sum(rho[:,:,:]*xnu[:,:,:,7]*dv,axis=(1,2)) / dm
    xne = np.sum(rho[:,:,:]*xnu[:,:,:,8]*dv,axis=(1,2)) / dm
    xmg = np.sum(rho[:,:,:]*xnu[:,:,:,9]*dv,axis=(1,2)) / dm
    xsi = np.sum(rho[:,:,:]*xnu[:,:,:,10]*dv,axis=(1,2)) / dm
    xs  = np.sum(rho[:,:,:]*xnu[:,:,:,11]*dv,axis=(1,2)) / dm
    xar = np.sum(rho[:,:,:]*xnu[:,:,:,12]*dv,axis=(1,2)) / dm
    xca = np.sum(rho[:,:,:]*xnu[:,:,:,13]*dv,axis=(1,2)) / dm
    xti = np.sum(rho[:,:,:]*xnu[:,:,:,14]*dv,axis=(1,2)) / dm
    xigav = np.sum(rho[:,:,:]*xig          *dv,axis=(1,2)) / dm

    vmid = np.sum(rho[:,:,:]*model.vex()      *dv,axis=(1,2)) / dm

    plt.close('all')
    plt.clf()

    f, ax = plt.subplots(1, sharex=True, figsize=(6,7))
    f.subplots_adjust(hspace=0.0,left=0.18,bottom=0.1,top=0.96,right=0.9)

    qc1,=ax.plot(model.xzn()/1e5,xhe,label='He')
    qc2,=ax.plot(model.xzn()/1e5,xc ,label='C')
    qc3,=ax.plot(model.xzn()/1e5,xo ,label='O')
    qc4,=ax.plot(model.xzn()/1e5,xne,label='Ne')
    qc5,=ax.plot(model.xzn()/1e5,xmg,label='Mg')
    qc6,=ax.plot(model.xzn()/1e5,xsi,label='Si')
    qc7,=ax.plot(model.xzn()/1e5,xca,label='Ca')
    qc8,=ax.plot(model.xzn()/1e5,xigav,label='Fe group')

    ax.legend(handles=[qc1,qc2,qc3,qc4,qc5,qc6,qc7,qc8],fontsize=11,handlelength=1.4,columnspacing=0.4,loc='upper left'\
,ncol=4)
#    ax.set_xlabel(r'$v_r\ [\mathrm{km}\, \mathrm{s}^{-1}]$')
    ax.set_xlabel(r'$r\ [\mathrm{km}]$')
    ax.set_ylabel(r'$\Delta M_i\ [M_\odot]$')

    #plt.savefig('composition.pdf')
    plt.close('all')

#    for i in range(nnuc):
#        xnu[:,:,:,i]=xnu[:,:,:,i]*rho[:,:,:]

    mul = np.cos(thl)
    mur = np.cos(thr)

    r3l = rl ** 3 / 3.0
    r3r = rr ** 3 / 3.0

    #cell centre coordinates in PROMETHEUS
    r3c  = 0.5 * (r3l+r3r)
    muc  = 0.5 * (mul+mur)
    thc  = 0.5 * (thl+thr)
    phic = 0.5 * (phil+phir)
    dphi = 2.0*np.pi/nph

    if (nph == 1):
        phil[:] = -0.0
        phir[:] =  0.0
        phic[:] =  0.0

    #cell spcaing in ARTIS
    y0 = x0
    y1 = x1

    dx = (x1-x0) / float(nx)
    dy = (y1-y0) / float(ny) # same as dx if we're in 3D, otherwise this should not be used anyway
    dz = (z1-z0) / float(nz)


    rhonew = np.zeros([nx,ny,nz])
    xnunew = np.zeros([nx,ny,nz,nnuc])
    xnunew [:,:,:,4]= 1e-50

    flipyy = 0
    
    # Map each cell
    # CHECK YINYANG DATA LAYOUT - ONE GRID PATCHED?

#    for kk in range(55,57):
    for kk in range(nph):

        km1 = max(kk-1,0)
        kp1 = min(kk+1,nph-1)

        for jj in range (nth):

            jm1 = max(jj-1,0)
            jp1 = min(jj+1,nth-1)
            
            if nph > 1:
                if jj == nth/2-1: jp1 = jj
                if jj == nth/2: jm1 = jj
                if jj >= nth/2:
                    flipyy = 1
                else:
                    flipyy = 0

#            print (jj,kk,yywt[jj,kk])
            if yywt[jj,kk]==0: continue                    

            for ii in range (nr):

                im1 = max(ii-1,0)
                ip1 = min(ii+1,nr-1)

# minmod slopes in each direction
                sr = (rho[ip1,jj ,kk ]-rho[ii , jj,kk ]) / (r3c[ip1] - r3c [ii] + 1e-80)
                sl = (rho[ii ,jj ,kk ]-rho[im1, jj,kk ]) / (r3c[ii ] - r3c [im1] + 1e-80)
                s1rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (rho[ii ,jp1,kk ]-rho[ii ,jj ,kk]) / (muc[jp1] - muc [jj] + 1e-80)
                sl = (rho[ii ,jj ,kk ]-rho[ii ,jm1,kk]) / (muc[jj ] - muc [jm1] + 1e-80)
                s2rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (rho[ii ,jj ,kp1]-rho[ii , jj,kk ]) / dphi
                sl = (rho[ii ,jj ,kk ]-rho[ii , jj,km1]) / dphi
                s3rho = min(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ip1,jj ,kk ,:]-xnu[ii , jj,kk ,:]) / (r3c[ip1] - r3c [ii] + 1e-80)
                sl = (xnu[ii ,jj ,kk ,:]-xnu[im1, jj,kk ,:]) / (r3c[ii ] - r3c [im1] + 1e-80)
                s1xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ii ,jp1,kk ,:]-xnu[ii ,jj ,kk,:]) / (muc[jp1] - muc [jj]  + 1e-80)
                sl = (xnu[ii ,jj ,kk ,:]-xnu[ii ,jm1,kk,:]) / (muc[jj ] - muc [jm1] + 1e-80)
                s2xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                sr = (xnu[ii ,jj ,kp1,:]-xnu[ii , jj,kk ,:]) / dphi
                sl = (xnu[ii ,jj ,kk ,:]-xnu[ii , jj,km1,:]) / dphi
                s3xnu = np.minimum(np.abs(sr),np.abs(sl))*np.sign(sr)*np.heaviside(sl*sr,0.)

                rho0 = rho[ii,jj,kk]
                xnu0 = xnu[ii,jj,kk,:]

# additional slope limiting for mass fraction
                s1xnu = s1xnu / (1.0 + np.abs(s1rho * (r3r[ii] - r3l[ii]) / rho0) / 12.0)
                s2xnu = s2xnu / (1.0 + np.abs(s2rho * (mur[jj] - mul[jj]) / rho0) / 12.0)
                s3xnu = s3xnu / (1.0 + np.abs(s3rho * (phir[kk] - phil[kk]) / rho0) / 12.0)

                delta_xnu0 = (s1rho * s1xnu * (r3r[ii] - r3l[ii]) ** 2 \
                    + s2rho * s2xnu * np.abs(mur[jj] - mul[jj]) ** 2 \
                    + s3rho * s3xnu * np.abs(phir[kk] - phil[kk]) ** 2) / 12.0
                xnu0 = xnu0 - delta_xnu0 / rho0
#                sys.exit()

                r3c0  = r3c[ii]
                muc0  = muc[jj]
                phic0 = phic[kk]

#debug
#                s1rho = 0.0
#                s2rho = 0.0
#                s3rho = 0.0

                def target_indices (vec):
                    x, y, z = vec
                    nonlocal x0, y0, z0, dx, dy, dz
                    m = int(np.floor ((x-x0) / dx))
                    n = int(np.floor ((y-x0) / dy))
                    o = int(np.floor ((z-z0) / dz))
                    return m, n, o

                def map_cube(Xl, Xr, Yl, Yr, Zl, Zr, tol):

                    nonlocal rho,rhonew,xnu,xnunew,rho0,xnu0,s1rho,s2rho,s3rho,s1xnu,s2xnu,s3xnu,nx,ny,nz,nsub,nsubmax,dx,dy,dz,r3c0,muc0,phic0,xigmap,flipyy

#                    print('map_cube',Xl,Xr,Yl,Yr,Zl,Zr)
# maybe doing all of the eight corners is overkill...
                    c1 = target_indices (to_cart(Xl, Yl, Zl, flipyy=flipyy))
                    c2 = target_indices (to_cart(Xl, Yl, Zr, flipyy=flipyy))
                    c3 = target_indices (to_cart(Xl, Yr, Zl, flipyy=flipyy))
                    c4 = target_indices (to_cart(Xl, Yr, Zr, flipyy=flipyy))
                    c5 = target_indices (to_cart(Xr, Yl, Zl, flipyy=flipyy))
                    c6 = target_indices (to_cart(Xr, Yl, Zr, flipyy=flipyy))
                    c7 = target_indices (to_cart(Xr, Yr, Zl, flipyy=flipyy))
                    c8 = target_indices (to_cart(Xr, Yr, Zr, flipyy=flipyy))

# Check whether we're overlapping with more than one target cell
                    m0 = min(c1[0], c2[0], c3[0], c4[0], c5[0], c6[0], c7[0], c8[0])
                    m1 = max(c1[0], c2[0], c3[0], c4[0], c5[0], c6[0], c7[0], c8[0])
                    n0 = min(c1[1], c2[1], c3[1], c4[1], c5[1], c6[1], c7[1], c8[1])
                    n1 = max(c1[1], c2[1], c3[1], c4[1], c5[1], c6[1], c7[1], c8[1])
                    o0 = min(c1[2], c2[2], c3[2], c4[2], c5[2], c6[2], c7[2], c8[2])
                    o1 = max(c1[2], c2[2], c3[2], c4[2], c5[2], c6[2], c7[2], c8[2])

                    Xc = 0.5 * (Xr+Xl)
                    Yc = 0.5 * (Yr+Yl)
                    Zc = 0.5 * (Zr+Zl)
                    mc, nc, oc = target_indices (to_cart(Xc, Yc, Zc, flipyy=flipyy))

                    if (nph==1): dv = abs(Xr-Xl) * abs(Yr-Yl) * dphi
                    # CHECK YINYANG WEIGHT FACTOR NEEDED?
                    if (nph>1): dv = abs(Xr-Xl) * abs(Yr-Yl) * abs(Zr-Zl) * yywt[jj,kk]
#                    print('dv',abs(Xr-Xl),abs(Yr-Yl), abs(Zr-Zl), yywt[jj,kk],flipyy)
                    if (ny==1):
                        dv0 = np.pi*(2*mc+1)*dx**2*dz
                    if (ny>1): dv0 = float(dx*dy*dz)

#                    if (nsub==1):
#                        print ('map_cube2',m0,m1,n0,n1,o0,o1,mc,nc,oc,dv,tol*dv0,nsub)
#                        print (np.cbrt(Xc*3),Yc,Zc,abs(Xr-Xl),abs(Yr-Yl),nph)
#                    print ('map_cube2',m0,m1,n0,n1,o0,o1,Xr-Xl,Yr-Yl,Zr-Zl)
                    if ((max(m0,m1)>=nx) | (max(n0,n1)>=ny) | (max(o0,o1)>=nz)): return
                    if ((min(m0,m1)<0) | (min(n0,n1)<0) | (min(o0,o1)<0)): return
                        
                    if ((m0==m1) & (n0==n1) & (o0==o1)) | (dv < tol*dv0):
                        #map entire cube to one destination cell
#                        print ('converged',mc,nc,oc,(Xc*3)**0.33333,Yc,Zc)
#                        if ((rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0))<0):
#                            print('rho<0',rho0,s1rho,s2rho,s3rho,Xc,Yc,Zc)
#                            sys.exit()
                        rhonew[mc,nc,oc]   += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * dv/dv0
#                        print (rho0,s1rho,dv,dv0,rhonew[mc,nc,oc])
#                        print (mc,nc,oc,rhonew[mc,nc,oc])
# 1st order
#                        xnunew[mc,nc,oc,:] += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * xnu0 * dv/dv0
# 2nd order
                        xnunew[mc,nc,oc,:] += (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * \
                            (xnu0 + s1xnu * (Xc-r3c0) + s2xnu * (Yc-muc0) + s3xnu * (Zc-phic0)) * \
                            dv/dv0
##                        xigmap +=  (rho0 + s1rho * (Xc-r3c0) + s2rho * (Yc-muc0) + s3rho * (Zc-phic0)) * np.sum(xnu0[14:nnuc-1]) * dv
                        return
                    else:
                        #divide into eight subcubes and map them
                        nsub = nsub + 1
                        nsubmax = max(nsub,nsubmax)
                        map_cube(Xl, Xc, Yl, Yc, Zl, Zc, tol)
                        map_cube(Xl, Xc, Yc, Yr, Zl, Zc, tol)
                        map_cube(Xc, Xr, Yl, Yc, Zl, Zc, tol)
                        map_cube(Xc, Xr, Yc, Yr, Zl, Zc, tol)
                        # print(nph)
                        if (nph > 1):
                            map_cube(Xl, Xc, Yl, Yc, Zc, Zr, tol)
                            map_cube(Xl, Xc, Yc, Yr, Zc, Zr, tol)
                            map_cube(Xc, Xr, Yl, Yc, Zc, Zr, tol)
                            map_cube(Xc, Xr, Yc, Yr, Zc, Zr, tol)

                        nsub = nsub - 1

                nsub  = 1
                nsubmax = 1
                map_cube(r3l[ii], r3r[ii], mul[jj], mur[jj], phil[kk], phir[kk], tol)

#                print ('Mapped zone',ii,jj,kk,nsubmax)
#            print('y zone = ', jj)
        print('z zone = ', kk)

#    rhonew = np.maximum(rhonew, 1e-50)
    rhonew = np.maximum(rhonew, 1e-5 * np.min(rho))
    for i in range(nnuc):
        xnunew[:,:,:,i]=xnunew[:,:,:,i]/rhonew[:,:,:]

    x = np.arange(x0,x1,dx)+0.5*dx
    if ny > 1:
        y = np.arange(y0,y1,dy)+0.5*dy
    z = np.arange(z0,z1,dz)+0.5*dz
    # r = np.sqrt(np.outer(x**2,z**0) + np.outer(x**0,z**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    # GENERALISE R TO 3D

    xignew=np.sum(xnunew[:,:,:,14:nnuc-1],axis=3)

    plt.contourf(z,x,np.log10(rhonew[:,0,:]),50)
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$z\ [\mathrm{cm}]$')
    plt.savefig('density.pdf')

    plt.contourf(z,x,xignew[:,0,:],50)
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$z\ [\mathrm{cm}]$')
    plt.savefig('iron_group.pdf')

    plt.contourf(z,x,xnunew[:,0,:,7],50)
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$z\ [\mathrm{cm}]$')
    plt.savefig('oxygen.pdf')

    plt.contourf(z,x,xnunew[:,0,:,4],50)
    plt.xlabel(r'$x\ [\mathrm{cm}]$')
    plt.ylabel(r'$z\ [\mathrm{cm}]$')
    plt.savefig('helium.pdf')

# Density and Detailed composition for ARTIS
# Assuming slightly proton-rich ejecta

    xartis = np.zeros([nx,ny,nz,1+30])

    xartis[:,:,:,0] = rhonew [:,:,:]

    xartis[:,:,:,1]  = xnunew [:,:,:,0] + xnunew [:,:,:,1] # n, p -> H
    xartis[:,:,:,2]  = xnunew [:,:,:,4] # He
    xartis[:,:,:,3]  = 0.0 # Li
    xartis[:,:,:,4]  = 0.0 # Be
    xartis[:,:,:,5]  = 0.0 # B
    xartis[:,:,:,6]  = xnunew [:,:,:,5] # C
    xartis[:,:,:,7]  = xnunew [:,:,:,6] # N
    xartis[:,:,:,8]  = xnunew [:,:,:,7] # O
    xartis[:,:,:,10] = xnunew [:,:,:,8] # Ne
    xartis[:,:,:,12] = xnunew [:,:,:,9] # Mg
    xartis[:,:,:,14] = xnunew [:,:,:,10] # Si
    xartis[:,:,:,16] = xnunew [:,:,:,11] # S
    xartis[:,:,:,18] = xnunew [:,:,:,12] # Ar
    xartis[:,:,:,20] = xnunew [:,:,:,13] + xignew * 3e-3 # Ca
    xartis[:,:,:,22] = xignew * 1.5e-3 # Ti
    xartis[:,:,:,24] = xignew * 2.5e-3 # Cr
    xartis[:,:,:,26] = xignew * 2.0e-3 # Fe52

    xdec = 0.5**(model.time()/(6.077*86400.0))
    xcr48 = xignew * 2.5e-3 # Cr
    xfe52 = xignew * 2.0e-3 # Fe52
    xco56 = xignew * 0.99392 * (1.0-xdec)
    xni56 = xignew * 0.99392 * xdec

    xartis[:,:,:,27] = xco56 # Co56
    xartis[:,:,:,28] = xni56 + xignew *  8e-5 # Ni56, Ni58

# Normalise mass fractions
    xsum = np.sum(xartis[:,:,:,1:],3)
    for i in range (30):
        xartis[:,:,:,i+1] = xartis[:,:,:,i+1] / xsum
    xcr48 = xcr48 / xsum
    xfe52 = xfe52 / xsum
    xco56 = xco56 / xsum
    xni56 = xni56 / xsum
    xignew = xignew / xsum

    print('Writing input files for ARTIS.')

    f = open('model.txt', 'tw')
    g = open('abundances.txt', 'tw')


# CHECK WITH STUART AND FINN HOW THEY SET UP 3D INPUT DATA
    ij = 0
    f.write ("%d %d \n" % (nx, nz))
    time = model.time()
    f.write ("%16.8e \n" % (time/86400.))
    f.write ("%16.8e \n" % (max(x1/time,z1/time)))
    for kk in range(nz):
        for jj in range(ny):
            for ii in range(nx):
                ij += 1
                f.write ("%d %16.8e %16.8e %16.8e %16.8e \n %16.8e %16.8e %16.8e %16.8e %16.8e \n" % (ij, x[ii], y[jj], z[kk], xartis[ii,jj,kk,0], xignew[ii,jj,kk], xni56[ii,jj,kk], xco56[ii,jj,kk], xfe52[ii,jj,kk], xcr48[ii,jj,kk]))
                out = " ".join(("%16.8e " %dat) for dat in tuple(xartis[ii,jj,kk,1:]))
            g.write(("%d " % ij) + out + " \n")

    f.closed
    g.closed

    with h5py.File('s3.5_mapped.h5', 'w') as hf:
        hf.create_dataset('model', data=xartis)


# Check mass of mapped model
    dv0 = dx*dy*dz
    print("Mass of mapped model:",np.sum(xartis[:,0,:,0]*dv0)/1.99e33,"M_sun")
    print("Iron group:          ",np.sum(xartis[:,0,:,0]*xignew[:,0,:]*dv0)/1.99e33,"M_sun")
#    print("Iron group:          ",np.sum(xartis[:,0,:,0]*xignew[:,0,:]*dv0)/1.99e33,xigmap/1.99e33,"M_sun")


    print('Done!')
    return rhonew,xnunew,x,y,z,r
