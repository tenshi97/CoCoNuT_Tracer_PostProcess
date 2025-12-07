#import h5reader
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from astropy.constants import c, G, u, m_e, m_n, M_sun
from astropy import units
from tqdm import tqdm
from phycon_and_nuc_table import *
from matplotlib import colors
from scipy.signal import find_peaks


@dataclass
class ShockRadius:
    min: list
    max: list
    mean: list
    vel: list
    coarse_r: list
    coarse_v: list
    coarse_a: list
    coarse_t: list
    enc_mass: list
    raw: list
    time: list


@dataclass
class ExplosionEnergy:
    e_expl: list
    mc: list
    time: list


@dataclass
class NeutrinoLuminosity:
    enu_lum: list
    enubar_lum: list
    heavy_lum: list
    enu_mean_e: list
    enubar_mean_e: list
    heavy_mean_e: list
    time: list


@dataclass
class BHMass:
    by_mass: list
    gv_mass: list
    cv_mass: list
    time: list

@dataclass
class NSMass:
    by_mass: list
    gv_mass: list
    time: list

@dataclass
class AccRate:
    acc: list
    time: list


@dataclass
class SonicPoint:
    contour: list


def bary_to_grav(mb):
    return (-1 + np.sqrt(1 + 4 * 0.084 * mb)) / (2 * 0.084)


def get_bh_time(data_handle):
    for i in range(data_handle.ngroups):
        data_handle.set_index(i)
        if data_handle.ah_radius() > data_handle.xzn()[0]:
            return data_handle.time()


def destep(time, data):
    # Remove step-like pattern in the data by removing adjacent duplicates
    sm_data = []
    sm_time = []
    prev = None
    for (t, v) in zip(time, data):
        if v == prev:
            continue
        else:
            sm_data.append(v)
            sm_time.append(t)

        prev = v

    return (sm_time, sm_data)


def get_shock_radii(data_handle, abort_at_bounce=False, abort_at_bh=False):

    shock_r = ShockRadius([], [], [], [], [], [], [], [], [], [], [])
    bounce = False
    prev_locs = None

    for i in range(0,data_handle.ngroups,1):
        data_handle.set_index(i)
        pre = data_handle.pre()
        den = data_handle.den()
        vex = data_handle.vex()
        xzn = data_handle.xzn()
        gac = data_handle.gac()
        pre_grad = np.gradient(pre, axis=0) / pre

        if bounce and abort_at_bounce:
            break

        
        # Detect the shock formation using a simple method first
        if not bounce:
            if not np.any(pre_grad < -0.38):
                continue
                sr = np.zeros_like(data_handle.yzn())
                shock_locs = sr

            bounce = True
            shock_locs = len(xzn) - 1 - np.argmin(np.flip(pre_grad, axis=0), axis=0)
            shock_locs = shock_locs.flatten()
            sr = xzn[shock_locs]
            jvel = [0,]
            
            proms = []
            for j in range(len(data_handle.yzn())):
                peaks, prop = find_peaks(-pre_grad[:,j,0], 
                                        prominence=0.005)
                jmax = np.argmax(prop['prominences'])
                proms.append(prop['prominences'][jmax])


        else:
            shock_locs = []
            proms = []
            jvel = []
            for j in range(len(data_handle.yzn())):
                peaks, prop = find_peaks(-pre_grad[:,j,0] * xzn**2, 
                                     prominence=0.01)
                csr = xzn[peaks]
                cp = prop['prominences']

                
                # Prevent looking too far from the previous location
                if xzn[prev_locs[j]] < 3e8:
                    sd = 70
                else:
                    sd = 20
                    
                window = ((peaks - prev_locs[j]) < sd) & ((peaks - prev_locs[j]) > -sd)
                    
                csr = csr[window]
                cp = cp[window]
                
                # p_ranks = np.empty_like(cp)
                # p_ranks[cp.argsort()] = np.arange(len(cp)) + 1
                
                # Compare each peak with the previous peak
                psr = xzn[prev_locs[j]]
                pp = prev_proms[j]
                
                
                #diff_fact = np.abs((csr - psr)**1 / psr**1 / cp**1 / p_ranks**1)
                diff_fact = 1 / cp * np.abs((csr - psr) / csr)**1
                try:
                    best = np.argmin(diff_fact)
                except ValueError:
                    shock_locs.append(0)
                    proms.append(0)
                else:
                    shock_locs.append(peaks[window][best])
                    proms.append(cp[best])
                

                # While we are here get the shock velocity along the ray
                r_ind = min(shock_locs[-1] + 5, len(data_handle.xzn())-1)
                l_ind = max(shock_locs[-1] - 5, 0)
                dl = den[l_ind,j,0]
                dr = den[r_ind,j,0]
                vl = vex[l_ind,j,0]
                vr = vex[r_ind,j,0]
                pl = pre[l_ind,j,0]
                pr = pre[r_ind,j,0]
                gl = gac[l_ind,j,0]
                gr = gac[r_ind,j,0]

                csr = np.sqrt(gr * pr / dr)
                csl = np.sqrt(gl * pl / dl)
                dmach = abs(vl / csl - vr / csr)

                if r_ind >= len(data_handle.xzn()):
                    # The shock has hit the grid boundary
                    # Things probably get inaccurate now?
                    continue

                #jvel.append((dr * vr - dl * vl) / (dr - dl))
                #jvel.append(np.sqrt((gl + 1) * pr / (2 * dl)) - vl)
                jvel.append((32 + 8 * dmach * (1 + gr) + dmach**2 * (1 + gr)**2) * csr / 32)

            shock_locs = np.array(shock_locs)
            sr = xzn[shock_locs.flatten()]

# Modification
#        if data_handle.index > 400 and i%20==0 and False:
        if data_handle.index > 400 and i%50==0:
            r_grid, theta_grid = np.meshgrid(data_handle.xzn()[100:350]/1e5, data_handle.yzn())
            divnorm = colors.TwoSlopeNorm(vcenter=0)
            plt.contourf(theta_grid, r_grid, np.transpose((pre_grad[100:350,:,0]**1)),
                         levels=255, cmap='RdBu', norm=divnorm, vmin=-0.1)
            plt.colorbar()
            plt.plot(data_handle.yzn(), [r / 1e5 for r in sr])

            plt.show()

        shock_r.min.append(np.min(sr) / 1e5)
        shock_r.max.append(np.max(sr) / 1e5)
        shock_r.mean.append(np.mean(sr) / 1e5)
        shock_r.raw.append(sr / 1e5)
        shock_r.time.append(data_handle.time())

        shock_r.vel.append(np.mean(jvel) / 1e5)

        if False:
            dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
            dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
            dv = dv_r[:,np.newaxis,np.newaxis] * \
                dv_theta[np.newaxis,:,np.newaxis] * \
                2 * np.pi * data_handle.phi()**6
            
            # Make a mask for the interior of the shock and use that to sum the mass
            # of material inside the shock (i.e. enclosed mass)
            interior = [np.pad(np.array([True,] * si), (0, len(data_handle.xzn()) - si),
                            mode='constant',
                            constant_values=False) for si in shock_locs]
            interior = np.transpose(np.array(interior))
            enc_m = np.sum(data_handle.den()[:,:,0] * dv[:,:,0], where=interior)
            shock_r.enc_mass.append(enc_m)

        prev_locs = shock_locs
        prev_proms = proms

        if data_handle.alpha()[0,0,0] < 0 and abort_at_bh:
            break

        if np.mean(sr) / 1e5 > 1e3:
            pass


    # Skip chunks to get some rough velocity and (maybe) acceleration data
    step = 10   # This can be changed to get smoother or more grainy results
    # shock_r.coarse_t = shock_r.time[::step]
    # shock_r.coarse_r = shock_r.mean[::step]
    # shock_r.coarse_v = np.gradient(shock_r.coarse_r, shock_r.coarse_t)
    # shock_r.coarse_a = np.gradient(shock_r.coarse_v, shock_r.coarse_t)

    return shock_r


def plot_shock_radii(sr: ShockRadius):
    t_min, d_min = destep(sr.time, sr.min)
    t_max, d_max = destep(sr.time, sr.max)
    t_avg, d_avg = destep(sr.time, sr.mean)
    t_min = [t - sr.time[0] for t in t_min]
    t_max = [t - sr.time[0] for t in t_max]
    t_avg = [t - sr.time[0] for t in t_avg]
    fig, ax = plt.subplots(figsize=(5,3.5))
    ax.plot(t_min, d_min, label='min', c='#861657')
    ax.plot(t_max, d_max, label='max', c='#00916E')
    ax.plot(t_avg, d_avg, label='mean', c='#D2BF55')
    ax.set_xlabel("Time after bounce [s]")
    ax.set_ylabel("Shock radius [km]")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    #ax.axvline(x=get_bh_time(data_handle) - sr.time[0], c='grey', alpha=0.5, ls='--')
    ax.axvline(0 - sr.time[0], c='grey', alpha=0.5, ls='--')
    ax.legend()
    plt.minorticks_on()
    plt.show()
    fig.savefig("shock_radius_2d.pdf", bbox_inches='tight')

#modification
def get_binding_energy(data_handle):

    data_handle.set_index(data_handle.ngroups-1)

    w = 1 / np.sqrt(1 - (data_handle.vex() ** 2 + data_handle.vey() ** 2 + data_handle.vez() ** 2) / pc_cl ** 2)
    mnuc = np.zeros(21)
    mnuc[0:20] = (pc_nuc[:, 3] + pc_nuc[:, 2] * (pc_mb * pc_cl ** 2) * pc_ergmev) / pc_nuc[:, 2]
    mnuc[20] = wc_me
    e_int = data_handle.ene() * pc_gc / pc_cl ** 2
    e_nuc = np.sum([data_handle.xnu()[:, :, :, s] * mnuc[s] * pc_meverg / pc_mb for s in range(21)], axis=0)
    e_cor = (wc_mn - 8.8) * pc_meverg / pc_mb
    e_int = (e_int + pc_cl ** 2 * (1 - w)) / w + data_handle.pre() * (1 - w ** 2) / (data_handle.den() * w ** 2)
    e_int -= e_nuc * w
    e_int += e_cor * w

    e_bind = data_handle.alpha() * (w * (c.cgs.value ** 2 + e_int + data_handle.pre() / data_handle.den()) - data_handle.pre() / data_handle.den() / w) - c.cgs.value ** 2
    return e_bind

def get_exp_energy(data_handle, t_bounce, skip):

    energy = ExplosionEnergy([], [], [])

    for i in range(0, data_handle.ngroups, skip):
        data_handle.set_index(i)

        if data_handle.time() < t_bounce:
            continue

        mnuc = np.zeros(21)
        mnuc[0:20] = (pc_nuc[:,3] + pc_nuc[:,2] * (pc_mb * pc_cl**2) * pc_ergmev) / pc_nuc[:,2]
        mnuc[20] = wc_me
        
        w = 1 / np.sqrt(1 - (data_handle.vex()**2 + data_handle.vey()**2 + data_handle.vez()**2) / pc_cl**2)
        e_int = data_handle.ene() * pc_gc / pc_cl**2
        e_nuc = np.sum([data_handle.xnu()[:,:,:,s] * mnuc[s] * pc_meverg / pc_mb for s in range(21)], axis=0)
        e_cor = (wc_mn - 8.8) * pc_meverg / pc_mb
        e_int = (e_int + pc_cl**2 * (1 - w)) / w + data_handle.pre() * (1 - w**2) / (data_handle.den() * w**2)
        e_int -= e_nuc * w
        e_int += e_cor * w

        e_bind = data_handle.alpha() * (w * (c.cgs.value**2 + e_int + data_handle.pre() / data_handle.den()) - data_handle.pre() / data_handle.den() / w) - c.cgs.value**2

        dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
        dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        dv_phi = data_handle.zzr() - data_handle.zzl()

        dv = dv_r[:,np.newaxis,np.newaxis] * \
             dv_theta[np.newaxis,:,np.newaxis] * \
             2 * np.pi * data_handle.phi()**6 

        dm = np.sum(data_handle.den() * dv, axis=(1,2))
        dphimod = G.cgs.value * dm / data_handle.xzn()
        phimod = np.cumsum(dphimod[::-1])[::-1]

        e_bind += phimod[:,np.newaxis,np.newaxis]
        e_bind = np.where(np.logical_and((e_bind) > 0, data_handle.vex() > 0), e_bind, 0) * dv * data_handle.den()
        #e_bind = e_bind * dv * data_handle.den()

        if i==400 and i==500:
            r_grid, theta_grid = np.meshgrid(data_handle.xzn()[:450]/1e5, data_handle.yzn())
            divnorm = colors.TwoSlopeNorm(vcenter=0)
            plt.contourf(theta_grid, r_grid, np.transpose((e_bind[:450,:,0])),
                         levels=255, cmap='RdBu', norm=divnorm)
            plt.colorbar()
            plt.show()

        e_exp = np.sum(e_bind, axis=(1,2))

        energy.e_expl.append(e_exp)
        energy.mc.append(np.gradient(np.mean(data_handle.alpha(), axis=(1,2)), data_handle.xzn()) * c.cgs.value**2 / G.cgs.value * data_handle.xzn()**2)
        energy.time.append(data_handle.time() - t_bounce)

    return energy


def plot_exp_energy(ee: ExplosionEnergy,t_bounce,data_handle):
    fig, ax = plt.subplots(figsize=(5,3.5))
    ax.plot(ee.time, [e / 1e50 for e in ee.e_expl], c='#861657')
    ax.set_xlabel("Time after bounce [s]")
    ax.set_ylabel(r"Explosion energy [$10^{50}$ erg]")
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    plt.minorticks_on()
    ax.axvline(x=get_bh_time(data_handle) - t_bounce, c='grey', alpha=0.5, ls='--')
    #plt.yscale('log')
    plt.show()
    fig.savefig("explosion_energy_2d.pdf", bbox_inches='tight')


def get_neutrino_luminosity(data_handle, skip=1):
    nu_lum = NeutrinoLuminosity([], [], [], [], [], [], [])

    for i in range(0, data_handle.ngroups, skip):
        data_handle.set_index(i)

        # Get the area of each cell on the outer interface
        da_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        da_phi = 2 * np.pi

        da = da_theta * da_phi * data_handle.xzn()[-1]**2

        enu = np.sum(data_handle.fnu()[-1,:,0,0] * da)

        nu_lum.enu_lum.append(np.sum(data_handle.fnu()[-1,:,0,0] * da))
        nu_lum.enubar_lum.append(np.sum(data_handle.fnu()[-1,:,0,1] * da))
        nu_lum.heavy_lum.append(np.sum(data_handle.fnu()[-1,:,0,2] * da))
        nu_lum.time.append(data_handle.time())

        # Energy in the gain region
        # dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
        # dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))

        # dv = dv_r[:,np.newaxis,np.newaxis] * \
        #      dv_theta[np.newaxis,:,np.newaxis] * \
        #      2 * np.pi * data_handle.phi()**6
        
        # nu_lum.e_gain.append(np.sum(np.where((data_handle.den() < 1e11) & (data_handle.qen() > 0) & (data_handle.ah_radius() < 13000), data_handle.qen() * dv, 0)))

        nu_lum.enu_mean_e.append(np.sum((data_handle.enu()[-1,:,0,0] / data_handle.dnu()[-1,:,0,0]) * da, axis=0) / np.sum(da))
        nu_lum.enubar_mean_e.append(np.sum((data_handle.enu()[-1,:,0,1] / data_handle.dnu()[-1,:,0,1]) * da, axis=0) / np.sum(da))
        nu_lum.heavy_mean_e.append(np.sum((data_handle.enu()[-1,:,0,2] / data_handle.dnu()[-1,:,0,2]) * da, axis=0) / np.sum(da))

    return nu_lum


def plot_neutrino_luminosity(nl: NeutrinoLuminosity):
    plt.plot(nl.time, nl.enu_lum, label=r'$\mu_{e}$')
    plt.plot(nl.time, nl.enubar_lum, label=r'$\bar{\mu}_{e}$')
    plt.plot(nl.time, nl.heavy_lum, label=r'$\mu_{x}$')
    plt.xlabel("Time [s]")
    plt.ylabel("Neutrino Luminosity [erg/s]")
    plt.legend()
    plt.show()


def get_bh_mass(data_handle, t_bh, skip=1):
    bh_mass = BHMass([], [], [], [])

    #bh_is_new = True
    ah_ind_prev = 0
    ah_ind = 0
    prev_time = 0

    def mass_in_shell(s_i):
        dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
        dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        dv = dv_r[:,np.newaxis,np.newaxis] * \
             dv_theta[np.newaxis,:,np.newaxis] * \
             2 * np.pi * data_handle.phi()**6
        vsq = data_handle.vex()**2 + data_handle.vey()**2 + data_handle.vez()**2
        W = 1 / np.sqrt(1 - vsq / c.cgs.value**2)
        
        return np.sum(data_handle.den()[s_i,:,:] * W[s_i,:,0] * dv[s_i,:,:])
    

    def mass_flux(s_i):
        da_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        da_varphi = 2 * np.pi
        da = data_handle.phi()[s_i,:,0]**6 * da_theta * da_varphi * \
            data_handle.xzn()[s_i]**2
        vsq = data_handle.vex()**2 + data_handle.vey()**2 + data_handle.vez()**2
        W = 1 / np.sqrt(1 - vsq / c.cgs.value**2)
        data_handle.set_index(i-1)
        t_prev = data_handle.time()
        data_handle.set_index(i)
        dt = data_handle.time() - t_prev
        return np.sum(-da * data_handle.alpha()[s_i,:,0] * 
                      (data_handle.vex()[s_i,:,0] / 
                       data_handle.phi()[s_i,:,0]**2 - 
                       data_handle.beta1()[s_i,:] / 
                       data_handle.alpha()[s_i,:,0]) * 
                      data_handle.den()[s_i,:,0] * W[s_i,:,0]) * dt


    for i in range(0, data_handle.ngroups, skip):
        data_handle.set_index(i)
        # if data_handle.time() < t_bh:
        #     continue     

        # Start with the gravitational/ADM mass
        r = data_handle.xzn()
        phi = data_handle.phi()[:,0,0]
        dphi = np.gradient(phi, r)
        ah = list(r).index(data_handle.ah_radius())
        fact = c.cgs.value**2 / G.cgs.value / M_sun.cgs.value

        bh_mass.gv_mass.append(-2 * r[ah]**2 * dphi[ah] * fact)
        # n = -4
        # bh_mass.gv_mass.append(2 * (phi[ah] * r[ah] - (phi[ah] * r[ah] + dphi[ah] * r[ah]**2 - r[ah]) / (n + 3) - r[ah]) * fact)

        # # Now estimate the baryonic mass
        # new_mass = 0
        ah_ind = list(data_handle.xzn()).index(data_handle.ah_radius())

        # # Get the initial mass and update if ah moves
        # if ah_ind_prev == 0 or ah_ind != ah_ind_prev and False:
        #     data_handle.set_index(i-1)
        #     new_mass += sum([mass_in_shell(s_i) for s_i in range(ah_ind_prev, ah_ind + 1)])
        #     #print("initial", new_mass / M_sun.cgs.value)
        #     data_handle.set_index(i)

        # # Add the mass flux
        # new_mass += mass_flux(ah_ind)

        # new_mass /= M_sun.cgs.value
        # #print(new_mass)

        # ah_ind_prev = ah_ind

        # if len(bh_mass.by_mass) == 0:
        #     bh_mass.by_mass.append(new_mass)
        # else:
        #     bh_mass.by_mass.append(bh_mass.by_mass[-1] + new_mass)

        # Estimate mass by baryonic conservation
        dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
        dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        dv = dv_r[:,np.newaxis,np.newaxis] * \
             dv_theta[np.newaxis,:,np.newaxis] * \
             2 * np.pi * data_handle.phi()**6
        vsq = data_handle.vex()**2 + data_handle.vey()**2 + data_handle.vez()**2
        W = 1 / np.sqrt(1 - vsq / c.cgs.value**2)
        
        mass_on_grid = np.sum((data_handle.den() * W * dv)[ah_ind:])

        # Subtract off mass flux on the outer boundary since this interfers with the conservation
        mass_on_grid += np.sum((dv_theta[np.newaxis,:,np.newaxis] * 2 * np.pi * r**2 * data_handle.vex() * data_handle.den() * (data_handle.time() - prev_time))[-1])

        bh_mass.cv_mass.append(mass_on_grid)
        

        bh_mass.time.append(data_handle.time())
        prev_time = data_handle.time()

    return bh_mass


def get_ns_mass(data_handle, t_bh, skip=1):
    ns_mass = NSMass([], [], [])

    for i in range(0, data_handle.ngroups, skip):
        data_handle.set_index(i)
        if data_handle.time() >= t_bh:
            # print("AH", data_handle.ah_radius())
            # input()
            continue

        dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
        dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
        dv_phi = data_handle.zzr() - data_handle.zzl()

        dv = dv_r[:,np.newaxis,np.newaxis] * \
            dv_theta[np.newaxis,:,np.newaxis] * \
            2 * np.pi

        vsq = data_handle.vex()**2 + data_handle.vey()**2 + data_handle.vez()**2
        W = 1 / np.sqrt(1 - vsq / c.cgs.value**2)

        dm = data_handle.den() * W * data_handle.phi()**6

        enclosed = np.sum(np.where(dm > 1e11, dm * dv, 0))
        ns_mass.by_mass.append(enclosed / M_sun.cgs.value)

        r = data_handle.xzn()
        phi = data_handle.phi()[:,0,0]
        dphi = np.gradient(phi, r)
        ns_int = np.where(np.mean(data_handle.den(), axis=(1,2)) > 1e11)[0]
        ns_r_index = 0 if len(ns_int) == 0 else np.max(ns_int) + 0
        fact = c.cgs.value**2 / G.cgs.value / M_sun.cgs.value

        # print(r[ns_r_index])
        ns_mass.gv_mass.append(-2 * r[ns_r_index]**2 * dphi[ns_r_index] * fact)
        n = -4
        # ns_mass.gv_mass.append(2 * (phi[ns_r_index] * r[ns_r_index] - (phi[ns_r_index] * r[ns_r_index] + dphi[ns_r_index] * r[ns_r_index]**2 - r[ns_r_index]) / (n + 3) - r[ns_r_index]) * fact)


        ns_mass.time.append(data_handle.time())

    return ns_mass


def get_acc_rate(data_handle, skip=1):
    acc_rate = AccRate([], [])
    y = data_handle

    da_theta = abs(np.cos(y.yzl())-np.cos(y.yzr()))
    da_phi = 2 * np.pi
    
    da = da_theta * da_phi

    for i in range(0, data_handle.ngroups, skip):
        data_handle.set_index(i)
    
        flux_index = np.searchsorted(y.xzn(), 1e7)
        w = 1 / np.sqrt(1 - (y.vex()[flux_index,:,0]**2 + y.vey()[flux_index,:,0]**2 + y.vez()[flux_index,:,0]**2) / c.cgs.value**2)
        alpha = y.alpha()[flux_index,:,0]
        vr = y.vex()[flux_index,:,0]
        br = y.beta1()[flux_index,:,0]
        den = y.den()[flux_index,:,0]
        phi = y.phi()[flux_index,:,0]
        mflux = alpha * (vr / phi**2 - br / alpha) * den * w * phi**6 * da * y.xzn()[flux_index]**2
        acc_rate.acc.append(np.sum(mflux))
        acc_rate.time.append(y.time())

    return acc_rate


def get_sonic_point(data_handle):
    # Determine the sonic point (as a contour) at the set timestep

    sonic = SonicPoint([])

    dv_r = 1./3. * (data_handle.xzr()**3-data_handle.xzl()**3)
    dv_theta = abs(np.cos(data_handle.yzl())-np.cos(data_handle.yzr()))
    dv_phi = data_handle.zzr() - data_handle.zzl()

    dv = dv_r[:,np.newaxis,np.newaxis] * \
        dv_theta[np.newaxis,:,np.newaxis] * \
        2 * np.pi * data_handle.phi()**6
        
    # Get the velocity
    vex = data_handle.vex()

    # Get the local speed of sound
    cs = np.sqrt(data_handle.gac() * data_handle.pre() / data_handle.den())

    # Find the line of intersection with the zero plane of the difference
    diff = vex[:,:,0] + cs[:,:,0]
    diff = np.average(diff, axis=1, weights=dv[:,:,0])
    r1 = diff[:-1]
    r2 = diff[1:]
    sflip = np.where(r1 * r2 < 0)[0]
    if len(sflip > 0):
        sonic.contour = [sflip[-1] for i in range(128)]
    else:
        sonic.contour = [0 for i in range(128)]
    # for j in range(len(data_handle.yzn())):
    #     r1 = diff[:-1,j]
    #     r2 = diff[1:,j]
    #     sflip = np.where(r1 * r2 < 0)[0]
    #     sonic.contour.append(sflip[-1])
    return sonic
        

def test():
    j = 100
    k = 0
    bounce = False
    prev_locs = None

    sr = []
    time = []
    v_shock = []
    v0 = []
    vp1 = []
    vm1 = []
    vp5 = []
    vm5 = []
    cs = []
    us = []
    us1 = []
    us2 = []

    Rp = []
    Rm = []
    Lw = []

    dmach = []

    for i in tqdm(range(y.ngroups)):
        y.set_index(i)
        pre = y.pre()
        vex = y.vex()
        xzn = y.xzn()
        pre_grad = np.gradient(pre, axis=0) / pre

        # Detect the shock formation using a simple method first
        if not bounce:
            if not np.any(pre_grad < -0.33):
                continue
                si = np.zeros_like(y.yzn())
                shock_locs = si
            else:
                bounce = True
                shock_locs = len(xzn) - 1 - np.argmax(np.flip(pre_grad < -0.33, axis=0), axis=0)
                si = xzn[shock_locs.flatten()]

        else:
            # Make a masked array of pre centered on the previous shock
            # location in each ray direction
            mask_band = np.ones_like(pre_grad)
            for jv in range(len(y.yzn())):
                s = prev_locs[j][0]
                w = 15
                mask_band[s-w:s+w,jv] = 0

            masked_pre_grad = np.ma.masked_array(pre_grad, mask=mask_band)
            shock_locs = np.argmin(masked_pre_grad, axis=0)
            si = shock_locs.flatten()[j]

            sr.append(xzn[si])
            time.append(y.time())

            v0.append(y.vex()[si,j,k])
            vp1.append(y.vex()[si+1,j,k])
            vm1.append(y.vex()[si-1,j,k])
            vp5.append(y.vex()[si+5,j,k])
            vm5.append(y.vex()[si-5,j,k])

            cs.append(np.sqrt(y.gac()[si+5,j,k] * y.pre()[si+5,j,k] / y.den()[si+5,j,k]))

            us.append(y.vex()[si+5,j,k] + np.sqrt(y.gac()[si+5,j,k] * y.pre()[si+5,j,k] / y.den()[si+5,j,k]) *
                      np.sqrt(1 + (y.gac()[si+5,j,k] + 1) / (2 * y.gac()[si+5,j,k]) *
                      (y.pre()[si-5,j,k] / y.pre()[si+5,j,k] - 1)))

            us1.append((y.den()[si-5,j,k] * y.vex()[si-5,j,k] -
                        y.den()[si+5,j,k] * y.vex()[si+5,j,k]) /
                       (y.den()[si-5,j,k] - y.den()[si+5,j,k]))
            us2.append((y.den()[si-5,j,k] * y.vex()[si-5,j,k]**2 + y.pre()[si-5,j,k] -
                        y.den()[si+5,j,k] * y.vex()[si+5,j,k]**2 - y.pre()[si+5,j,k]) /
                       (y.den()[si-5,j,k] * y.vex()[si-5,j,k] -
                        y.den()[si+5,j,k] * y.vex()[si+5,j,k]))

            den = y.den()[:,j,k]
            pre = y.pre()[:,j,k]
            gac = y.gac()[:,j,k]
            vex = y.vex()[:,j,k]
            r = y.xzn()

            den_c1 = (np.log10(den[si+5]) - np.log10(den[si+10])) / (np.log10(r[si+5]) - np.log10(r[si+10]))
            den_c2 = (np.log10(den[si+5]) * np.log10(r[si+10]) - np.log10(den[si+10]) * np.log10(r[si+5])) / \
                     (np.log10(r[si+10]) - np.log10(r[si+5]))

            den_pl = 10**(den_c1 * np.log10(r[si-5]) + den_c2)


            pre_c1 = (np.log10(pre[si+5]) - np.log10(pre[si+10])) / (np.log10(r[si+5]) - np.log10(r[si+10]))
            pre_c2 = (np.log10(pre[si+5]) * np.log10(r[si+10]) - np.log10(pre[si+10]) * np.log10(r[si+5])) / \
                     (np.log10(r[si+10]) - np.log10(r[si+5]))

            pre_pl = 10**(pre_c1 * np.log10(r[si-5]) + pre_c2)

            dJp = vex[si-5] + 2 * np.sqrt(gac[si-5] * pre[si-5] / den[si-5]) / (gac[si-5] - 1) - \
                  vex[si+5] - 2 * np.sqrt(gac[si+5] * pre_pl / den_pl) / (gac[si+5] - 1)
            Rp.append(np.sqrt(den_pl * np.sqrt(gac[si+5] * pre_pl / den_pl) *
                                                  4 * np.pi * r[si]**2) * dJp)

            dJm = vex[si-5] - 2 * np.sqrt(gac[si-5] * pre[si-5] / den[si-5]) / (gac[si-5] - 1) - \
                  vex[si+5] + 2 * np.sqrt(gac[si+5] * pre_pl / den_pl) / (gac[si+5] - 1)
            Rm.append(np.sqrt(den_pl * np.sqrt(gac[si+5] * pre_pl / den_pl) *
                                                  4 * np.pi * r[si]**2) * dJm)

            Lw.append(1 / 4 * (Rp[-1]**2 - Rm[-1]**2))

            dmach.append(vex[si-5] / np.sqrt(gac[si-5] * pre[si-5] / den[si-5]) -
                         vex[si+5] / np.sqrt(gac[si+5] * pre[si+5] / den[si+5]))

        prev_locs = shock_locs

    dtime, dsr = destep(time, sr)
    cs = np.array(cs)
    vp5 = np.array(vp5)
    vm5 = np.array(vm5)
    # plt.plot(time, dmach, label='dmach')
    plt.plot(time, Rp, label='Rp')
    plt.plot(time, Rm, label='Rm')
    #plt.plot(time, Lw, label="Lw")
    plt.axhline(y=0)
    # plt.plot(dtime, np.gradient(dsr, dtime), label="gradient")
    # plt.plot(time, v0, label="v0")
    # plt.plot(time, vp1, label="vp1")
    # plt.plot(time, vm1, label="vm1")
    # plt.plot(time, cs - vp5, label="vp5")
    # plt.plot(time, cs - vm5, label="vm5")
    # plt.plot(time, cs, label="cs")
    # plt.plot(time, us, label='us')
    # plt.plot(time, us1, label='us1')
    # plt.plot(time, us2, label='us2')
    plt.legend()
    plt.show()


# shock_r = get_shock_radii()
# #test()
# plot_shock_radii(shock_r)
# #
# exp_energy = get_exp_energy(shock_r.time[0] + 0.001)
# plot_exp_energy(exp_energy, shock_r.time[0] + 0.001)
#
# nu_lum = get_neutrino_luminosity()
# plot_neutrino_luminosity(nu_lum)
