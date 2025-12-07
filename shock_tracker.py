import numpy as np
from scipy.signal import find_peaks
from tqdm import tqdm

def get_shock_radii(data_handle):
    time = []
    shock_ind = []
    
    bounce = False
    prev_locs = None

    for i in tqdm(range(data_handle.ngroups)):
        data_handle.set_index(i)
        pre = data_handle.pre()
        xzn = data_handle.xzn()
        yzn = data_handle.yzn()
        zzn = data_handle.zzn()
        pre_grad = np.gradient(pre, axis=0) / pre


        # Detect the shock formation using a simple method first
        if not bounce:
            if not np.any(pre_grad < -0.38):
                continue

            #print("Bounce detected")
            bounce = True
            shock_locs = len(xzn) - 1 - np.argmin(np.flip(pre_grad, axis=0), axis=0)
            
            proms = np.empty_like(shock_locs)
            for j in range(len(yzn)):
                for k in range(len(zzn)):
                    peaks, prop = find_peaks(-pre_grad[:,j,k], prominence=0.005)
                    proms[j,k] = np.max(prop['prominences'])


        else:
            shock_locs = np.empty((len(yzn), len(zzn)), dtype=np.int16)
            proms = np.empty((len(yzn), len(zzn)))
            for j in range(len(yzn)):
                for k in range(len(zzn)):
                    peaks, prop = find_peaks(-pre_grad[:,j,k], 
                                         prominence=0.015)
                    csr = xzn[peaks]
                    cp = prop['prominences']

                    # Prevent looking too far from the previous location
                    try: 
                        if xzn[prev_locs[j,k]] < 3e8:
                            sd = 70
                        else:
                            sd = 20
                    except IndexError:
                        #print(prev_locs[j,k])
                        quit()

                    window = ((peaks - prev_locs[j,k]) < sd) & ((peaks - prev_locs[j,k]) > -sd)

                    csr = csr[window]
                    cp = cp[window]

                    # Compare each peak with the previous peak
                    psr = xzn[prev_locs[j,k]]
                    pp = prev_proms[j,k]

                    diff_fact = 1 / cp * np.abs((csr - psr) / csr)**1
                    
                    try:
                        best = np.argmin(diff_fact)
                    except ValueError:
                        #print(j,k)
                        # Generally this means the shock has vanished
                        # Either fallen onto the remnent or weakened too much somehow
                        shock_locs[j,k] = 0
                        proms[j,k] = 0
                    else:
                        shock_locs[j,k] = int(peaks[window][best])
                        proms[j,k] = cp[best]
            
            
        time.append(data_handle.time())
        shock_ind.append(shock_locs)
    
        prev_locs = shock_locs
        prev_proms = proms

    return (time, shock_ind)