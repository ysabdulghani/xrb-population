from xspec import *
import argparse
import os
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool, cpu_count, TimeoutError
from tqdm import tqdm
from xspec_simulations import *
import random
import urllib
import h5py
import subprocess

def idx_of_value_from_grid(grid,value,atol=1e-08,verbose=False):
    """
    Finds the index of a specified value in a grid array with a specified absolute tolerance.
    If no exact match is found, the tolerance is increased iteratively to find a close match.

    """
    index, = np.where(np.isclose(grid,value,rtol=1e-05, atol=atol))
    atol_tmp = atol
    while len(index) == 0:
        atol_tmp *= 10
        index, = np.where(np.isclose(grid,value,rtol=1e-05, atol=atol_tmp))
    if len(index) != 1:
        index = index + (index[-1]-index[0]) // 2
        if verbose==True:
            print("Warning: found more than grid value (for GR correction) that are equally close to the user specified value of %s. Taking the median value %s as the one to be closest." % (value,grid[index[0]]))
    return index

def get_total_correction_GR_and_Rin_Rg_ratio(inc,a,limb_dark=True,verbose=True):
    """
    Calculates the total correction factor and the Rin/Rg ratio based on the provided inclination angle and spin parameter.
    The correction factors are derived using General Relativity (GR) correction values from Greg Salvesen's repository.
    The function retrieves and uses data from the file 'gGR_gNT_J1655.h5' for these calculations.

    Args:
        inc (array-like): Inclination angles in degrees. This should contain at least the minimum and maximum of the desired inclination range.
        a (float, optional): Spin parameter, a dimensionless quantity representing the angular momentum of the black hole. Defaults to 0.
        verbose (bool, optional): If True, enables printing of warnings and additional information. Defaults to True.
        inc_sample_method (str, optional): Method for sampling inclination values, can be 'uniform' or 'isotropic'. Defaults to 'uniform'.
        i_indicies (list of int, optional): Pre-calculated inclination indices. If None, indices are calculated based on the inclination sample method.

    Returns:
        tuple: A tuple containing the following elements:
            - The GR correction factors multiplied by the cosine of inclination and limb darkening.
            - Rin/Rg ratio value.
            - Selected spin value from the spin grid.
            - Array of selected inclination angles.
            - Indices of the selected inclination angles in the inclination grid.

    References:
        Greg Salvesen's repository for GR correction values: https://github.com/gregsalvesen/bhspinf
    """
    if not os.path.isfile('gGR_gNT_J1655.h5'):  
        urllib.request.urlretrieve("https://raw.githubusercontent.com/gregsalvesen/bhspinf/main/data/GR/gGR_gNT_J1655.h5","gGR_gNT_J1655.h5") ## Thanks Greg!!
    ## Copied from gregsalvesen/bhspinf/
    fh5      = 'gGR_gNT_J1655.h5'
    f        = h5py.File(fh5, 'r')
    a_grid   = f.get('a_grid')[:]  # [-]
    r_grid   = f.get('r_grid')[:]  # [Rg]
    i_grid   = f.get('i_grid')[:]  # [deg]
    gGR_grid = f.get('gGR')[:,:]   # gGR(r[Rg], i[deg])
    gNT_grid = f.get('gNT')[:]     # gNT(r[Rg])
    f.close()
    ##

    atol = 1e-08
    a_index, = np.where(np.isclose(a_grid,a,rtol=1e-05, atol=atol))
    atol_tmp = atol
    while len(a_index) == 0:
        atol_tmp *= 10
        a_index, = np.where(np.isclose(a_grid,a,rtol=1e-05, atol=atol_tmp))
    if len(a_index) != 1:
        a_index = a_index + (a_index[-1]-a_index[0]) // 2
        if verbose==True:
            print("Warning: found more than one spin grid value (for GR correction) that are equally close to the user specified spin value. Taking the median value %s as the one to be closest." % (a_grid[a_index[0]]))
    Rin_ratio = r_grid[a_index]

    inc_index = idx_of_value_from_grid(i_grid,inc,verbose=True)
    
    GR_correction = gGR_grid[a_index[0],inc_index[0]] * gNT_grid[a_index[0]]
    limb_darkening = (1/2 + (3/4)*np.cos(i_grid[inc_index[0]]*(np.pi/180)))
    cos_i = np.cos(i_grid[inc_index[0]]*(np.pi/180))
    selected_i = i_grid[inc_index[0]]

    # GR_correction = 1 ##Uncomment for removing GR corrections
    if not limb_dark:
      limb_darkening = 1

    return (GR_correction*cos_i*limb_darkening,Rin_ratio[0],a_grid[a_index[0]],selected_i)

def to_norm(d,mass,a,inc,limb_dark=True):

    G = 6.6743e-11  # SI units
    c = 2.998e8  # SI units
    kappa = 1.7
    m_sun = 1.989e30  # SI units

    # GR correction values
    corr, Rin_ratio, _, _= get_total_correction_GR_and_Rin_Rg_ratio(inc, a, limb_dark=limb_dark, verbose=True)

    norm = (((Rin_ratio * G * m_sun * 1e-2) / ((kappa ** 2) * (c ** 2)))**2) * corr * (mass/d)**2

    return norm

def to_d(norm,mass,a,inc,limb_dark=True):

    G = 6.6743e-11  # SI units
    c = 2.998e8  # SI units
    kappa = 1.7
    m_sun = 1.989e30  # SI units

    # GR correction values
    corr, Rin_ratio, _, _= get_total_correction_GR_and_Rin_Rg_ratio(inc, a, limb_dark=limb_dark, verbose=True)

    d = ((Rin_ratio * G * m_sun * 1e-2) / ((kappa ** 2) * (c ** 2))) * np.sqrt(1 / norm) * np.sqrt(corr) * (mass)

    return d

def scale_powerlaw_norm(gamma,temp,ezdiskbb_norm,ratio_pl_to_disk):
   
   AllModels.clear()
   AllData.clear()

   model = "(po+ezdiskbb)"
   m = Model(model,setPars={1:gamma,3:temp,4:ezdiskbb_norm})
   
   m.setPars({2: 0})
   AllModels.calcFlux("2.0 20.0")
   ezdiskbb_flux = m.flux[0]
   m.setPars({2: 1},{4: 0})
   AllModels.calcFlux("2.0 20.0")
   powerlaw_flux = m.flux[0]
   pl_norm = (ezdiskbb_flux * ratio_pl_to_disk) / powerlaw_flux

   return pl_norm

def run_simulation(arguments):
    Xset.seed = random.randint(0, 10000)
    nH_value, d, args, iteration= arguments

    ezdiskbb_norm = to_norm(d,args.mass,args.a,args.inc,limb_dark=True)

    powerlaw_norm = scale_powerlaw_norm(args.gamma,args.temp,ezdiskbb_norm,(1-args.ratio_disk_to_tot)/args.ratio_disk_to_tot)

    gamma_fit_range = "2.3,,1.7,1.7,3.0,3.0"

    result = {"nH": nH_value, "d": d, "red_chi_squared": None, "gamma": None, "power_norm_fake": powerlaw_norm, "power_norm_fit": None, "temp": None, "disk_norm_fake": ezdiskbb_norm, "disk_norm_fit": None, "error_disk_norm_low": None, "error_disk_norm_up": None, "d_fit": None, "error_d_low": None , "error_d_up": None, "frac_uncert": None}

    sim1 = simulation("tbabs*(po+ezdiskbb)",'maxi',{1: nH_value, 2:args.gamma, 3: powerlaw_norm, 4: args.temp, 5: ezdiskbb_norm},{1: str(nH_value) + ",0", 2: gamma_fit_range, 4: ',,0.1,0.1'})
    m = sim1.run(id=iteration,exposure=args.exposure,backExposure=args.exposure)

    try:
        result.update({"red_chi_squared": Fit.statistic / Fit.dof, "gamma": m.powerlaw.PhoIndex.values[0], "power_norm_fit": m.powerlaw.norm.values[0], "temp": m.ezdiskbb.T_max.values[0], "disk_norm_fit": m.ezdiskbb.norm.values[0], "error_disk_norm_low": m.ezdiskbb.norm.error[0], "error_disk_norm_up": m.ezdiskbb.norm.error[1], "d_fit": to_d(m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True),"error_d_low": to_d(m.ezdiskbb.norm.error[1],args.mass,args.a,args.inc,limb_dark=True),"error_d_up": to_d(m.ezdiskbb.norm.error[0],args.mass,args.a,args.inc,limb_dark=True), "frac_uncert": ((to_d(m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True)- to_d(m.ezdiskbb.norm.error[1],args.mass,args.a,args.inc,limb_dark=True)) + (to_d(m.ezdiskbb.norm.error[0],args.mass,args.a,args.inc,limb_dark=True)- to_d(m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True)) / 2) / (to_d(m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True))})
    except:
            pass

    return result

if __name__ == "__main__":

    Xset.parallel.leven = 2
    Xset.parallel.error = 2
    Xset.chatter = 0
    Xset.logChatter = 0

    
    parser = argparse.ArgumentParser(description='Tests nH effect for different distance on fake spectra given a gamma, temp, spin (a), mass, inclination (inc), ratio_disk_to_pl, and spectrum exposure')
    parser.add_argument('gamma', type=float)
    parser.add_argument('temp', type=float)
    parser.add_argument('a', type=float)
    parser.add_argument('mass', type=float)
    parser.add_argument('inc', type=float)
    parser.add_argument('ratio_disk_to_tot', type=float)
    parser.add_argument('exposure', type=float)

    # Parse the argument
    args = parser.parse_args()

    d_list = [1,2,3,4,5,6,8,12,18,26]
    nH_list = [0.1]

    # d_list = [2]
    # nH_list = [1.0]

    start_time = time.perf_counter()

    all_args = []
    for nH_value in nH_list:
        for d in d_list: 
            for iteration in range(300):
                all_args.append((nH_value,d,args,iteration))

    with Pool(int(cpu_count()/2) - 1) as pool:  # Use all but one CPU core   
        results = []
        try:
            it = pool.imap(run_simulation, all_args, chunksize=1)
            for _ in tqdm(range(len(all_args)), desc="Running simulations", position=0, leave=True):
                try:
                    result = it.next(timeout=30)  # Timeout set to 30 sec per task
                    results.append(result)
                except TimeoutError:
                    print("A task has timed out. Skipping this task.")
                    continue  # Skip the hanging task and move on to the next one
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            pool.terminate()
            pool.join()
    
    
    df_full = pd.DataFrame(results)
    df_full.to_csv("/disk/data/youssef/scripts/xrb-population/results/table_g"+str(args.gamma)+"_T"+str(args.temp)+"_a"+str(args.a)+"_m"+str(args.mass)+"_i"+str(args.inc)+"_r"+str(args.ratio_disk_to_tot)+"_e"+str(args.exposure)+"_full_1nH.csv", index=False)

    table_red = []

    for nH_value in nH_list:
        for d in d_list: 
            filtered_results = [res for res in results if res["d"] == d and res["nH"] == nH_value]
            df = pd.DataFrame(filtered_results)
            # table_red.append({"nH": nH_value, "red_chi_squared": df["red_chi_squared"].median(), "gamma": df["gamma"].median(), "power_norm_fake": filtered_results[0]["power_norm_fake"], "power_norm_fit": df["power_norm_fit"].median(), "temp": df["temp"].median(), "disk_norm_fake": filtered_results[0]["disk_norm_fake"], "disk_norm_fit": df["disk_norm_fit"].median(), "error_disk_norm": df["disk_norm_fit"].median() - filtered_results[0]["disk_norm_fake"],"d": d,"d_fit": df["d_fit"].median(),"error_d": df["d_fit"].median() - filtered_results[0]["d"], "frac_uncert": (df["d_fit"].median() - filtered_results[0]["d"]) / filtered_results[0]["d"], "med_frac_uncert": df["frac_uncert"].median()})
            table_red.append({
                "nH": nH_value,
                "red_chi_squared": df["red_chi_squared"].median() if 'red_chi_squared' in df.columns else None,
                "gamma": df["gamma"].median() if 'gamma' in df.columns else None,
                "power_norm_fake": filtered_results[0]["power_norm_fake"] if filtered_results else None,
                "power_norm_fit": df["power_norm_fit"].median() if 'power_norm_fit' in df.columns else None,
                "temp": df["temp"].median() if 'temp' in df.columns else None,
                "disk_norm_fake": filtered_results[0]["disk_norm_fake"] if filtered_results else None,
                "disk_norm_fit": df["disk_norm_fit"].median() if 'disk_norm_fit' in df.columns else None,
                "error_disk_norm": (df["disk_norm_fit"].median() - filtered_results[0]["disk_norm_fake"]) if 'disk_norm_fit' in df.columns and filtered_results else None,
                "d": d,
                "d_fit": df["d_fit"].median() if 'd_fit' in df.columns else None,
                "error_d": (df["d_fit"].median() - filtered_results[0]["d"]) if 'd_fit' in df.columns and filtered_results else None,
                "frac_uncert": ((df["d_fit"].median() - filtered_results[0]["d"]) / filtered_results[0]["d"]) if 'd_fit' in df.columns and filtered_results else None,
                "med_frac_uncert": df["frac_uncert"].median() if 'frac_uncert' in df.columns else None
            })

    df_red = pd.DataFrame(table_red)
    df_red.to_csv("/disk/data/youssef/scripts/xrb-population/results/table_g"+str(args.gamma)+"_T"+str(args.temp)+"_a"+str(args.a)+"_m"+str(args.mass)+"_i"+str(args.inc)+"_r"+str(args.ratio_disk_to_tot)+"_e"+str(args.exposure)+"_1nH.csv", index=False)

    command = f'rm -rf gGR_gNT_J1655.h5'
    process = subprocess.Popen(command, shell=True)
    process.wait()

    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(f"The script took {total_time} seconds to complete.")