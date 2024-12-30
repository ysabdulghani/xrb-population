from xspec import *
import argparse
import os
import pandas as pd
import numpy as np
import time
from multiprocessing import Pool, cpu_count, TimeoutError, set_start_method
from tqdm import tqdm
from xspec_simulations import *
import random
import urllib
import h5py
import subprocess
import tempfile
import sys
import gc
from datetime import datetime

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


def correction_file():
    """
    The function retrieves and uses data from the file 'gGR_gNT_J1655.h5' for these calculations.
    References:
        Greg Salvesen's repository for GR correction values: https://github.com/gregsalvesen/bhspinf
    """
    if not os.path.isfile('gGR_gNT_J1655.h5'):  
        urllib.request.urlretrieve("https://raw.githubusercontent.com/gregsalvesen/bhspinf/main/data/GR/gGR_gNT_J1655.h5","gGR_gNT_J1655.h5") ## Thanks Greg!!
    ## Copied from gregsalvesen/bhspinf/
    fh5      = 'gGR_gNT_J1655.h5'
    with h5py.File(fh5, 'r') as f:
        a_grid   = f['a_grid'][:]    # spin values
        r_grid   = f['r_grid'][:]    # radius in Rg
        i_grid   = f['i_grid'][:]    # inclination angles
        gGR_grid = f['gGR'][:,:]     # 2D correction factor
        gNT_grid = f['gNT'][:]       # 1D correction factor

        # Return a dictionary of arrays (or any other structure you prefer)
        data = {
        'a_grid': a_grid,
        'r_grid': r_grid,
        'i_grid': i_grid,
        'gGR_grid': gGR_grid,
        'gNT_grid': gNT_grid
        }
    
    return data


def get_total_correction_GR_and_Rin_Rg_ratio(data,inc,a,limb_dark=True,verbose=True):
    """
    Calculates the total correction factor and the Rin/Rg ratio based on the provided inclination angle and spin parameter.
    The correction factors are derived using General Relativity (GR) correction values from Greg Salvesen's repository.

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
    a_grid   = data['a_grid']  # [-]
    r_grid   = data['r_grid']  # [Rg]
    i_grid   = data['i_grid']  # [deg]
    gGR_grid = data['gGR_grid']  # gGR(r[Rg], i[deg])
    gNT_grid = data['gNT_grid']     # gNT(r[Rg])
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

def to_norm(f,d,mass,a,inc,limb_dark=True):

    G = 6.6743e-11  # SI units
    c = 2.998e8  # SI units
    kappa = 1.7
    m_sun = 1.989e30  # SI units

    # GR correction values
    corr, Rin_ratio, _, _= get_total_correction_GR_and_Rin_Rg_ratio(f,inc, a, limb_dark=limb_dark, verbose=True)

    norm = (((Rin_ratio * G * m_sun * 1e-2) / ((kappa ** 2) * (c ** 2)))**2) * corr * (mass/d)**2

    return norm

def to_d(f,norm,mass,a,inc,limb_dark=True):

    G = 6.6743e-11  # SI units
    c = 2.998e8  # SI units
    kappa = 1.7
    m_sun = 1.989e30  # SI units

    # GR correction values
    corr, Rin_ratio, _, _= get_total_correction_GR_and_Rin_Rg_ratio(f,inc, a, limb_dark=limb_dark, verbose=True)

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

   AllModels.clear()
   AllData.clear()

   return pl_norm

def run_simulation(arguments):
    Xset.seed = random.randint(0, 10000)
    nH_value, d, args, iteration, tmp_dir, f = arguments

    ezdiskbb_norm = to_norm(f,d,args.mass,args.a,args.inc,limb_dark=True)

    powerlaw_norm = scale_powerlaw_norm(args.gamma,args.temp,ezdiskbb_norm,(1-args.ratio_disk_to_tot)/args.ratio_disk_to_tot)

    gamma_fit_range = "2.3,,1.7,1.7,3.0,3.0"

    result = {"nH": nH_value, "d": d, "red_chi_squared": None, "gamma": None, "power_norm_fake": powerlaw_norm, "power_norm_fit": None, "temp": None, "disk_norm_fake": ezdiskbb_norm, "disk_norm_fit": None, "error_disk_norm_low": None, "error_disk_norm_up": None, "d_fit": None, "error_d_low": None , "error_d_up": None, "frac_uncert": None}

    AllModels.clear()
    AllData.clear()
    
    sim1 = simulation("tbabs*(po+ezdiskbb)",args.instrument,{1: nH_value, 2:args.gamma, 3: powerlaw_norm, 4: args.temp, 5: ezdiskbb_norm},{1: str(nH_value) + ",0", 2: gamma_fit_range, 4: ',,0.1,0.1'})
    m = sim1.run(id=iteration,spec_dir=tmp_dir,exposure=args.exposure,backExposure=args.exposure)

    try:
        result.update({"red_chi_squared": Fit.statistic / Fit.dof, "gamma": m.powerlaw.PhoIndex.values[0], "power_norm_fit": m.powerlaw.norm.values[0], "temp": m.ezdiskbb.T_max.values[0], "disk_norm_fit": m.ezdiskbb.norm.values[0], "error_disk_norm_low": m.ezdiskbb.norm.error[0], "error_disk_norm_up": m.ezdiskbb.norm.error[1], "d_fit": to_d(f,m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True),"error_d_low": to_d(f,m.ezdiskbb.norm.error[1],args.mass,args.a,args.inc,limb_dark=True),"error_d_up": to_d(f,m.ezdiskbb.norm.error[0],args.mass,args.a,args.inc,limb_dark=True), "frac_uncert": ((to_d(f,m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True)- to_d(f,m.ezdiskbb.norm.error[1],args.mass,args.a,args.inc,limb_dark=True)) + (to_d(f,m.ezdiskbb.norm.error[0],args.mass,args.a,args.inc,limb_dark=True)- to_d(f,m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True)) / 2) / (to_d(f,m.ezdiskbb.norm.values[0],args.mass,args.a,args.inc,limb_dark=True))})
    except:
            pass

    return result

# Helper function to create a new pool and imap iterator
def new_pool_and_iterator(processes,start_index):
    pool_ = Pool(processes=processes)
    # Only map the remaining tasks (from start_index onward)
    it_ = pool_.imap(run_simulation, all_args[start_index:], chunksize=1)
    return pool_, it_

def main(all_args):
    max_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    processes = max_cores - 2 # e.g., up to 100, or just use max_cores
    # processes = 200
    print(f"Starting simulations with up to {processes} processes")

    results = []
    timed_out_tasks = []
    errored_tasks = []

    timeout_counter = 0
    max_consecutive_timeouts = 10

    idx = 0  # We'll manually track the index over `all_args`.



    # Initialize the first pool & iterator
    pool, it = new_pool_and_iterator(processes,idx)

    with tqdm(total=len(all_args), desc="Running simulations") as pbar:

        while idx < len(all_args):

            if timeout_counter > max_consecutive_timeouts:
                err_msg = "Maximum consecutive timeouts exceeded. Stopping."
                print(err_msg)
                script_call = " ".join(sys.argv)
                current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                with open("error_log.log", "a") as error_file:
                        error_file.write(f"{current_time}: {script_call}: {err_msg} \n")
                break

            try:
                # Attempt to get the next result with a timeout
                result = it.next(timeout=30)
                results.append(result)
                # Successfully got a result => reset timeout counter
                timeout_counter = 0
                pbar.update(1)
                idx += 1

            except TimeoutError:
                timed_out_tasks.append(all_args[idx][:-1])
                timeout_counter += 1
                pbar.update(1)
                print(f"Timeout #{timeout_counter} at index={idx}. Terminating pool and restarting...")
                idx += 1  # Skip the timed-out task
                
                pool.terminate()  
                pool.join()
                print("Pool is terminated")

                try:
                    print("Forcing grabage collection before attempting to restart pool")
                    gc.collect()
                    # Re-create the pool and iterator for remaining tasks
                    print("Attempting to restart pool...")
                    pool, it = new_pool_and_iterator(processes,idx)
                except OSError as e:
                    err_msg = f"Failed to restart pool due to OSError: {e}"
                    print(err_msg)
                    script_call = " ".join(sys.argv)
                    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                    with open("error_log.log", "a") as error_file:
                            error_file.write(f"{current_time}: {script_call}: {err_msg} \n")

            except Exception as e:
                errored_tasks.append((all_args[idx], str(e)))
                print(f"Task at index={idx} errored out: {e}. Skipping.")
                timeout_counter = 0  # reset for non-timeout errors
                pbar.update(1)
                idx += 1
                # We do NOT terminate the pool here, but you could if desired

    # Cleanly close any remaining pool after finishing
    if timeout_counter > max_consecutive_timeouts:
        pool.terminate()
    else:
        pool.close()
    pool.join()

    print("Loop finished. Returning results.")
    return results, timed_out_tasks, errored_tasks


if __name__ == "__main__":

    # set_start_method('spawn')
    # random.seed(42)

    Xset.parallel.leven = 1
    Xset.parallel.error = 1
    Xset.chatter = 0
    Xset.logChatter = 0
    Xset.allowPrompting = False

    
    parser = argparse.ArgumentParser(description='Tests nH effect for different distance on fake spectra given a gamma, temp, spin (a), mass, inclination (inc), ratio_disk_to_pl, and spectrum exposure')
    parser.add_argument('gamma', type=float)
    parser.add_argument('temp', type=float)
    parser.add_argument('a', type=float)
    parser.add_argument('mass', type=float)
    parser.add_argument('inc', type=float)
    parser.add_argument('ratio_disk_to_tot', type=float)
    parser.add_argument('exposure', type=float)
    parser.add_argument('instrument', type=str)

    # Parse the argument
    args = parser.parse_args()

    d_list = [1,2,3,4,5,6,8,12,18,26]
    nH_list = [0.1,0.5,5,10]

    # d_list = [2]
    # nH_list = [1.0]

    start_time = time.perf_counter()

    all_args = []
    
    tmp_dir_name = tempfile.mkdtemp(prefix="tmp_", dir="/dev/shm")
    print(f"Created temporary directory: {tmp_dir_name}")

    counter = 0  # Initialize a counter

    f = correction_file()


    all_args = []
    for nH_value in nH_list:
        for d in d_list: 
            for iteration in range(300):
                unique_iteration = counter  # Use the counter as a unique identifier
                all_args.append((nH_value, d, args, unique_iteration, tmp_dir_name,f))
                counter += 1

    results, timeouts, errors = main(all_args)

    os.makedirs("results/"+str(args.instrument)+"_results", exist_ok=True)
    
    df_full = pd.DataFrame(results)
    df_full.to_csv("results/"+str(args.instrument)+"_results/table_g"+str(args.gamma)+"_T"+str(args.temp)+"_a"+str(args.a)+"_m"+str(args.mass)+"_i"+str(args.inc)+"_r"+str(args.ratio_disk_to_tot)+"_e"+str(args.exposure)+"_full.csv", index=False)
    
    table_red = []

    for nH_value in nH_list:
        for d in d_list: 
            filtered_results = [res for res in results if res["d"] == d and res["nH"] == nH_value]
            df = pd.DataFrame(filtered_results)
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
    df_red.to_csv("results/"+str(args.instrument)+"_results/table_g"+str(args.gamma)+"_T"+str(args.temp)+"_a"+str(args.a)+"_m"+str(args.mass)+"_i"+str(args.inc)+"_r"+str(args.ratio_disk_to_tot)+"_e"+str(args.exposure)+".csv", index=False)

    command = f'rm -rf '+tmp_dir_name
    process = subprocess.Popen(command, shell=True)
    process.wait()

    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(f"The script took {total_time} seconds to complete.")
    print("Timed Out:", timeouts)
    print("Errored:", errors)