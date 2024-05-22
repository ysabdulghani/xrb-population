from xspec import *
import argparse
from data_read import findSavedModel
import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import time
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from xspec_simulations import *
import random


# def spin_to_rin():

#     return spin

# def to_norm(d,spin,inc,limb_dark=True):

#     return norm

# def scale_powerlaw_norm(ezdiskbb_flux,ratio_disk_to_pl):

#    return pl_norm

def run_simulation(arguments):
    Xset.seed = random.randint(0, 100) # Need to fix or make sure you generate from a larger pool of int
    d_ratio, nH_value, iteration, args = arguments

    os.chdir(args.sourceSpectrumDir)
    XspecSettings.restore(saveModelFile, saveModelFile)
    os.chdir(initDir)
    m = AllModels(1)

    ezdiskbb_new_norm = ezdiskbb_origin_norm * (d_ratio**2)
    m.setPars({1: 0}, {3: 0}, {5: ezdiskbb_new_norm})
    AllModels.calcFlux("2.0 20.0")
    s1 = AllData(1)
    ezdiskbb_new_flux = s1.flux[0]
    try:
        powerlaw_new_norm = powerlaw_origin_norm * ((ezdiskbb_new_flux * ratio_flux) / powerlaw_flux)
    except:
        powerlaw_new_norm = 0

    gamma = "2.3,,1.7,1.7,3.0,3.0"

    sim1 = simulation("tbabs*(po+ezdiskbb)",'maxi',{1: nH_value, 3: powerlaw_new_norm, 5: ezdiskbb_new_norm},{1: str(nH_value) + ",0", 2: gamma, 4: ',,0.1,0.1'})
    m = sim1.run(id=iteration,exposure=1500,backExposure=1500)

    result = {"nH": nH_value, "d": d_origin / d_ratio, "red_chi_squared": None, "gamma": None, "power_norm_fake": powerlaw_new_norm, "power_norm_fit": None, "temp": None, "disk_norm_fake": ezdiskbb_new_norm, "disk_norm_fit": None, "error_disk_norm_low": None, "error_disk_norm_up": None, "frac_uncert": None}

    try:
        result.update({"red_chi_squared": Fit.statistic / Fit.dof, "gamma": m.powerlaw.PhoIndex.values[0], "power_norm_fit": m.powerlaw.norm.values[0], "temp": m.ezdiskbb.T_max.values[0], "disk_norm_fit": m.ezdiskbb.norm.values[0], "error_disk_norm_low": m.ezdiskbb.norm.error[0], "error_disk_norm_up": m.ezdiskbb.norm.error[1], "frac_uncert": (((m.ezdiskbb.norm.values[0] - m.ezdiskbb.norm.error[0]) + (m.ezdiskbb.norm.error[1] - m.ezdiskbb.norm.values[0])) / 2) / (m.ezdiskbb.norm.values[0])})
    except:
            pass

    return result

if __name__ == "__main__":

    Xset.parallel.leven = 2
    Xset.parallel.error = 2
    Xset.chatter = 0
    Xset.logChatter = 0

    parser = argparse.ArgumentParser(description='Tests nH effect on fake spectra')
    parser.add_argument('sourceSpectrumDir', type=str)

    # Parse the argument
    args = parser.parse_args()

    initDir = os.getcwd()
    os.chdir(args.sourceSpectrumDir)
    print(f"Using source spectrum in: {os.getcwd()}")

    saveModelFile = findSavedModel()
    if saveModelFile:
        XspecSettings.restore(saveModelFile,saveModelFile)

    os.chdir(initDir)
    print(os.getcwd())

    m = AllModels(1)
    ezdiskbb_origin_norm = m.ezdiskbb.norm.values[0]
    powerlaw_origin_norm = m.powerlaw.norm.values[0]
    nH_origin = m.TBabs.nH.values[0]

    m.show()
    m.setPars({1:0},{5:0})
    AllModels.calcFlux("2.0 20.0")
    s1 = AllData(1)
    powerlaw_flux = s1.flux[0]
    m.setPars({3:0},{5:ezdiskbb_origin_norm})
    AllModels.calcFlux("2.0 20.0")
    s1 = AllData(1)
    ezdiskbb_flux = s1.flux[0]
    ratio_flux = powerlaw_flux/ezdiskbb_flux


    d_origin = 8.13
    nH_origin = 0.5  # Define nH_origin or assign it a meaningful value

    d_list = [x for x in range(1, 25, 1)]
    nH_list = [x * 0.5 for x in range(2, 21)]

    d_list.append(d_origin)
    nH_list.append(nH_origin)

    # d_list = [2,8]
    # nH_list = [1.0,10]

    table_full = []
    table_red = []

    start_time = time.perf_counter()

    all_args = []
    for d_ratio in [d_origin / x for x in d_list]:
        for nH_value in nH_list:
            for iteration in range(300):
                all_args.append((d_ratio, nH_value, iteration, args))

    with Pool(int(cpu_count()/2) - 1) as pool:  # Use all but one CPU core
        results = list(tqdm(pool.imap(run_simulation, all_args,chunksize=1), total=len(all_args), desc="Running simulations",position=0, leave=True))

    for result in results:
        table_full.append(result)

    df_full = pd.DataFrame(table_full)
    df_full.to_csv("/disk/data/youssef/scripts/xrb-population/results/table_full_test_latest.csv", index=False)

    for d_ratio in [d_origin / x for x in d_list]:
        for nH_value in nH_list:
            filtered_results = [res for res in results if res["d"] == d_origin / d_ratio and res["nH"] == nH_value]
            df = pd.DataFrame(filtered_results)
            table_red.append({"nH": nH_value, "d": d_origin / d_ratio, "red_chi_squared": df["red_chi_squared"].median(), "gamma": df["gamma"].median(), "power_norm_fake": filtered_results[0]["power_norm_fake"], "power_norm_fit": df["power_norm_fit"].median(), "temp": df["temp"].median(), "disk_norm_fake": filtered_results[0]["disk_norm_fake"], "disk_norm_fit": df["disk_norm_fit"].median(), "abs_l1_error_disk_norm": df["disk_norm_fit"].median() - filtered_results[0]["disk_norm_fake"], "frac_uncert": (df["disk_norm_fit"].median() - filtered_results[0]["disk_norm_fake"]) / filtered_results[0]["disk_norm_fake"], "med_frac_uncert": df["frac_uncert"].median()})

    df_red = pd.DataFrame(table_red)
    df_red.to_csv("/disk/data/youssef/scripts/xrb-population/results/table_test_latest.csv", index=False)

    end_time = time.perf_counter()
    total_time = end_time - start_time
    print(f"The loop took {total_time} seconds to complete.")