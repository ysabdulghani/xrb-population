from xspec import *
import argparse
from data_read import findSavedModel
import os
import subprocess
import pandas as pd

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

m.setPars({1:nH_origin,3:powerlaw_origin_norm},{5:ezdiskbb_origin_norm})

print(f"flux ratio is {ratio_flux}")

d_origin = 8.13
d_list = [d_origin,2,4,6,8,10,12,14,16,18,20]
nH_list = [nH_origin,1.0, 5.0, 10.0]
table = []

for d_ratio in [d_origin/x for x in d_list]:
    for nH_value in nH_list:  # Example nH values

        AllModels.clear()
        AllData.clear()

        os.chdir(args.sourceSpectrumDir)
        XspecSettings.restore(saveModelFile,saveModelFile)
        os.chdir(initDir)

        ezdiskbb_new_norm = ezdiskbb_origin_norm*(d_ratio**2)
        m.setPars({1:0},{3:0},{5:ezdiskbb_new_norm})
        AllModels.calcFlux("2.0 20.0")
        s1 = AllData(1)
        ezdiskbb_new_flux = s1.flux[0]
        powerlaw_new_norm =  powerlaw_origin_norm * ((ezdiskbb_new_flux*ratio_flux) / powerlaw_flux)
        m.setPars({1:nH_value},{3:powerlaw_new_norm},{5:ezdiskbb_new_norm})


        fake_settings = FakeitSettings(fileName = "fakeit_tmp.pha")
        AllData.fakeit(1,fake_settings, applyStats=True)

        command = f'ftgrouppha fakeit_tmp.pha backfile=fakeit_tmp_bkg.pha fakeit_tmp_binned.pha snmin 3'         
        # Execute the command
        process = subprocess.Popen(command, shell=True)
        process.wait()
        print(f'Executed command: {command}')  

        AllData.clear()
        gamma = "2.3,,1.7,1.7,3.0,3.0"
        s1 = Spectrum("fakeit_tmp_binned.pha")
        AllData.ignore("bad")
        s1.ignore("**-2.0,20.0-**")
        AllModels.clear()
        Model("tbabs*(po+ezdiskbb)",setPars={1:str(nH_value)+",0",2:gamma,4:',,0.1,0.1'})
        Fit.query = "yes"
        Fit.perform()
        
        try:
            Fit.error("1-5")
            Plot.xAxis = "KeV"
            Plot("ldata delchi")
            Plot()
            table.append({"nH: ":nH_value, "d": d_origin/d_ratio, "red_chi_squared": Fit.statistic/Fit.dof, "disk_norm": m.ezdiskbb.norm.values[0],"error_disk_norm_low": m.ezdiskbb.norm.error[0],"error_disk_norm_up": m.ezdiskbb.norm.error[1],"frac_uncert": ((m.ezdiskbb.norm.error[0]+m.ezdiskbb.norm.error[1])/2)/(m.ezdiskbb.norm.values[0])})                 
        except:
            table.append({"nH: ":nH_value, "d": d_origin/d_ratio, "red_chi_squared": Fit.statistic/Fit.dof, "disk_norm": m.ezdiskbb.norm.values[0],"error_disk_norm_low": None,"error_disk_norm_up": None,"frac_uncert": None})


df = pd.DataFrame(table)

# Write the DataFrame to a CSV file
df.to_csv("/disk/data/youssef/scripts/xrb-population/results/table_"+args.sourceSpectrumDir+".csv", index=False)
