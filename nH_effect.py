from xspec import *
import argparse
from data_read import findSavedModel
import os
import subprocess

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

d = 1

print(f"flux ratio is {ratio_flux}")


ezdiskbb_new_norm = ezdiskbb_origin_norm/(d**2)
m.setPars({3:0},{5:ezdiskbb_new_norm})
AllModels.calcFlux("2.0 20.0")
s1 = AllData(1)
ezdiskbb_new_flux = s1.flux[0]
powerlaw_new_norm =  powerlaw_origin_norm * ((ezdiskbb_new_flux*ratio_flux) / powerlaw_flux)
print(f"norm ratio calculated from flux is {((ezdiskbb_new_flux*ratio_flux) / powerlaw_flux)}")
m.setPars({1:10},{3:powerlaw_new_norm},{5:ezdiskbb_new_norm})

m.show()

# AllData.fakeit(1, applyStats=True, noWrite=True)

# for nH_value in [1.0, 5.0, 10.0]:  # Example nH values
nH_value = 10.0
m.TBabs.nH = nH_value
gamma = "2.3,,1.7,1.7,3.0,3.0"
fake_settings = FakeitSettings(fileName = "fakeit_tmp.pha")
AllData.fakeit(1,fake_settings, applyStats=True)

command = f'ftgrouppha fakeit_tmp.pha backfile=fakeit_tmp_bkg.pha fakeit_tmp_binned.pha snmin 3'         
# Execute the command
process = subprocess.Popen(command, shell=True)
process.wait()
print(f'Executed command: {command}')  


AllData.clear()
s1 = Spectrum("fakeit_tmp_binned.pha")
AllData.ignore("bad")
s1.ignore("**-2.0,20.0-**")
AllModels.clear()
Model("tbabs*(po+ezdiskbb)",setPars={1:str(nH_value)+",0",2:gamma,4:',,0.1,0.1'})
Fit.query = "yes"
Fit.perform()
Fit.error("1-5")
Plot.xAxis = "KeV"
Plot("ldata delchi")
# xs.Plot.setRebin(minSig=25, maxBins=10, groupNum=-1, errType="quad")
Plot()
print(f"nH: {nH_value}, red_chi_squared: {Fit.statistic/Fit.dof}, Disk Norm: {m.ezdiskbb.norm.values[0]}, Error (90%) in Disk Norm: {m.ezdiskbb.norm.error[0:2]}")

