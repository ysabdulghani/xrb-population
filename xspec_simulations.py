
from xspec import *
import subprocess
import heasoftpy as hsp

class simulation:
    
    def __init__(self,model_def,instrument,simulation_params_dic,fit_params_dic):
        '''
        Initiates the simulation object

        Arguments:
        model_def: defines the xspec model used in simulation and fitting (can easily edit to account for different models)
        instrument: will apply energy range, response file and background file specific to the instrument
        simulation_params_dic: the dictionary (in PyXspec notation) that specifics the parameters of the faked spectrum (example: {1:'0.5,0',5:',,0.1,0.1'})
        fit_params_dic : similar to simulation_params_dic but for the model to be fitted

        '''
        self.model =  model_def
        self.sim_params_dic = simulation_params_dic 
        self.fit_params_dic = fit_params_dic
        if instrument == 'maxi': # Will need to change this part according to your need
            self.energyRange_low = '2.0'
            self.energyRange_high= '20.0'
            self.responseFilename =  "sim_files/gx339-4_g_low.rsp" # Need to be short name and in same directory where your run simulation
            self.backgroundFilename = "sim_files/gx339-4_g_low_bgd.pi"
        elif instrument == 'xrt':
            self.energyRange_low = '0.7'
            self.energyRange_high= '10.0'
            self.responseFilename =  "sim_files/swxwt0to2s6_20131212v015.rmf" # Need to be short name and in same directory where your run simulation
            self.backgroundFilename = "sim_files/00010627114bgd_wt.pha"
        else:
            raise ValueError('Only maxi or xrt allowed.')

    def run(self,id='',**kwargs):

        '''
        Perform a simulation run (fake a spectrum and then fit it)

        Arguments:
        id: needed for the temp fake it file (especially when running multiprocessing simulation)
        **kwargs: This is passed to the FakeitSettings object from PyXspec. Needed to change exposure of the faked spectrum for example

        Output:
        A fitted model object. If fitting fails it will return the unfitted model object 
 
        '''
        
        AllModels.clear()
        

        Model(self.model,setPars=self.sim_params_dic)

        fake_settings = FakeitSettings(response=self.responseFilename,background=self.backgroundFilename,fileName="fakeit_tmp_"+str(id)+".pha",**kwargs)
        AllData.fakeit(1, fake_settings, applyStats=True)

        with hsp.utils.local_pfiles_context():
            hsp.ftgrouppha(infile='fakeit_tmp_'+str(id)+'.pha',backfile='fakeit_tmp_'+str(id)+'_bkg.pha',outfile='fakeit_tmp_'+str(id)+'_binned.pha',grouptype='snmin',groupscale='3')

        # command = f'ftgrouppha fakeit_tmp_'+str(id)+'.pha backfile=fakeit_tmp_'+str(id)+'_bkg.pha fakeit_tmp_'+str(id)+'_binned.pha snmin 3'
        # process = subprocess.Popen(command, shell=True)
        # process.wait()
        AllData.clear()
        s1 = Spectrum("fakeit_tmp_"+str(id)+"_binned.pha")
        AllData.ignore("bad")
        s1.ignore("**-"+self.energyRange_low+","+self.energyRange_high+"-**")
        AllModels.clear()
        fitModel = Model(self.model,setPars=self.fit_params_dic)
        Fit.query = "yes"

        try:
            Fit.perform()
            try:
                Fit.error("1-"+str(fitModel.nParameters))
            except:
                pass
        except:
            pass
        
        command = f'rm -rf fakeit_tmp_'+str(id)+'.pha fakeit_tmp_'+str(id)+'_bkg.pha fakeit_tmp_'+str(id)+'_binned.pha snmin 3'
        process = subprocess.Popen(command, shell=True)
        process.wait()

        return fitModel
    




















