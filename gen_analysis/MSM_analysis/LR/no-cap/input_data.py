__author__ = 'JFR'

import numpy as np
import math
from copy import copy, deepcopy

class InputData(object):
    """
        Feed in input data here and whatnot...
    """

    def __init__(self): 
        """
        """
        Nproc = 16 
        self.DEBUG = False  # For special debugging options, see below

        # File paths
        data_path = '../'
        # traj data
        self.LR_const_params = data_path+'LR_const_params.npz'

        # basic samp and opt options
        # model options
        self.flag_ind_res = False
        self.flag_cap = False
        self.flag_fix_inner_cap = False
        # step variables
        self.NSweep_p_Samp = 200 # sample in multiples of this number of sweeps
        self.NSamp_eq = 100 # throw out this many minimization sweeps before collecting statistics
        self.NSamp_max = 1000 # max number of samples before quitting, or moving to next lambda for slow_growth

        # output
        self.outfnm = 'output_variables'

        if (self.DEBUG):
            self.NSweep_p_Samp = 2
            self.NSweep_eq = 0
            self.dlamb0 = 0.10
            self.Ejumpfrac = 100.0
            self.Estuckfrac = 0.01
            self.Econvfrac = 0.50
            self.EWtol = 0.01

    def _save_output(self):
        np.savez(self.outfnm+'_ID', bin_ctrs_fnm=self.bin_ctrs_fnm, bin_ctrs_CG_fnm=self.bin_ctrs_CG_fnm, dtraj_AA_fnm=self.dtraj_AA_fnm, dtraj_CG_fnm=self.dtraj_CG_fnm, tau=self.tau, prior=self.prior, nts=self.nts, nmfpt=self.nmfpt, nmfptb=self.nmfptb, EQ_ref=self.EQ_ref, NSweep_p_Samp=self.NSweep_p_Samp, NSweep_eq=self.NSweep_eq, lamb0=self.lamb0, dlamb0=self.dlamb0, beta0C=self.beta0C)





