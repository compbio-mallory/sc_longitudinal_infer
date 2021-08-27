#!/usr/bin/env python3

import numpy as np
import math
from scipy.stats import beta, truncnorm
import bottleneck as bn

try:
    from libs.CRP import CRP
except ImportError:
    from CRP import CRP


# ------------------------------------------------------------------------------
# LEARNING ERROR RATES - NO NANs
# ------------------------------------------------------------------------------

class CRP_errors_learning(CRP):
    def __init__(self, data, DP_alpha=1, param_beta=[1, 1], \
                FP_mean=0.001, FP_sd=0.0005, FN_mean=0.25, FN_sd=0.05):
        super().__init__(data, DP_alpha, param_beta, FN_mean, FP_mean)
        # Error rate prior
        FP_trunc_a = (0 - FP_mean) / FP_sd
        FP_trunc_b = (1 - FP_mean) / FP_sd
        self.FP_prior = truncnorm(FP_trunc_a, FP_trunc_b, FP_mean, FP_sd)
        # self.FP_prior = beta(1, (1 - FP_mean) / FP_mean)
        self.FP_sd = np.array([FP_sd * 0.5, FP_sd, FP_sd * 1.5])
        # print("FP_mean: {}\n".format(FP_mean))
        # print("FN_mean: {}\n".format(FN_mean))

        FN_trunc_a = (0 - FN_mean) / FN_sd
        FN_trunc_b = (1 - FN_mean) / FN_sd
        self.FN_prior = truncnorm(FN_trunc_a, FN_trunc_b, FN_mean, FN_sd)
        # self.FN_prior = beta(1, (1 - FN_mean) / FN_mean)
        self.FN_sd = np.array([FN_sd * 0.5, FN_sd, FN_sd * 1.5])

       # our additon
       # need to decides true value on MR_mean and  MR_sd after experiment
       # Following trun norm for missing rate as that appliled to FP and FN
        MR_mean = 0.001 
        MR_sd = 0.0005
        MR_trunc_a = (0 - MR_mean) / MR_sd
        MR_trunc_b = (1 - MR_mean) / MR_sd
        self.MR_sd = np.array([MR_sd * 0.5, MR_sd, MR_sd * 1.5])
        self.MR_prior  = truncnorm(MR_trunc_a, MR_trunc_b, MR_mean, MR_sd)

    def __str__(self):
        out_str = '\nDPMM with:\n' \
            f'\t{self.cells_total} cells\n\t{self.muts_total} mutations\n' \
            f'\tlearning errors\n' \
            '\n\tPriors:\n' \
            f'\tparams.:\tBeta({self.p},{self.q})\n' \
            f'\tCRP a_0:\tGamma({self.DP_a_gamma[0]:.2f},{self.DP_a_gamma[1]})\n' \
            f'\tFP:\t\ttrunc norm({self.FP_prior.args[2]},{self.FP_prior.args[3]})\n' \
            f'\tFN:\t\ttrunc norm({self.FN_prior.args[2]},{self.FN_prior.args[3]})\n'
        return out_str


    def get_lprior_full(self):
        return super().get_lprior_full() \
            + self.FP_prior.logpdf(self.FP) + self.FN_prior.logpdf(self.FN) + self.MR_prior.logpdf(self.MR)

    # ------------------------------------------------------------------------------
    # our additon to get likelihood breaking down into 6 different equations
    # ------------------------------------------------------------------------------
    # def get_ll(self, FP, FN, MR, x, theta):
    #     if(math.isnan(x) and theta == 0):
    #         return MR
    #     elif(math.isnan(x) and theta == 1):
    #         return MR
    #     elif(x == 0 and theta == 0):
    #         return 1 - FN - MR
    #     elif(x == 0 and theta == 1):
    #         return FP
    #     elif(x == 1 and theta == 0):
    #         return FN
    #     else:
    #         return 1 - FP - MR 


    def update_error_rates(self):
        # print("calling update_error_rates !!")
        self.FP, FP_count = self.MH_error_rates('FP')
        self.FN, FN_count = self.MH_error_rates('FN')
        self.MR, err_count = self.MH_error_rates('MR')
        # print("n: {}\n".format(FP_count))
        # print("FN_count: {}\n".format(FN_count))
        return FP_count, FN_count


    def get_ll_full_error(self, FP, FN, MR):
    # def get_ll_full_error(self, FP, FN):
        sum = 0
        # print("calling get_ll_full_error !!")
        par = self.parameters[self.assignment]
        # print("parameter", par[0])
        # print("par", type(par))
        # ll_FN = par * (1 - FN) ** self.data * FN ** (1 - self.data)
        ll_FN = par * (1 - FN - MR) ** self.data * FN ** (1 - self.data)
        # ll_FP = (1 - par) * (1 - FP) ** (1 - self.data) * FP ** self.data
        ll_FP = (1 - par) * (1 - FP - MR) ** (1 - self.data) * FP ** self.data
        ll_full = np.log(ll_FN + ll_FP)
        # adding MR part
        for i in range(self.muts_total):
            for j in range(self.cells_total):
                if math.isnan(self.data[i][j]):
                    sum = sum + math.log(MR)
        return bn.nansum(ll_full) + sum
        # return bn.nansum(ll_full) 

    # ------------------------------------------------------------------------------
    # our additon to get likelihood breaking down into 6 different equations
    # ------------------------------------------------------------------------------

    # def get_ll_full_error(self, FP, FN, MR):
    #     temp_ll = 1
    #     par = self.parameters[self.assignment]
    #     for i in range(self.muts_total):
    #         for j in range(self.cells_total):
    #            temp_ll = temp_ll * math.log(self.get_ll(FP, FN, MR, self.data[i][j], par[i][j]))

    #     return temp_ll

    def MH_error_rates(self, error_type):
        # Set error specific values
        if error_type == 'FP':
            old_error = self.FP
            prior = self.FP_prior
            stdevs = self.FP_sd
        
        # our additon
        elif error_type == 'MR':
            old_error = self.MR
            prior = self.MR_prior
            stdevs = self.MR_sd

        else:
            old_error = self.FN
            prior = self.FN_prior
            stdevs = self.FN_sd

        # Get new error from proposal distribution
        std = np.random.choice(stdevs)
        a = (0 - old_error) / std
        b = (1 - old_error) / std
        try:
            new_error = truncnorm.rvs(a, b, loc=old_error, scale=std)
        except FloatingPointError:
            new_error = truncnorm.rvs(a, np.inf, loc=old_error, scale=std)

        # Calculate transition probabilitites
        new_p_target = truncnorm \
            .logpdf(new_error, a, b, loc=old_error, scale=std)
        a_rev, b_rev = (0 - new_error) / std, (1 - new_error) / std
        old_p_target = truncnorm \
            .logpdf(old_error, a_rev, b_rev, loc=new_error, scale=std)

        # Calculate likelihood
        if error_type == 'FP':
            # new_ll = self.get_ll_full_error(new_error, self.FN)
            # old_ll = self.get_ll_full_error(old_error, self.FN)
            new_ll = self.get_ll_full_error(new_error, self.FN, self.MR)
            # print("new_ll", new_ll)
            old_ll = self.get_ll_full_error(old_error, self.FN, self.MR)
            # print("old_ll", old_ll)
        # our additon
        elif error_type == 'MR':
            new_ll = self.get_ll_full_error(self.FP, self.FN, new_error)
            old_ll = self.get_ll_full_error(self.FP, self.FN, old_error)

        else:
            # new_ll = self.get_ll_full_error(self.FP, new_error)
            # old_ll = self.get_ll_full_error(self.FP, old_error)
            new_ll = self.get_ll_full_error(self.FP, new_error, self.MR)
            # print("new_ll", new_ll)
            old_ll = self.get_ll_full_error(self.FP, old_error, self.MR)
            # print("old_ll", old_ll)

        # Calculate priors
        new_prior = prior.logpdf(new_error)
        old_prior = prior.logpdf(old_error)

        # Calculate MH decision treshold
        A = new_ll + new_prior - old_ll - old_prior + old_p_target - new_p_target
        # print("new_error", new_error)
        # print("old_error", old_error)
        if np.log(np.random.random()) < A:
            return new_error, [1, 0]

        return old_error, [0, 1]


if __name__ == '__main__':
    print('Here be dragons...')