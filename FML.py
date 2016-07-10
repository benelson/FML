import numpy as np
import rebound

import compute
import create

from scipy.stats import truncnorm
from math import erf

from collections import OrderedDict
from ipywidgets import FloatProgress
from IPython.display import display


class FML:
    samples = None
    obs = None
    nPlanets = 0.
    nImportSamp = 0
    scale = 0.

    logAvg = 0.
    f_MCMC = 0.
    logFML = 0.


def lnlike(theta, obs):
    '''
    The log-likelihood function. Observational model assume RVs have independent Gaussian error bars.
        
        Parameters
        ----------
        theta : 1-d list of Keplerian orbital parameters (p, K, e, w, M)
                for every modeled planet
        obs   : observations object with radial velocity measurement (.rv)
                and uncertainty (.err) attributes
    '''
    model = compute.rv(theta, obs)
    inv_sigma2 = 1.0/(obs.err * obs.err)
    return -0.5*( np.sum( (obs.rv - model)*(obs.rv - model)*inv_sigma2 - np.log(inv_sigma2) ) )

def lnprior(theta):
    '''
    The log-prior function.
        
         Parameters
         ----------
         theta : 1-d list of Keplerian orbital parameters (p, K, e, w, M)
                 for every modeled planet
    '''
    planets = [theta[x:x+5] for x in range(0, len(theta), 5)]
    nPlanets = len(planets)
    thetaT = np.transpose(planets)
    p, K, e, w, M = thetaT[0], thetaT[1], thetaT[2], thetaT[3], thetaT[4]
    
    P0, K0 = 1., 1.
    lnp = 0.
    
    for i in range(nPlanets):
        if p[i] < 99999.and K[i] < 99999.and 0. <= e[i] < 1. \
        and 0. <= abs(w[i]) <= 2*np.pi and 0. <= abs(M[i]) <= 2*np.pi:
            lnp += -np.log(1.+p[i]/P0) - np.log(1+K[i]/K0) - 2.*np.log(2.*np.pi)
        else:
            return -np.inf
    return lnp
    
def lnpost(theta, obs):
    '''
    The log-posterior function. Sum of log-likelihood and log-prior.
        
        Parameters
        ----------
        theta : 1-d list of Keplerian orbital parameters (p, K, e, w, M)
                for every modeled planet
        obs   : observations object with radial velocity measurement (.rv)
                and uncertainty (.err) attributes
    '''
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, obs)

class computeFML(FML):
    def __init__(self, samples, obs, nPlanets=0, nImportSamps=10000, scale=1.0):
        """
        Computes fully marginalized likelihood

        Parameters
        ----------
        samples      : posterior samples from emcee
        nPlanets     : number of planets
        nImportSamps : number of importance samples
        scale        : sets the scale for the truncation of the multivariate normal, measured in "sigmas"
        """

        self.samples = samples
        self.nPlanets = nPlanets
        self.nImportSamps = nImportSamps
        self.scale = scale

        param_keys, param_IS_keys = create.dict_keys(self.nPlanets)
        postSamp, nPostSamples = create.posterior_samples_from_emcee(self.samples, param_keys)
        postSamp_pKhkl = compute.pKewM_to_pKhkl(postSamp, param_IS_keys, self.nPlanets)

        self.mediansG, self.covMatrixG, self.choleskyDecomp, self.logDetSigmaG = compute.matrix_info(postSamp_pKhkl)

        nParams = len(param_keys)
        random_values = [ truncnorm.rvs(-self.scale, self.scale, size=nParams) for i in range(self.nImportSamps) ]

        samples = [ [] for i in range(self.nImportSamps) ]
        g_samples = [ [] for i in range(self.nImportSamps) ]
        loggs = [ 0. for i in range(self.nImportSamps) ]

        # visualizes a progress bar
        f = FloatProgress(min=0, max=self.nImportSamps)
        display(f)

        print("Drawing importance samples...")

        for x in range(self.nImportSamps):
            dispersion = np.dot( self.choleskyDecomp, np.transpose(random_values[x]) )
            samples[x] = self.mediansG + dispersion
            g_samples[x] = list(samples[x])

            logg = -0.5 * (nParams*np.log(2.*np.pi) + self.logDetSigmaG + \
                    np.dot( np.transpose(np.subtract(samples[x],self.mediansG)), \
                    np.linalg.solve(self.covMatrixG, np.subtract(samples[x],self.mediansG) ) ) ) - \
                    nParams*np.log(erf(self.scale/np.sqrt(2.)))
            loggs[x] = logg
            f.value = x

        print("Done drawing importance samples!")
        print("")
    
        g_samples_T = np.transpose(g_samples)
        importSamp_dict = OrderedDict()

        for i, item in enumerate(g_samples_T):
            importSamp_dict[param_IS_keys[i]] = item


        importSamp_pKhkl_dict = compute.pKhkl_to_pKewM(importSamp_dict, param_keys, nPlanets)
        importSamp_pKewM = np.transpose([ vals for key, vals in importSamp_pKhkl_dict.items() ])

        f = FloatProgress(min=0, max=self.nImportSamps)
        display(f)

        print("Evaluating lnpost at importance samples...")

        logPosteriors = np.array([ np.nan for i in range(self.nImportSamps) ])
        for i in range(nImportSamps): 
            logPosteriors[i] = lnpost(importSamp_pKewM[i], obs)
            f.value = i
    
        print("Done evaluating lnpost!")
        print("")


        logSum = -(9.**99.)

        for i in range(self.nImportSamps):    
            diff = logPosteriors[i] - loggs[i]
            logSum = np.logaddexp(logSum, diff)
    
        self.logAvg = logSum - np.log(self.nImportSamps)
        self.f_MCMC = 0.

        postSamp_wo_keys = []
        for key in postSamp_pKhkl:
            postSamp_wo_keys.append(postSamp_pKhkl[key])
    
        postSamp_wo_keys = np.transpose(np.array(postSamp_wo_keys))
        diff = postSamp_wo_keys-self.mediansG

        for j in range(nPostSamples):

            z = np.linalg.solve(self.choleskyDecomp, diff[j])

            if all([abs(k)<=scale for k in z]):
                self.f_MCMC += 1.
            else:
                self.f_MCMC += 0.
        
        self.f_MCMC = self.f_MCMC/nPostSamples
        self.logFML = self.logAvg - np.log(self.f_MCMC)

        print("FML computed!")
        print("Done!")
