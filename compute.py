import numpy as np
import rebound
import copy
from collections import OrderedDict

def within_range(x, lower, upper):
    x = [ i - 2.*np.pi if i>upper else i for i in x ]
    x = [ i + 2.*np.pi if i<lower else i for i in x ]
    return x

def sma_and_mass(theta, mStar):
    planets = [theta[x:x+5] for x in range(0, len(theta), 5)]
    mSum = mStar
    MetersPerSecondtoAUperTwoPiYear = 1./( 4740.57581 * 2 * np.pi)

    planets_am = []

    for planet in planets:
        Pinyears = planet[0] / 365.2425
        K = planet[1]
        e = planet[2]
        
        sma = 0.
        mass = 0.
        
        for aa in range(6):
            sma = (Pinyears * Pinyears * (mSum + mass)) ** (0.333333)
            mass = MetersPerSecondtoAUperTwoPiYear * K * np.sqrt( mStar * sma * (1 - e*e) )

        mSum += mass
        
        planets_am.append({"a": sma, "m": mass, "e": e, "omega": planet[3], "M": planet[4]})

    return planets_am


def rv(theta, obs):
    AUperTwoPiYearToMetersPerSecond = 29805.672
    times = obs.t * 2 * np.pi / 365.2425

    sim = rebound.Simulation()
    sim.add(m=obs.mStar)
    for planet in sma_and_mass(theta, obs.mStar):
        sim.add(primary=sim.particles[0],**planet)
    sim.move_to_com()
    
    rv = np.zeros(len(times))
    for i, t in enumerate(times):
        sim.integrate(t)
        rv[i] = sim.particles[0].vx * AUperTwoPiYearToMetersPerSecond

    return rv


def pKewM_to_pKhkl(postSamp, keys_tf, nPlanets):
    Pmax = 999999.
    Kmax = 999999.

    dic = OrderedDict()
    for i in range(len(keys_tf)):
        dic[keys_tf[i]] = []

    keys = [ key for key in postSamp ]

    for i in range(nPlanets):
        dic[keys_tf[i*5+0]] = postSamp[keys[i*5+0]]
        dic[keys_tf[i*5+1]] = postSamp[keys[i*5+1]]
        
        dic[keys_tf[i*5+2]] = np.array(postSamp[keys[i*5+2]])\
                                * np.sin(np.array(postSamp[keys[i*5+3]])) #esinw
        dic[keys_tf[i*5+3]] = np.array(postSamp[keys[i*5+2]])\
                                * np.cos(np.array(postSamp[keys[i*5+3]])) #ecosw
        
        dic[keys_tf[i*5+4]] = within_range( np.array(postSamp[keys[i*5+3]]) + np.array(postSamp[keys[i*5+4]]),\
                                          0., 2.*np.pi) #w+M

    for i in range(len(keys_tf) - nPlanets*5):
        dic[keys_tf[nPlanets*5+i]] = postSamp[keys[nPlanets*5+i]]
    
    return dic


def pKhkl_to_pKewM(importSamp, keys_tf, nPlanets):
    dic = OrderedDict()
    for i in range(len(keys_tf)):
        dic[keys_tf[i]] = []
    keys = [ key for key in importSamp ]

    for i in range(nPlanets):
        esinw = np.array(importSamp[keys[i*5+2]])
        ecosw = np.array(importSamp[keys[i*5+3]])
        e = np.sqrt(esinw*esinw + ecosw*ecosw)

        omegas = np.arctan(esinw/ecosw)
        for j,w in enumerate(omegas):
            if ecosw[j]<0.: omegas[j] = np.pi + w
            if(esinw[j]<0. and ecosw[j]>0.): omegas[j] = 2.*np.pi + w

        omegas = within_range(omegas, 0., 2.*np.pi)

        dic[keys_tf[i*5+0]] = importSamp[keys[i*5]]
        dic[keys_tf[i*5+1]] = importSamp[keys[i*5+1]]
        dic[keys_tf[i*5+2]] = e
        dic[keys_tf[i*5+3]] = omegas
        tmp = np.array(importSamp[keys[i*5+4]]) - np.array(omegas)
        dic[keys_tf[i*5+4]] = within_range(tmp, 0., 2.*np.pi)

    for i in range(len(keys_tf) - nPlanets*5):
        dic[keys_tf[nPlanets*5+i]] = importSamp[keys[nPlanets*5+i]]

    return dic


def pKhkl_to_amewM(importSamp, keys_tf, nPlanets, mStar):
    MetersPerSecondtoAUperTwoPi = 1./( 4740.57581 * 2 * np.pi)

    dic = OrderedDict()
    for i in range(len(keys_tf)):
        dic[keys_tf[i]] = []
    keys = [ key for key in importSamp ]

    mSum = mStar

    for i in range(nPlanets):
        sma = 0.
        mass = 0.

        Pinyears = np.array(importSamp[keys[i*5+0]]) / 365.2425
        K = np.array(importSamp[keys[i*5+1]])
        esinw = np.array(importSamp[keys[i*5+2]])
        ecosw = np.array(importSamp[keys[i*5+3]])
        e = np.sqrt(esinw*esinw + ecosw*ecosw)

        omegas = np.arctan(esinw/ecosw)*180./np.pi
        for j,w in enumerate(omegas):
            if ecosw[j]<0.: omegas[j] = 180. + w
            if(esinw[j]<0. and ecosw[j]>0.): omegas[j] = 360. + w

        for aa in range(6):
            sma = (Pinyears * Pinyears * (mSum + mass)) ** (0.333333)
            mass = MetersPerSecondtoAUperTwoPi * K * np.sqrt( mStar * sma * (1 - e*e) )

        mSum += mass

        dic[keys_tf[i*5+0]] = sma
        dic[keys_tf[i*5+1]] = mass

        dic[keys_tf[i*5+2]] = e
        dic[keys_tf[i*5+3]] = omegas
        dic[keys_tf[i*5+4]] = importSamp[keys[i*5+4]] - omegas

    for i in range(len(keys_tf) - nPlanets*5):
        dic[keys_tf[nPlanets*5+i]] = importSamp[keys[nPlanets*5+i]]

    return dic


def pKhkl_to_amewl_old(postSamp, keys_tf, nPlanets, mStar):
    MetersPerSecondtoAUperTwoPi = 1./( 4740.57581 * 2 * np.pi)
    
    dic = OrderedDict()
    for i in range(len(keys_tf)):
        dic[keys_tf[i]] = []

    keys = [ key for key in postSamp ]
    
    mSum = mStar
    
    for i in range(nPlanets):
        sma = 0.
        mass = 0.
        
        Pinyears = np.array(postSamp[keys[i*5+0]]) / 365.2425
        K = np.array(postSamp[keys[i*5+1]])
        esinw = np.array(postSamp[keys[i*5+2]])
        ecosw = np.array(postSamp[keys[i*5+3]])
        e = np.sqrt(esinw*esinw + ecosw*ecosw)
        
        omegas = np.arctan(esinw/ecosw)*180./np.pi
        for j,w in enumerate(omegas):
            if ecosw[j]<0.: omegas[j] = 180. + w
            if(esinw[j]<0. and ecosw[j]>0.): omegas[j] = 360. + w
        
        for aa in range(6):
            sma = (Pinyears * Pinyears * (mSum + mass)) ** (0.333333)
            mass = MetersPerSecondtoAUperTwoPi * K * np.sqrt( mStar * sma * (1 - e*e) )
            
        mSum += mass
        
        dic[keys_tf[i*5+0]] = sma
        dic[keys_tf[i*5+1]] = mass
        
        dic[keys_tf[i*5+2]] = e
        dic[keys_tf[i*5+3]] = omegas
        dic[keys_tf[i*5+4]] = postSamp[keys[i*5+4]]
        
    for i in range(len(keys_tf) - nPlanets*5):
        dic[keys_tf[nPlanets*5+i]] = postSamp[keys[nPlanets*5+i]]
    
    return dic


def matrix_info(dic):
    vals = []
    for x in dic:
        vals.append(dic[x])
    covMatrixG = np.cov(vals)
    
    s,logDetSigmaG = np.linalg.slogdet(covMatrixG)
    
    matrixA = np.linalg.cholesky(covMatrixG)

    mediansG = []
    for a in dic:
        mediansG.append(np.median(dic[a]))
        
    return np.array(mediansG), covMatrixG, matrixA, logDetSigmaG


def physical_model(mStar, planets, times):
    times = times * 2 * np.pi / 365.2425
    AUperTwoPiYearToMetersPerSecond = 4743.72 * 2 * np.pi
    
    sim = rebound.Simulation()
    sim.add(m=mStar)
    for p in planets:
        sim.add(a=p['a'], m=p['m'], e=p['e'], omega=p['w'], M=p['M'])
    sim.move_to_com()
    
    rvs_true = np.zeros(len(times))    
    for i, t in enumerate(times):
        sim.integrate(t)
        rvs_true[i] = -sim.particles[0].vy * AUperTwoPiYearToMetersPerSecond

    return rvs_true


def simulated_obs(mStar, planets, times):
    rvs_err = 1.0 + 1.0 * np.random.rand(len(times))
    rvs_obs = physical_model(mStar, planets, times) + rvs_err * np.random.randn(len(times))
    return rvs_obs, rvs_err


def model_evaluation(mStar, planets, data):
    times = data['times'] * 2 * np.pi / 365.2425
    rvs_true = physical_model(mStar, planets, times)
    rvs_obs = copy.deepcopy(data['rvs'])
    rv_errs = data['errs']
    rv_errs_cov = np.diag(rv_errs*rv_errs)
    ind = data['obs_indices']

    diff = rvs_real - rvs_obs
    chisq = np.dot( diff, np.linalg.solve(rv_errs_cov, diff) )
    return -0.5 * chisq   


'''
def model_evaluation(mStar, planets, offs, data):
    AUperTwoPiYearToMetersPerSecond = 4743.72 * 2 * np.pi
    sim = rebound.Simulation()
    sim.add(m=mStar)
    for p in planets:
        sim.add(a=p['a'], m=p['m'], e=p['e'], omega=p['w'], M=p['M'])
    sim.move_to_com()
    times = data['times'] * 2 * np.pi / 365.2425
    
    rvs_real = np.zeros(len(times))
    rvs_obs = copy.deepcopy(data['rvs'])
    rv_errs = data['errs']
    rv_errs_cov = np.diag(rv_errs*rv_errs)
    ind = data['obs_indices']
    
    for i, t in enumerate(times):
        sim.integrate(t)
        rvs_real[i] = -sim.particles[0].vy * AUperTwoPiYearToMetersPerSecond
        rvs_obs[i] = rvs_obs[i] - offs[ind[i]]
        
    diff = rvs_real - rvs_obs
    chisq = np.dot( diff, np.linalg.solve(rv_errs_cov, diff) )
    return -0.5 * chisq
'''
