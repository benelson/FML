import numpy as np
import rebound

class Observation:
    t = None
    rv = None
    err = None
    mStar = None
    nPoints = 0
    error = 0

class FakeObservation(Observation):
    def __init__(self, planets, nPoints=30, error=1., tmax=1., mStar=1.):
        """
            Generates fake observations. 
        """
        AUperTwoPiYearToMetersPerSecond = 29805.672
        self.nPoints = nPoints
        self.error = error
        self.mStar = mStar
        sim = rebound.Simulation()
        sim.add(m=self.mStar)

        for planet in planets:
#            input = {"a": planet[0], "m": planet[1], "e": planet[2], "omega": planet[3], "M": planet[4]}
#            sim.add(sim.particles[0], **input)
            sim.add(primary=sim.particles[0],**planet)
        sim.move_to_com()
        
        self.t = np.sort(np.random.uniform(0.,tmax,self.nPoints))
        self.rv = np.zeros(nPoints)
        self.err = np.array([error] * nPoints)
        for i, t in enumerate(self.t * 2. * np.pi / 365.2425):
            sim.integrate(t)
            self.rv[i] = sim.particles[0].vx * AUperTwoPiYearToMetersPerSecond + np.random.normal(0., error)
