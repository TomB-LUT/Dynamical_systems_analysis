from MonteCarloSim import MonteCarloSim
import numpy as np

class RooslerSystem(MonteCarloSim):

    def __init__(self):
        self.par = self.set_par()
        self.tspan = self.t_span()
        self.f = self.__class__.f 
        self.init_cond = self.set_IC()

    def set_par(self):
        a = 0.2
        b = 9.0
        c = 48.0 
        return (a,b,c)
    
    def par_to_save(self):
        return None
    
    @staticmethod
    def f(y, t, par):   
        a,b,c = par
        return np.array(
                        [-y[1]-y[2],
                        y[0]+a*y[1],
                        b+y[2]*(y[0]-c)] )
    
    def t_span(self):
        t0 = 0 
        tf = 300
        dt = 0.01
        return (t0,tf,dt)
    
    def set_IC(self):
        rng = np.random.default_rng()
        x1 = rng.uniform(-85,85)
        x2 = rng.uniform(-80,60)
        x3 = rng.uniform(-30,270)
        return [x1, x2, x3]