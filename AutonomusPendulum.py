from MonteCarloSim import MonteCarloSim
import numpy as np
from math import sin


class AutonomusPendulum(MonteCarloSim):

    def __init__(self):
        self.par = self.set_par()
        self.tspan = self.t_span()
        self.f = self.__class__.f 
        self.init_cond = self.set_IC()

    def set_par(self):
        rng = np.random.default_rng()
        omega = 1
        return (omega,)
    
    def par_to_save(self):
        return None
    
    @staticmethod
    def f(y, t, par):   
        omega, = par
        return np.array(
                        [y[1],
                        -omega*sin(y[0])] )
    
    def t_span(self):
        t0 = 0 
        tf = 300
        dt = 0.01
        return (t0,tf,dt)
    
    def set_IC(self):
        rng = np.random.default_rng()
        x1 = rng.uniform(-1,1)
        x2 = rng.uniform(-1,1)
        return [x1, x2]