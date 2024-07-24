from MonteCarloSim import MonteCarloSim
import numpy as np

class DuffingFixedPoint(MonteCarloSim):
    
    def __init__(self):
        self.par = self.set_par()
        self.tspan = self.t_span()
        self.f = self.__class__.f 
        self.init_cond = self.set_IC()
    
    def set_par(self):
        p1 = 0.1
        return (p1,)
    
    @property
    def raw_f(self):
        return '''
p1, = par        
dydt = [y[1], -p1*y[1]+y[0]-pow(y[0],3)]  
                '''
    
    @staticmethod
    def f(y, t, par):   
        p1, = par
        return np.array(
                        [y[1],
                        -p1*y[1]+y[0]-pow(y[0],3)] )

    def t_span(self):
        t0 = 0 
        tf = 3000
        dt = 0.01
        return (t0,tf,dt)
    
    def set_IC(self):
        rng = np.random.default_rng()
        x1 = rng.uniform(-2,2)
        #x1 = 1.6
        x2 = rng.uniform(-2,2)
        return [x1, x2]
    
