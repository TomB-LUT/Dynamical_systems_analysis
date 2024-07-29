from MonteCarloSim import MonteCarloSim
import numpy as np
from math import pi, sqrt, cos

class LiNonDim2023(MonteCarloSim):

    def __init__(self):
        self.par = self.set_par()
        self.tspan = self.t_span()
        self.f = self.__class__.f 
        self.init_cond = self.set_IC()
        self.parToSave = self.par_to_save()
    
    def set_par(self):
        rng = np.random.default_rng()
        m1 = 0.2
        m2 = 0.05
        k1 = 2_000
        k2 = 2000
        k3 = 2000
        l1 = 0.08
        l2 = 0.07
        l3 = 0.04
        h = 0.072
        #h = rng.uniform(0.9*h, 1.1*h)
        c1 = 0.1
        #c1 = rng.uniform(0.9*c1, 1.1*c1)
        c2 = 0.1
        #c2 = rng.uniform(0.9*c2, 1.1*c2)
        #amp = 6 
        amp = rng.uniform(1,10)
        induc = 0.005
        res = 10
        alpha = 0.8
        #alpha = rng.uniform(0.9*alpha_r, 1.1*alpha_r)
        omega1 = sqrt(k1/m1)
        omega2 = sqrt(k2/m2)
        omega3 = sqrt(k3/m2)
        lambd = m1/m2
        mi1 = k1/k2
        mi2 = k1/k3
        ni1_r = h/l1
        ni1 = rng.uniform(0.95*ni1_r, 1.05*ni1_r)
        ni2_r = h/l3
        ni2 = rng.uniform(0.95*ni2_r, 1.05*ni2_r)
        phi1_r = c1/(m1*omega1)
        phi1 = rng.uniform(0.95*phi1_r, 1.05*phi1_r)
        phi2_r = c2/(m2*omega1)
        phi2 = rng.uniform(0.95*phi2_r, 1.05*phi2_r) 
        #omega = omegaHz/omega1
        omega = rng.uniform(0,1)
        P = amp/(l1*omega1*omega1)
        #P = rng.uniform(0,1)
        ro_r = alpha/(m1*l1*omega1*omega1)
        ro = rng.uniform(0.95*ro_r, 1.05*ro_r)
        theta = res/(induc*omega1)
        epsilon_r = alpha*l1/(induc)
        epsilon = rng.uniform (0.95*epsilon_r, 1.05*epsilon_r)

        return (omega, P, lambd, ni1, ni2,  phi1, phi2, mi1, mi2, theta, epsilon, ro)
    
    def par_to_save(self):
        return (0,1)

    @staticmethod
    def f(y, t, par):  
        omega, P,  lambd, ni1, ni2,  phi1, phi2, mi1, mi2, theta, epsilon, ro = par
        return np.array( 
                        [ y[1], 
                        -y[0]*(1 - 1/sqrt( y[0]*y[0]+ni1*ni1*(1-y[2])*(1-y[2]) ) )-phi1*y[1] + (P*cos(omega*t)) - ro*y[4], # <----Jakie omega tu jest?
                        y[3],
                        -lambd/mi1*y[2] - 2*lambd/mi2*y[2]*( 1- 1/( sqrt(1+ni2*ni2*y[2]*y[2]) ) ) + lambd*(1-y[2])*( 1 - 1/( sqrt(y[0]*y[0]+ni1*ni1*(1-y[2])*(1-y[2])) ) ) - phi2*y[3], #- (-amp*cos(omega*t)),
                        -theta*y[4]+epsilon*y[1] ] )
    

    def t_span(self):
        period_of_system = 2*pi/self.par[0]
        step = 1500
        periods_no = 500
        t0 = 0 
        return (t0,periods_no,step,period_of_system)
    
    def set_IC(self):   
        rng = np.random.default_rng()
        x10 = rng.uniform(-1, 1)
        x20 = rng.uniform(-1, 1)
        return [x10, x20, 0.0, 0.0, 0.0]




#sys = (LiNonDim2023())
#print(sys.__dict__)
    