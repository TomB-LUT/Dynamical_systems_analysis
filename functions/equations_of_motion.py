import numpy as np
from math import sqrt, cos, sin
#import warnings

def run_time_error_wrapper(func):
    def wrap_func(*args, **kwargs):
        try:
            results = func(*args, **kwargs)
        except RuntimeWarning:
            print(f'Run time warning encountered in the {__name__}')
        return results 
    return wrap_func

def derivs_duff(y, t:float, par, omega = 10): 

    omega, m, kx, lo, h, amp, induc, alpha, c, res = par
    dydt =[ y[1], 
            -2*kx/m* y[0] + 2*kx/m*lo*(1/sqrt(y[0]*y[0]+h*h))*y[0] + alpha/m*y[3] - c/m*y[1] + amp*cos(omega*t),
            y[3], 
            -alpha/induc*y[1] - res/induc*y[3] ]

    return np.array(dydt)


def derivs_duff_assym(y, t:float, par): 

    omega, m, kx, lo, h1, h2, amp, induc, alpha, c, res = par

    dydt =[ y[1], 
            -2*kx/m* y[0] + kx/m*lo*( (1/sqrt(y[0]*y[0]+h1*h1)) + (1/sqrt(y[0]*y[0]+h2*h2)) )*y[0] + alpha/m*y[3] - c/m*y[1] + amp*cos(omega*t),
            y[3], 
            -alpha/induc*y[1] - res/induc*y[3] ]

    return np.array(dydt)

def derivs_duff_assym_nonDim(y, t, par): 

    omega, omega1, beta, sigma1, sigma2, bigP, theta, phi, epsilon, lambd = par
    tal = t*omega
    dydt =[ y[1], 
            -2*beta* y[0] + beta*y[0]*( (1/sqrt(y[0]*y[0]+sigma1*sigma1)) +  ( (1/sqrt(y[0]*y[0]+sigma2*sigma2)) ) )  + theta*y[2] - phi*y[1] + bigP*cos(tal),
            -epsilon*y[1] - lambd*y[2] ]

    return np.array(dydt)

def derivs_duff_assym_nonDim_2(y, t, par): 

    BigOmega, sigma1, sigma2, P, theta, phi, epsilon, lambd = par
    dydt =[ y[1], 
            -2* y[0] + y[0]*( (1/sqrt(y[0]*y[0]+sigma1*sigma1)) +  ( (1/sqrt(y[0]*y[0]+sigma2*sigma2)) ) )  + theta*y[2] - phi*y[1] + P*cos(BigOmega*t),
            -epsilon*y[1] - lambd*y[2] ]

    return np.array(dydt)

def derivs_Li_MSSP_2023(y, t, par):
    omega, m1, m2,  k1, k2, k3, l1, l2, l3, h, c1, c2, induc, alpha, res, amp = par
    dydt =[ y[1], 
           -k1/m1*y[0]*(1 - l1/sqrt( y[0]*y[0]+(h-y[2])*(h-y[2]) ) )-c1/m1*y[1] - (-amp*cos(omega*t)) - alpha/m1*y[5],
           y[3],
           -k2/m2*y[2] - 2*k3/m2*y[2]*( 1- l3/( sqrt(l3*l3+y[2]*y[2]) ) ) + k1/m2*(h-y[2])*( 1 - l1/( sqrt(y[0]*y[0]+(h-y[2])*(h-y[2])) ) ) - c2/m2*y[3], #- (-amp*cos(omega*t)),
           y[5],
           -res/induc*y[5]+alpha/induc*y[1] ]
    return np.array(dydt) 

#@run_time_error_wrapper
def derivs_Li_MSSP_2023_nondim(y, t, par):

    omega, P,  lambd, ni1, ni2,  phi1, phi2, mi1, mi2, theta, epsilon, ro = par
    dydt =[ y[1], 
        -y[0]*(1 - 1/sqrt( y[0]*y[0]+ni1*ni1*(1-y[2])*(1-y[2]) ) )-phi1*y[1] + (P*cos(omega*t)) - ro*y[4], # <----Jakie omega tu jest?
        y[3],
        -lambd/mi1*y[2] - 2*lambd/mi2*y[2]*( 1- 1/( sqrt(1+ni2*ni2*y[2]*y[2]) ) ) + lambd*(1-y[2])*( 1 - 1/( sqrt(y[0]*y[0]+ni1*ni1*(1-y[2])*(1-y[2])) ) ) - phi2*y[3], #- (-amp*cos(omega*t)),
        -theta*y[4]+epsilon*y[1] ]
    
    return np.array(dydt) 

def derivs_Costa_2024(y, tal, par):

    omega, sigma, c1, c2, a1, a2, b1, b2, ksi1, ksi2, k1, k2, ro, omegaS, phi1, phi2 = par

    #exct = sigma*sin(omega*tal)
    #exct_dot = sigma*omega*cos(omega*tal)
    exct_dot_dot = -sigma*omega*omega*sin(omega*tal)

    dydt = [ y[1],
            -2*c1*y[1]+2*c2*(y[3]-y[1])-(1+a1)*y[0]-b1*y[0]*y[0]*y[0]+ro*omegaS*omegaS*(y[2]-y[0])+ksi1*y[4]-ksi2*y[5] - exct_dot_dot,
            y[3],
            (-2*c2*(y[3]-y[1])-a2*y[2]-b2*y[2]*y[2]*y[2]-ro*omegaS*omegaS*(y[2]-y[0])+ksi2*y[5] - exct_dot_dot   )/ro,
            -phi1*y[4]-k1*y[1],
            -phi2*y[5]-k2*(y[3]-y[1]) ]

    return np.array(dydt) 

def derivs_Wang_NonLinDyn_2024(y, t, par):

    omega, c1, c2, k, ksi, a, b, gamma, f, d= par 

    g_fun = 2*c2*y[1] + k*(y[0]+d) if y[0] <= -d else 0 
    #print('t: {}, g_fun: {}'.format(t, g_fun))

    dydt = [ y[1], 
            -2*c1*y[1] + y[0] - b*y[0]*y[0] - gamma*y[0]*y[0]*y[0] - g_fun + ksi*ksi*y[2] + f*cos(omega*t),
            -a*y[2] - y[1]    ]

    return np.array(dydt)

def derivs_Schen_2023(y, t, par):

    c1, c2, k1, k3, a, lamb, ke, f1, f3, ua, r = par

    dydt = [ y[1], 
            -2*c1*y[1] - k1*y[0] - k3*y[0]*y[0]*y[0] - 2*c2*a*lamb*(y[1]-y[3]) - a*lamb*lamb*(y[0]-y[2]) - ke*y[4] + f1*ua*y[1] + f3/ua*y[1]*y[1]*y[1],
            y[3], 
            -2*c2*lamb*(y[3]-y[1]) - lamb*lamb*(y[2]-y[0]), 
            y[1] - y[4]/r
    ]

def derivs_Li_MSSP_2023_fixed_point(y, t, par):
    omega, m1, m2,  k1, k2, k3, l1, l2, l3, h, c1, c2, induc, alpha, res, amp = par
    dydt =[ y[1], 
           -k1/m1*y[0]*(1 - l1/sqrt( y[0]*y[0]+(h-y[2])*(h-y[2]) ) )-c1/m1*y[1] - alpha/m1*y[5],
           y[3],
           -k2/m2*y[2] - 2*k3/m2*y[2]*( 1- l3/( sqrt(l3*l3+y[2]*y[2]) ) ) + k1/m2*(h-y[2])*( 1 - l1/( sqrt(y[0]*y[0]+(h-y[2])*(h-y[2])) ) ) - c2/m2*y[3], #- (-amp*cos(omega*t)),
           y[5],
           -res/induc*y[5]+alpha/induc*y[1] ]
    return np.array(dydt)

def derivs_duff_2dof(y, t:float, omega = 10): 

    #Do tupli to i rozsyłaj?
    m = 0.2
    kx = 1500
    lo = 0.113738
    h = 0.961587*0.113738
    amp = 6.2
    induc = 1.463
    alpha = 0#30
    c = 0.01 * 2*sqrt(kx*m) # X % of critical damping
    res = 2200
    print(c)

    #dydt = np.zeros((4,1))
    dydt = [0,0,0,0]
    dydt =[ y[1], 
            -2*kx/m* y[0] + 2*kx/m*lo*(1/sqrt(y[0]*y[0]+h*h))*y[0] - alpha/m*y[3] - c/m*y[1] + amp*cos(omega*t),
            y[3], 
            alpha/induc*y[1] - res/induc*y[3] ]

    return np.array(dydt)

def derivs_lorenz(y:list, t:float)->list: 

    sigma = 10
    beta = 8/3 
    rho = 28 

    #dydt = np.zeros((4,1))
    dy = [0,0,0]

    dy[0] = sigma * (y[1] - y[0])
    dy[1] = y[0] * (rho - y[2]) - y[1] 
    dy[2] = y[0] * y[1] - beta *y[2]
    

    #print(dydt)

    return np.array(dy)

def derivs_duff_simple(y, t):

    dydt = [0,0]
    dydt = [y[1],
            -0.1*y[1]+y[0]-pow(y[0],3)]

    
    return np.array(dydt)
    

