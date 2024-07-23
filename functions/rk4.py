from typing import Callable
import numpy as np


def rk4(fun:Callable, x_k:list, t_k:float, dt:float)->list: 

    f1 = fun(x_k, t_k)    
    x_arg = list(map(lambda x,y: x + (dt/2) * y, x_k, f1))
    f2 = fun(x_arg, t_k+dt/2)   
    x_arg = []
    x_arg = list(map(lambda x,y: x + (dt/2) * y, x_k, f2))
    f3 = fun(x_arg, t_k+dt/2)     
    x_arg = []
    x_arg = list(map(lambda x,y: x+dt*y, x_k, f3))
    f4 = fun(x_arg,t_k+dt) 
    rk4_list = list(map(lambda x, y, z, u: x + 2*y + 2*z + u, f1, f2, f3, f4))   
    x_k1 = list(map(lambda x, y: x + (dt/6) * y, x_k, rk4_list))

    return x_k1

def rk4_np(fun:Callable, x_k:list, t_k:float, dt:float, par:tuple, tf)->list: 

    x_k = np.array(x_k)
    f1 = fun(x_k, t_k,par)
    f2 = fun(x_k + (dt/2.0) * f1, t_k + dt/2.0, par) 
    f3 = fun(x_k + (dt/2.0) * f2, t_k + dt/2.0, par) 
    f4 = fun(x_k + dt * f3, t_k + dt, par) 
    x_k1 = x_k + (dt/6.0) * (f1 + 2.0*f2 + 2.0*f3 + f4)

    return x_k1

def rk4_np_2(fun:Callable, x_k:list, t_k:float, dt:float)->list: 

    x_k = np.array(x_k)
    f1 = dt*fun(x_k, t_k) 
    f2 = dt*fun(x_k + f1/2, t_k + dt/2.0) 
    f3 = dt*fun(x_k + f2/2, t_k + dt/2.0) 
    f4 = dt*fun(x_k + f3, t_k + dt) 
    x_k1 = x_k + (f1 + 2.0*f2 + 2.0*f3 + f4)/6.0
    
    return x_k1

def rk4_sd(fun:Callable, x_k:list, t_k:float, dt:float, par:tuple, tf,  tol=1e-6):
    print('full: ')
    y_full = rk4_np(fun, x_k, t_k, dt, par, tf)
    print('half: ')
    y_half = rk4_np(fun, x_k, t_k, dt/2, par, tf)
    print('half, half: ')
    y_half = rk4_np(fun, y_half, t_k+dt/2, dt/2, par, tf)
    print('------------')

    error = np.linalg.norm(y_half - y_full) /15.0

    return y_full, dt

   # if error < tol:
        #t += h
        #y = y_double_step
        #times.append(t)
        #ys.append(y)
        #h = h * min(2, (tol/error)**0.20)
    #else:
        #h = h * max(0.5, (tol/error)**0.20)

def rk_PP_np(fun, x_k, t_k, dt):
    x_k = np.array(x_k)
    f1 = fun(x_k, t_k) 
    f2 = fun(x_k + (dt/2.0) * f1, t_k + dt/2.0) 
    f3 = fun(x_k + (dt/2.0) * f2, t_k + dt/2.0) 
    f4 = fun(x_k + dt * f3, t_k + dt) 
    x_k1 = x_k + (dt/6.0) * (f1 + 2.0*f2 + 2.0*f3 + f4)

    return x_k1


def rk_45(fun, x_k, t_k, dt, par, tf,  tol=1e-4):
    
    if t_k+dt == t_k:
        print('Time step too small (Computer treat it as 0)')
    
    x_k = np.array(x_k)

    #a2, a3, a4, a5, a6 = 1/4, 3/8, 12/13, 1, 1/2
    #b21 = 1/4
    #b31, b32 = 3/32, 9/32
    #b41, b42, b43 = 1932/2197, -7200/2197, 7296/2197
    #b51, b52, b53, b54 = 439/216, -8, 3680/513, -845/4104
    #b61, b62, b63, b64, b65 = -8/27, 2, -3544/2565, 1859/4104, -11/40
    #c1, c3, c4, c5, c6 = 16/135, 6656/12825, 28561/56430, -9/50, 2/55
    #ch1, ch3, ch4, ch5 = 25/216, 1408/2565, 2197/4104, -1/5

    k1 = dt * fun(x_k, t_k,  par)
    k2 = dt * fun(x_k + 1/4*k1, t_k + 1/4*dt,  par)
    k3 = dt * fun(x_k + 3/32*k1 + 9/32*k2, t_k + 3/8*dt, par)
    k4 = dt * fun(x_k + 1932/2197*k1 + -7200/2197*k2 + 7296/2197*k3,t_k + 12/13*dt,  par)
    k5 = dt * fun(x_k + 439/216*k1 + -8*k2 + 3680/513*k3 + -845/4104*k4, t_k + 1*dt,  par)
    k6 = dt * fun(x_k + -8/27*k1 + 2*k2 + -3544/2565*k3 + 1859/4104*k4 + -11/40*k5, t_k + 1/2*dt,  par)

    x_k1_5th = x_k + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 + -9/50*k5 + 2/55*k6

    x_k1_4th = x_k + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 + -1/5*k5

    err = np.linalg.norm(x_k1_5th - x_k1_4th, ord = np.inf) # Error estimate
    blendTol = tol*(1 + np.linalg.norm(x_k1_5th, ord = np.inf))

    if err < blendTol:
        return x_k1_5th, dt
    else:      
        q = 0.8*pow((blendTol/err),(1/3)) # q to jest o ile zwiększymy/zmniejszymy time step, nie time step 
        q = min(q,4)
        dt_old = dt
        dt = min(q*dt, tf-t_k)
        #print('At {} Needed to change step size from {} to {}'.format(t_k, dt_old, dt))
        #print('Error estimate: {:.6f}'.format(err))
        return rk_45(fun, x_k, t_k, dt, par, tf, tol=1e-6) # W rekursji musi być return przy przywoływaniu drugi raz funkcji 


    

    #return x_k1_5th
