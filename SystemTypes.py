import time
import numpy as np
from functions.rk4 import rk4_np 
import config2 as cfg
from math import sqrt
from PoincareMatrix import PoincareMatrix
from julia import Main
Main.include('functions/rk4.jl')

class LimitSet():

    def __init__(self, DynSys, use_julia ):
        self.init_cond = DynSys.init_cond
        self.par = DynSys.par
        self.RHS = DynSys.f
        self.sys_dim = len(self.init_cond) 
        self._marker_list = {}
        self.parToSave = DynSys.par_to_save()
        self._t0 = DynSys.tspan[0] 
        self._tf = DynSys.tspan[1] 
        self._dt = DynSys.tspan[2]
        self._marker_list['t0'] = self._t0
        for i, ic in enumerate(self.init_cond):
            self._marker_list[f'y_{i+1}0'] = ic
        if self.parToSave != None:
            for i in self.parToSave:
                self._marker_list[f'par_{i}'] = self.par[i]
        self.use_julia = use_julia

    @property
    def t0(self):
        return self._t0
        
    @t0.setter
    def t0(self, value):
        self._t0 = value

    @property
    def tf(self):
        return self._tf
        
    @tf.setter
    def tf(self, value):
        self._tf = value

    @property
    def dt(self):
        return self._dt
        
    @dt.setter
    def dt(self, value):
        self._dt = value
    
    @property
    def marker_list(self):
        return self._marker_list


class FixedPoint(LimitSet):
    '''
    For fixed point the tspan list is: 
    tspan[0] = t0
    tspan[1] = tf
    tspan[2] = dt
    '''
    
    def integrate_fixed_step(self):

        if self.use_julia:
            self.integrate_julia()
            return
        
        t_sim = np.arange(self._t0, self._t0 + self._tf , self._dt)
        y_arr = np.zeros((len(t_sim),self.sys_dim))
        y = self.init_cond

        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self._dt, self.par, self._tf)
            y = y_out

        self._y_arr = np.hstack((t_sim[:,np.newaxis],y_arr))
        for i, ic in enumerate(self._y_arr[-1,1:]):
            self._marker_list[f'y_{i+1}_last'] = round(ic)
    
    def integrate_julia(self):
        
        Main.par = self.par
        Main.init_cond = self.init_cond
        Main.t0 = self.t0
        Main.tf = self.t0+self.tf
        Main.dt = self.dt

        self._y_arr = Main.eval("integration(par, init_cond, t0, tf, dt)")
        for i, ic in enumerate(self._y_arr[-1,1:]):
            self._marker_list[f'y_{i+1}_last'] = round(ic)

    @property
    def y_arr(self):
        print(self._marker_list)
        return self._y_arr
    
class PeriodicAutonomus(LimitSet):

    def __init__(self, DynSys, use_julia):
        super().__init__(DynSys, use_julia)
        self.poincare_matrix = PoincareMatrix(self.sys_dim)
        self.period = -1
        self.period_of_system = None
        self._marker_list['empty_marker'] = 0
        self._marker_list['period'] = self.period
        self._marker_list['well'] = -1
        self.attractor = 0 
        for i in range(self.sys_dim):
            self._marker_list[f'max_y{i+1}'] = -1000
            self._marker_list[f'min_y{i+1}'] = 1000
        self._marker_list['attractor'] = 0

    def integrate_fixed_step(self):
        
        if self.use_julia:
            self.integrate_julia()
        
        t_sim = np.arange(self.t0, self.t0+self.tf ,self.dt)
        y_arr = np.zeros((len(t_sim),self.sys_dim))
        y = self.init_cond
        
        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self.dt, self.par, self.tf)
            if cfg.error_check and y_out[0] > 1_000:
                if self.parToSave != None:
                    print(f'Possible RuntimeWarning for {[self.par[x] for x in self.parToSave]}, check time step size')
                else:
                    print(f'Possible RuntimeWarning for {self.par}')

                break

            if t > cfg.poincare_t*self.tf and PoincareMatrix.flag(self, y_out, y, i):
                if self.period_of_system == None:
                    dt_p = self.dt*(y[0]/(y[0]-y_out[0]))
                    y_out_p = rk4_np(self.RHS, y, t, dt_p, self.par, self.tf)
                    self.poincare_matrix.push([t+self.dt/2,*y_out_p])
                else: 
                    self.poincare_matrix.push([t,*y])

            if self.poincare_matrix.periodicity_found(self):
                if cfg.no_of_sampels == 1:
                    self.poincare_matrix.save_to_txt()
                self.marker(rk4_np, y_out, t+self.dt, self.dt)
                t_sim = t_sim[0:i]
                y_arr = y_arr[0:i]
                self._marker_list['period'] = self.period
                break

            y = y_out
            
        self._y_arr = np.hstack((t_sim[:,np.newaxis],y_arr))
        
        if cfg.no_of_sampels == 1:
            self.poincare_matrix.save_to_txt()

    def integrate_julia(self):
        Main.par = self.par
        Main.init_cond = self.init_cond
        Main.t0 = self.t0
        Main.tf = self.t0+self.tf
        Main.dt = self.dt
        res = Main.eval("integration(par, init_cond, t0, tf, dt)")
        self.t0, *self.init_cond = res[-1,:]     

    @property
    def y_arr(self):
        print(self._marker_list)
        return self._y_arr
    
    def marker(self, integration_algorithm, yM, t, dtM):

        time_M = np.arange(t, t+(self.period_of_system*self.period), dtM)
        yM_arr = np.zeros((len(time_M),self.sys_dim))

        for i, tM in enumerate(time_M):
            yM_arr[i,:] = yM
            y_outM = integration_algorithm(self.RHS, yM, tM, dtM, self.par, self.tf)
            yM = y_outM 

            self._marker_list['attractor'] += y_outM[0]
            
            for i in range(self.sys_dim):
                if y_outM[i] > self._marker_list[f'max_y{i+1}']:
                    self._marker_list[f'max_y{i+1}'] = y_outM[i]
                if y_outM[i] < self._marker_list[f'min_y{i+1}']:
                    self._marker_list[f'min_y{i+1}'] = y_outM[i]

        yM_arr = np.hstack((time_M[:,np.newaxis],yM_arr))
        if cfg.no_of_sampels == 1:
            np.savetxt('results\\one_period.txt', yM_arr, delimiter=' ', fmt='%.6f')

class PeriodicNonAutonomus(PeriodicAutonomus):

    def __init__(self, DynSys, use_julia):
        super().__init__(DynSys, use_julia)
        self.period_of_system = DynSys.tspan[3] 
        self._marker_list['empty_marker'] = 0
        self._marker_list['well'] = -1
        self.step = self._dt
        self._dt = self.period_of_system/self._dt
        self._tf = self._tf*self.period_of_system 
