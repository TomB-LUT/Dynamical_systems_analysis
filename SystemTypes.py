import time
import numpy as np
from functions.rk4 import rk4_np 
import config2 as cfg
from math import sqrt
from special_containers import PoincareMatrix
from julia import Main
Main.include('functions/rk4.jl')

class LimitSet():

    def __init__(self, DynSys):
        self.init_cond = DynSys.init_cond
        self.par = DynSys.par
        self.RHS = DynSys.f
        self.sys_dim = len(self.init_cond) 
        self._marker_list = []
        self.parToSave = DynSys.par_to_save()
        self._t0 = DynSys.tspan[0] 
        self._tf = DynSys.tspan[1] 
        self._dt = DynSys.tspan[2]
        self._marker_list.append(self._t0)
        self._marker_list.extend([x for x in self.init_cond])
        if self.parToSave != None:
            self._marker_list.extend([self.par[x] for x in self.parToSave])
    
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
        t_sim = np.arange(self._t0, self._t0 + self._tf , self._dt)
        y_arr = np.zeros((len(t_sim),self.sys_dim))
        y = self.init_cond

        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self._dt, self.par, self._tf)
            y = y_out

        self._y_arr = np.hstack((t_sim[:,np.newaxis],y_arr))
        self._marker_list.extend([round(x) for x in y_arr[-1,:]])
    
    def integrate_julia(self):
        
        Main.par = self.par
        Main.init_cond = self.init_cond
        Main.t0 = self.t0
        Main.tf = self.t0+self.tf
        Main.dt = self.dt

        self._y_arr = Main.eval("integration(par, init_cond, t0, tf, dt)")
        self._marker_list.extend([round(x) for x in self._y_arr[-1,1:]])

    @property
    def y_arr(self):
        print(self._marker_list)
        return self._y_arr


class Periodic_NA(LimitSet):
    '''
    For Periodic_NA the tspan list is: 
    tspan[0] = t0
    tspan[1] = number of periods to calculate
    tspan[2] = time step
    tspan[3] = period of system
    '''
    def __init__(self, DynSys):
        super().__init__(DynSys)
        self.period_of_system = DynSys.tspan[3] 
        self.poincare_matrix = PoincareMatrix(self.sys_dim)
        self.period = -1
        self.well = -1
        self.attractor = 0 #To poprostu może stwórz klase macierz markerów i bedzie dodawana. A może instacje markera zrobić jak obserwatorów
        for i in range(self.sys_dim):
            setattr(self, str(f'max_y{i}'), -1000)
            setattr(self, str(f'min_y{i}'), 1000)

    @property
    def tf_NA(self):
        return self._tf*self.period_of_system    
    
    @property
    def dt_NA(self):
        return self.period_of_system/self._dt
    
    @property
    def step(self):
        return self._dt

    def integrate_fixed_step(self, start_time = None, start_y_arr = None):

        if start_time == None:
            t_sim = np.arange(self.t0, self.t0+self.tf_NA ,self.dt_NA)
        else: 
            t_sim = np.arange(start_time, start_time+round(self.tf_NA/20) ,self.dt_NA)

        if not isinstance(start_y_arr,np.ndarray):
            y = self.init_cond
        else: 
            y = start_y_arr

        y_arr = np.zeros((len(t_sim),self.sys_dim))
        
        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self.dt_NA, self.par, self.tf_NA)
            if cfg.error_check and y_out[0] > 1_000:
                print(f'Possible RuntimeWarning dla za mały czas całkowania dla {[self.par[x] for x in self.parToSave]}, sprawdź czy nie za mała wartość threshold, RuntimeWarning są wyłączone')
                break

            if t > cfg.poincare_t*self.tf_NA and i%self.step == 0:
                self.poincare_matrix.push([t,*y])

            if self.poincare_matrix.periodicity_found(self):
                self.marker(rk4_np, y_out, t+self.dt_NA, self.dt_NA)
                t_sim = t_sim[0:i]
                y_arr = y_arr[0:i]
                break

            y = y_out
            
        self._y_arr = np.hstack((t_sim[:,np.newaxis],y_arr))
        
        self._marker_list.append(self.period)
        self._marker_list.append(self.well)
        for i in range(self.sys_dim):
            self._marker_list.append( getattr(self, str(f'max_y{i}')) )
            self._marker_list.append( getattr(self, str(f'min_y{i}')) )
        self._marker_list.append(self.attractor)
        
        if cfg.no_of_sampels == 1:
            self.poincare_matrix.save_to_txt()

    def integrate_julia(self):
        Main.par = self.par
        Main.init_cond = self.init_cond
        Main.t0 = self.t0
        Main.tf = self.t0+self.tf_NA
        Main.dt = self.dt_NA

        res = Main.eval("integration(par, init_cond, t0, tf, dt)")
        self.integrate_fixed_step(res[-1,0], res[-1,1:])

    @property
    def y_arr(self):
        print(self._marker_list)
        return self._y_arr

    def marker(self, integration_algorithm, yM, t, dtM):

        time_M = np.arange(t, t+(self.period_of_system*self.period), dtM)
        y_ref = yM[0]
        yM_arr = np.zeros((len(time_M),self.sys_dim))

        for i, tM in enumerate(time_M):
            yM_arr[i,:] = yM
            y_outM = integration_algorithm(self.RHS, yM, tM, dtM, self.par, self.tf_NA)
            yM = y_outM 

            self.attractor += y_outM[0]

            cond = y_ref*y_outM[0] < 0 or self.well == 2
            self.well = 2 if cond else 1 
            
            #inst_energy = y_outM[2]*y_outM[2]* self.par[5] * dtM + y_outM[2]*y_outM[2]* self.par[5] * dtM # I^2 * R * delta_t
            #self.total_energy += inst_energy

            for i in range(self.sys_dim):
                if y_outM[i] > getattr(self, str(f'max_y{i}')):
                    setattr(self, str(f'max_y{i}'), y_outM[i])
                if y_outM[i] < getattr(self, str(f'min_y{i}')):
                    setattr(self, str(f'min_y{i}'), y_outM[i])
            
        yM_arr = np.hstack((time_M[:,np.newaxis],yM_arr))
        if cfg.no_of_sampels == 1:
            np.savetxt('results\\one_period.txt', yM_arr, delimiter=' ', fmt='%.6f')




    

