import numpy as np
from functions.rk4 import rk4_np 
import config2 as cfg
from math import sqrt
from julia import Main
Main.include('functions/rk4.jl')

class Sample():
    pass

class FixedPoint():
    '''
    For fixed point the tspan list is: 
    tspan[0] = t0
    tspan[1] = tf
    tspan[2] = dt
    '''

    def __init__(self, DynSys, integrator=None):
        self.init_cond = DynSys.init_cond
        self.par = DynSys.par
        self.RHS = DynSys.f
        self.sys_dim = len(self.init_cond) 
        self._marker_list = []
        self.parToSave = DynSys.parToSave
        self._t0 = DynSys.tspan[0] 
        self._tf = DynSys.tspan[1] 
        self._dt = DynSys.tspan[2]

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

    def integrate_fixed_step(self):
        t_sim = np.arange(self._t0, self._t0 + self._tf , self._dt)
        y_arr = np.zeros((len(t_sim),self.sys_dim))
        y = self.init_cond

        self._marker_list.append(self._t0)
        self._marker_list.extend([x for x in self.init_cond])
        self._marker_list.extend([self.par[x] for x in self.parToSave])

        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self._dt, self.par, self._tf)
            y = y_out

        self._y_arr = np.hstack((t_sim[:,np.newaxis],y_arr))
        self._marker_list.extend([int(x) for x in y_arr[-1,:]])
    
    def integrate_brute_force_julia(self):
        #timeJ = np.arange(self.t0, self.t0+self.tf ,self.dt)
        #y_arrJ = np.zeros((len(timeJ), self.sys_dim))
        
        self._marker_list.append(self.t0)
        self._marker_list.extend([x for x in self.init_cond])
        
        Main.par = self.par
        Main.init_cond = self.init_cond
        Main.t0 = self.t0
        Main.tf = self.t0+self.tf
        Main.dt = self.dt

        res = Main.eval("integration(par, init_cond, t0, tf, dt)")
        self._marker_list.extend(res[-1,1:])
        #y_arrJ =  res[0]
        #timeJ = np.array(res[1])

    @property
    def y_arr(self):
        print(self._marker_list)
        return self._y_arr
    
    @property
    def marker_list(self):
        return self._marker_list


class Periodic_NA(FixedPoint):
    '''
    For Periodic_NA the tspan list is: 
    tspan[0] = t0
    tspan[1] = periods_no
    tspan[2] = step
    tspan[3] = period_of_system
    '''
    def __init__(self, DynSys):
        super().__init__(DynSys)
        self.period_of_system = DynSys.tspan[3] 
        self.poincare_matrix = tuple([[] for _ in range(self.sys_dim+1)])
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

    def integrate_fixed_step(self):
        t_sim = np.arange(self._t0, self._t0 + self.tf_NA , self.dt_NA)
        y_arr = np.zeros((len(t_sim),self.sys_dim))
        y = self.init_cond

        self._marker_list.append(self._t0)
        self._marker_list.extend([x for x in self.init_cond])
        self._marker_list.extend([self.par[x] for x in self.parToSave])

        for i,t in enumerate(t_sim):
            y_arr[i,:] = y
            y_out = rk4_np(self.RHS, y, t, self.dt_NA, self.par, self.tf_NA)
            if cfg.error_check and y_out[0] > 1_000:
                print(f'Possible RuntimeWarning dla za mały czas całkowania dla {[self.par[x] for x in self.parToSave]}, sprawdź czy nie za mała wartość threshold, RuntimeWarning są wyłączone')
                break
            self.poincare(t, y, i)
            if self.check_periodicity():
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
            self.save_poincare()

    def poincare(self, t_curr, y_p,  i):
        if t_curr > cfg.poincare_t*self.tf_NA and i%self.step == 0:
            if len(self.poincare_matrix[0]) < cfg.Maxmaps: 
                for i in range(self.sys_dim):
                    self.poincare_matrix[i].append(y_p[i])
                else:
                    self.poincare_matrix[-1].append(t_curr)
    
    def check_periodicity(self):
        if len(self.poincare_matrix[0]) == cfg.Maxmaps:
            if self.find_period(self.poincare_matrix[0]) == self.find_period(self.poincare_matrix[0], ref = 1) and \
                self.find_period(self.poincare_matrix[1]) == self.find_period(self.poincare_matrix[1], ref = 1) and \
                self.find_period(self.poincare_matrix[0]) != -1:
                
                self.period = max(self.find_period(self.poincare_matrix[1]), self.find_period(self.poincare_matrix[0]))
                return True #if self.period !=-1 else False
            else: 
                for i in range(self.sys_dim+1): # +1 bo jeszcze czas jest na końcu w poincare_matrix
                    try:
                        self.poincare_matrix[i].pop(0)
                    except:
                        pass
                return False
        else:
            return False
        
    def find_period(self, poincareMap,  ref=0):
        suma = 0 
        a_old = 0 

        for i in range(ref+1,len(poincareMap)):
            distance = sqrt( (poincareMap[ref] - poincareMap[i])*(poincareMap[ref] - poincareMap[i]) ) 
            if distance < cfg.poincare_distance :
                a = i-ref 
                ref = i
                #print(f'suma: {suma} ')
                #print(f'distance: {distance}, cfg.distance: {cfg.poincare_distance}, {distance < cfg.poincare_distance}')
                #print(f'a: {a}, a_old: {a_old}')
                if a == a_old: 
                    suma +=1
                else: 
                    suma = 0
                a_old = a

            if suma == 2 : 
                return a
        return -1 
    
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

    def save_poincare(self):
        with open('results\\poincare_map.txt', 'w') as p_file:
            for i in range(len(self.poincare_matrix)): 
                for j in range(len(self.poincare_matrix[i])):
                    p_file.writelines(str(self.poincare_matrix[i][j]) +' ')
                p_file.writelines('\n')


    
