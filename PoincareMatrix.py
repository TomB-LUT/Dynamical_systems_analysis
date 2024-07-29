
from math import sqrt
import config2 as cfg
import numpy as np


class LengthError(Exception):
    pass

class PoincareMatrix:

    def __init__(self, dimension, max_maps=20):
        self.dimension = dimension
        self.max_maps = max_maps
        self.p_mat = [[] for _ in range(self.dimension+1)]

    def __len__(self):
        return len(self.p_mat[0])

    def __str__(self):
        return f'{[ print(["{:10.5f}".format(e) for e in row]) for row in self.p_mat ]}'
    
    @staticmethod
    def flag(sample, state_vars, state_vars_old, i):
        if sample.period_of_system == None:
            return state_vars[0]*state_vars_old[0] < 0 and state_vars[1] > 0
        else: 
            return i%sample.step == 0

    def push(self,vector):
        if len(vector) != len(self.p_mat):
            raise LengthError('State vector must be the same length as poincare matrix rows number')
        
        if not self.is_complete:
            for element, row in zip(vector, self.p_mat):
                row.append(element)

    def remove_first(self):
        for row in self.p_mat:
            row.pop(0)
    
    @property
    def is_complete(self):
        return len(self) == self.max_maps

    def periodicity_found(self, sample):
        if not self.is_complete:
            return False

        counter = 0 
        ref = 0 
        period = 1
        np_p_mat = np.array(self.p_mat[1:])
        for _ in range(cfg.max_period):
            if all( np.absolute(np_p_mat[:,ref] - np_p_mat[:,ref+period]) < cfg.poincare_distance ):
                counter +=1
                ref = period
            else: 
                period +=1
                counter = 0 
                
            if counter == cfg.poincare_iter:
                sample.period = period
                if sample.period_of_system == None:
                    sample.period_of_system = np.average(np.diff(self.p_mat[0]))
                return True
            
        else:
            self.remove_first()
            return False
        
    def save_to_txt(self):
        with open('results\\poincare_matrix.txt', 'w') as p_file:
            for i in range(len(self.p_mat)): 
                for j in range(len(self.p_mat[i])):
                    p_file.writelines(str(self.p_mat[i][j]) +' ')
                p_file.writelines('\n')