
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

    def find_period(self, row, ref=0):
        counter = 0 
        position_old = 0 

        for i in range(ref+1,len(row)):
            distance = sqrt( (row[ref] - row[i])*(row[ref] - row[i]) ) 
            if distance < cfg.poincare_distance :
                position = i-ref 
                ref = i
                if position == position_old: 
                    counter +=1
                else: 
                    counter = 0
                position_old = position
            
            if counter == 2 : 
                #print(row)
                #print(position)
                return position
        
        #print(row)
        #print(position)
        return -1 
    
    def periodicity_found_2(self):
        if not self.is_complete:
            return False

        temp_check = [self.find_period(row) == self.find_period(row, ref=1) for row in self.p_map[1:] ]
        #print(temp_check)
        
        if all(temp_check):
            
            period = max([self.find_period(row) for row in self.p_map[1:]])
            
            if period != -1:
                return period  
            
            # Trick do przetestowania
            #if max([self.find_period[row] for row in self.p_map[1:]])
        else:
            self.remove_first()
            return False

    def periodicity_found(self, sample):
        #print(sample.period.value)
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

class Marker:
    
    def __init__(self, name, value, marker_list):
        self._name = name
        self._value = value
        marker_list.append_marker(self)

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self,value):
        self._value = value

class MarkerList:
    def __init__(self):
        self.marker_list = []
    
    def append_marker(self,marker_object):
        self.marker_list.append(marker_object)

    def return_values(self):
        return [marker._value for marker in self.marker_list]
    
    def return_names(self):
        return [marker._name for marker in self.marker_list]
    
#Parameter też jako klas, ma name i value, klasa nazywa się par, w równaniach wtedy par.m1, name jako property które zwraca value


            









'''
pm_2 = PoincareMatrix(3)
for _ in range(10):
    pm_2.push([1,9,1,7])
    pm_2.push([1,2,4,6])
    pm_2.push([1,2,2,9])
    pm_2.push([1,2,4,6])
print(pm_2)

#[print(row) for row in pm_2.p_map]
#print(pm_2.check_periodicity())
#print(pm_2.periodicity_found_2())


'''