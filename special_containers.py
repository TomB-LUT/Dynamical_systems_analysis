
from math import sqrt
import config2 as cfg
import numpy as np


class LengthError(Exception):
    pass


class PoincareMap:


    def __init__(self, dimension, max_maps=20):
        self.dimension = dimension
        self.max_maps = max_maps
        self.p_map = [[] for _ in range(self.dimension+1)]

    #def __repr__(self):
    #    [print(row) for row in self.p_map]
    #    return None

    def push(self,vector):
        if len(vector) != len(self.p_map):
            raise LengthError('State vector must be the same length as poincare map')
        
        if not self.is_complete:
            for element, row in zip(vector, self.p_map):
                row.append(element)

    def remove_first(self):
        for row in self.p_map:
            row.pop(0)
    
    @property
    def is_complete(self):
        return len(self.p_map) == self.max_maps

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
    
    def check_periodicity_2(self):
        #if not self.is_complete:
        #    return False

        temp_check = [self.find_period(row) == self.find_period(row, ref=1) for row in self.p_map[1:] ]
        print(temp_check)
        
        if all(temp_check):
            
            period = max([self.find_period(row) for row in self.p_map[1:]])
            
            if period != -1:
                return period  
            
            # Trick do przetestowania
            #if max([self.find_period[row] for row in self.p_map[1:]])
        else:
            self.remove_first()
            return False

    def check_periodicity(self):
        #if not self.is_complete:
        #    return False


        counter = 0 
        ref = 0 
        period = 1
        max_period = 5
        np_p_map = np.array(self.p_map)
        for _ in range(max_period):
            if all( (np_p_map[:,ref] - np_p_map[:,ref+period]) < cfg.poincare_distance ):
                counter +=1
                ref = period
            else: 
                period +=1
                counter = 0 
                
            if counter == cfg.poincare_iter:
                return period
            
        else:
            self.remove_first()
            return False






pm_2 = PoincareMap(3)
for _ in range(10):
    pm_2.push([1,9,1,7])
    pm_2.push([1,2,4,6])
    pm_2.push([1,2,2,9])
    pm_2.push([1,2,4,6])


#[print(row) for row in pm_2.p_map]
#print(pm_2.check_periodicity())
print(pm_2.check_periodicity_2())


