import numpy as np 
import matplotlib.pyplot as plt
import sys
import matplotlib.pyplot as plt 
import matplotlib.patches as ptch
import matplotlib.colorbar as cbar
import pylab as pl
import concurrent.futures
import os
from operator import le, ge


def print_wrap(func):
    def func_with_wrap(*args, **kwargs):
        print('-'*10)
        results = func(*args,**kwargs)
        print('-'*10)
        return results
    return func_with_wrap

def concurrent_wrap(func):
    def func_to_multiP(*args, **kwargs):
        #Tutaj kod do concurent 
        results = 0 
        return results
    return func_to_multiP

#Zrób wrap na concurent futures?



class BoA:
    
    def __init__(self, path):
        if isinstance(path,str):
            self.data = np.loadtxt(path,skiprows=1)
        else:
            self.data = path 
        self.chaosCol = 11

    @print_wrap
    def find_attractors(self,atrCol):
        self.attractors = np.unique(self.data[:,atrCol])
        self.basin_stab = []
        print('I found following attractors: {}'.format(list(self.attractors)))
        for a in self.attractors:
            self.basin_stab.append( np.shape(self.data[self.data[:,atrCol] == a,:])[0] / len(self.data[:,1]) )
        self.atrCol = atrCol
    
    @print_wrap
    def eliminate_minorities(self): #Elimnate minorities może być dekoratorem albo wrapperem ?????
        print('I deleted following attractors: ')
        for a in self.attractors:
            if np.shape(self.data[self.data[:,self.atrCol] == a,:])[0] / len(self.data[:,1]) < 0.01:
                rows_to_del = np.where(self.data[:,self.atrCol] == a)
                self.data = np.delete(self.data, rows_to_del[0], axis= 0)
                print(a)
            else:
                print(None)
        self.attractors = np.unique(self.data[:,self.atrCol])
    
    def set_chaos_col(self, chaosCol):
        self.chaosCol = chaosCol


    @print_wrap
    def show_basin_info(self):
        try: 
            chaosy = self.data[self.data[:,self.chaosCol]== -1]
            chaosyLen = np.shape(chaosy)[0]
            chaosyPercentage = chaosyLen / len(self.data[:,1])
            print('In this basin there are {0} chaotic solutions which is {1:.2} % of all samples'.format(chaosyLen, chaosyPercentage*100 ))
        except IndexError as e: 
            print('Nie mam aż tylu markerów żeby sprawdzić chaos: %s' % e)
        except:
            print(sys.stderr)

        print('In this basin there are following attractors with its basin stability: ')
        for a, bs in zip(self.attractors, self.basin_stab):
            print('Attractor {0} - {1:.1f}'.format(a,bs*100))

    @print_wrap
    def delete_chaos(self):
        chaos_rows = np.where(self.data[:,self.chaosCol] == -1)
        if len(chaos_rows[0]) > 0:
            print('I found {} chaotic solutions and deleted it'.format(len(chaos_rows[0])))
            self.data = np.delete(self.data, chaos_rows[0], axis= 0)
        else:
            print('No chaotic solutions found')
        self.find_attractors(self.atrCol)


    
    def print_attractors(self):
        print('Now, we have following attractors: {}'.format(self.attractors))    


    def plot_BoA(self,cor1,cor2):
        fig, ax = plt.subplots()
        for a in np.flip(self.attractors):
            temp_data = self.data[self.data[:,self.atrCol]==a]
            ax.scatter(temp_data[:,cor1],temp_data[:,cor2],s=2, label=a)
        plt.legend()
        plt.show()
    
    def plot_BoA_separate(self,cor1,cor2):
        for a in np.flip(self.attractors):
            fig, ax = plt.subplots()
            temp_data = self.data[self.data[:,self.atrCol]==a]
            ax.scatter(temp_data[:,cor1],temp_data[:,cor2],s=2, label=a)
            plt.legend()
            plt.show()

    def combine_similar_attractors(self):
        self.data[:,self.atrCol] = np.round(self.data[:,self.atrCol],2)
        print(self.data[:,self.atrCol])
        self.find_attractors(self.atrCol)

    @staticmethod
    def join_marker(paths,f_name):
        res = np.loadtxt(paths[0],skiprows=1)
        paths = paths[1:]

        for p in paths:
            a = np.loadtxt(p,skiprows=1)
            res = np.append(res,a,axis=0)
        np.savetxt(f_name, res, fmt='%.6f')

    @staticmethod
    def slice_data(path, f_name, func, col, value):
        data_to_slice = np.loadtxt(path)
        #print(np.shape(data_to_slice))
        #new = data_to_slice[func(data_to_slice[:, col] , value)]
        #print(np.shape(new))

        np.savetxt(f_name, data_to_slice[func(data_to_slice[:, col] , value)], fmt='%.6f')

class ParametricData1D:

    def __init__(self, path, par, marker):
        self.data = np.loadtxt(path, skiprows=1)
        self._par = par
        self._marker = marker

    def plot(self):
        for i in np.unique(self.data[:,self._par]):
            #temp_data = self.data[self.data[:,self.atrCol]==a]
            temp = self.data[self.data[:,self._par]== i]
            a = BoA(temp)
            a.find_attractors(self._marker)
            a.plot_BoA(1,2)
            one = temp[temp[:,self._marker] == 1]

            print(len(one)/len(temp))


class ParametricData2D:

    def __init__(self, path, par1, par2, color_indic, userXmax = None):
        self.data = np.loadtxt(path,skiprows=1)
        self._par1 = par1
        self._par2 = par2
        self._color_indic = color_indic
        #self.xMax = userXmax
        if (userXmax):
            self.xMax = userXmax
        else:
            self.xMax = max(self.data[:,self._par1])
        self.xMin = min(self.data[:,self._par1]) 
        self.yMax = max(self.data[:,self._par2])
        self.yMin = min(self.data[:,self._par2])

#    @property
#    def xMax(self):
#        return self._xMax
    
#    @xMax.setter
#    def xMax(self, value):
#        if not value:
#            self._xMax = max(self.data[:,self._par1])
#        else:
#            self._xMax = value

    def show_data_info(self, chaosCol):
        try: 
            chaosy = self.data[self.data[:,chaosCol]== -1]
            chaosyLen = np.shape(chaosy)[0]
            chaosyPercentage = chaosyLen / len(self.data[:,1])
            print('In this basin there are {0} chaotic solutions which is {1:.2f} percentage of all samples'.format(chaosyLen, chaosyPercentage*100 ))
            print('In this basin there are {0} samples'.format(np.shape(self.data)[0]))
        except IndexError as e: 
            print('Nie mam aż tylu markerów żeby sprawdzić chaos: %s' % e)
        except:
            print(sys.stderr)

        #print('In this basin there are following attractors with its basin stability: ')
        #for a, bs in zip(self.attractors, self.basin_stab):
        #    print('Attractor {0} - {1:.1f}'.format(a,bs*100))

    def make_mesh_divByNum(self, DivNumX, DivNumY):

        yRange = np.linspace(self.yMin,self.yMax,DivNumY,endpoint=False)
        xRange = np.linspace(self.xMin,self.xMax,DivNumX, endpoint=False)

        self.yWall = (self.yMax-self.yMin)/DivNumY
        self.xWall = (self.xMax-self.xMin)/DivNumX

        self.rectList = [(x,y) for x in xRange for y in yRange]
        print('In this basin there are {} boxes'.format(len(self.rectList)))

    def make_mesh_divByWalSize(self, xWallSize, yWallSize):

        self.xWall = xWallSize
        self.yWall = yWallSize

        yRange = np.arange(self.yMin,self.yMax,self.yWall)
        xRange = np.arange(self.xMin,self.xMax,self.xWall)

        self.rectList = [(x,y) for x in xRange for y in yRange]

    def make_mesh(self):

        xMax = max(self.data[:,self._par1])
        xMin = min(self.data[:,self._par1]) 
        yMax = max(self.data[:,self._par2])
        yMin = min(self.data[:,self._par2])

        factorY = 3
        factorX = 6
        noY = int( (yMax-yMin)*factorY)
        noX = int( (xMax-xMin)*factorX)

        yRange = np.linspace(yMin,yMax,noY+1)
        xRange = np.linspace(xMin,xMax,noX+1)

        self.yWall = (yMax-yMin)/noY
        self.xWall = (xMax-xMin)/noX


        self.xMax = xMax
        self.xMin = xMin
        self.yMax = yMax
        self.yMin = yMin

        self.rectList = [(x,y) for x in xRange[:-1] for y in yRange[:-1]]

    def calculate_prob_concurent(self):
        with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()-8) as executor:
            results = [executor.submit(self.calculate_single_prob, rect) for rect in self.rectList]

        self.probabilities = [x.result() for x in results]

    def calculate_prob(self):
        self.probabilities = []
        for d in self.rectList: 
            self.probabilities.append(self.calculate_single_prob(d))
        
    def calculate_single_prob(self, rect):
        xLB, yLB = rect
        xUpperBound = xLB+self.xWall
        yUpperBound = yLB+self.yWall

        xPointsInRect = [x for x in np.where( np.logical_and(self.data[:,self._par1]>xLB , self.data[:,self._par1]<=xUpperBound))  ]
        yPointsInRect = [y for y in np.where( np.logical_and(self.data[:,self._par2]>yLB , self.data[:,self._par2]<=yUpperBound))  ]
        
        rowsRect = np.intersect1d(xPointsInRect, yPointsInRect) 
        
        currBasin = self.data[rowsRect,:]
        currBasinBig = currBasin[currBasin[:,self._color_indic]==2]
        #print(f'no of all: {np.shape(currBasin)[0]}' )
        #print(f'no of double well: {np.shape(currBasinBig)[0]}' )
        #print(f'Probability of double well: {np.shape(rowsCurrBasinBig)[1]/np.shape(currBasin)[0]}')
        try:
            probability = np.shape(currBasinBig)[0]/np.shape(currBasin)[0]
            sampels_in = np.shape(currBasin)[0]
        except ZeroDivisionError as e :
            print('Dzielenie przez ZERO przy liczeniu probability!!!')
            probability = 0 
            sampels_in = np.shape(currBasin)[0]
            print(f'no of all: {sampels_in}' )


        return probability, sampels_in

    def compare_prob(self, other):
        self.make_mesh_divByNum(40,40)
        self.calculate_prob()
        other.make_mesh_divByNum(40,40)
        other.calculate_prob()
        p1, n1 =(zip(*self.probabilities))
        p2, n2 =(zip(*other.probabilities))
        diff = np.array(p1) - np.array(p2)
        print(f'Max prob diff: {max(diff)}')
        self.probabilities = tuple(zip(abs(diff),n1))
        self.plot_res(titleName='Probability comparison', norm_range=(0,1))
        return None 

    
    def plot_raw_points(self):
        fig = plt.figure()
        ax = plt.axes()
        ax.scatter(self.data[:,self._par1], self.data[:,self._par2])
        plt.show()

    def plot_time_characteristics(self, tf_col):
        fig = plt.figure()
        ax = plt.axes()
        ax.scatter( np.arange(1,len(self.data[:,tf_col])+1 ), self.data[:,tf_col] )
        plt.show()

    
    def plot_res(self, titleName='Some title',
                  ylabel = '', xlabel = '', plot_num = False, userXmin = None, 
                  userYmin = None, userXmax = None, userYmax= None, norm_range = None):
        
        if (userXmin):
            self.xMin = userXmin
        if (userYmin):
            self.yMin = userYmin
        if (userXmax):
            self.xMax = userXmax
        if (userYmax):
            self.yMax = userYmax

        probabilities, no_of_all = (zip(*self.probabilities))

        normal = pl.Normalize(min(probabilities), max(probabilities)) if norm_range == None else pl.Normalize(*norm_range)
        colors = pl.cm.jet(normal(probabilities))

        fig, ax = plt.subplots()
        for rect,color, amount in zip(self.rectList,colors,no_of_all):
            rectToPlot = ptch.Rectangle(rect, self.xWall, self.yWall, facecolor = color, edgecolor = 'black',linewidth = 0.1)
            ax.add_patch(rectToPlot)
            if plot_num == True:
                ax.text(rect[0], rect[1],str(amount), fontsize='x-small')
            
        #ax.set_xlim(np.floor(self.xMin),np.ceil(self.xMax))
        #ax.set_xlim(np.floor(self.xMin),0.1)
        ax.set_xlim(self.xMin,self.xMax)
        #ax.set_ylim(np.floor(self.yMin),np.ceil(self.yMax))
        ax.set_ylim(self.yMin,self.yMax)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(titleName)
        cax, _ = cbar.make_axes(ax) 
        cb2 = cbar.ColorbarBase(cax, cmap=pl.cm.jet,norm=normal) 
        cb2.set_label('Corss-well probability')
        plt.show()



        
        

if __name__ == "__main__":
 
    a = ParametricData1D(r'results\marker.txt', -3, -2)
    a.plot()