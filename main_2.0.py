from Li_2023_nonDim import LiNonDim2023
from DuffingFixedPoint import DuffingFixedPoint
from SystemTypes import FixedPoint, Periodic_NA
import matplotlib.pyplot as plt
import config2 as cfg
import concurrent.futures 
import numpy as np
import time
import os

completed = 0 
def save_marker(data, my_class=None):    
    file_handler = open('results\\marker.txt', 'w')
    for i in range(len(data)):
        for j in range(len(data[i])):
            file_handler.writelines(str(data[i][j]) + ' ')
        file_handler.writelines('\n')
    file_handler.close()


def save_traj(data):
    np.savetxt('results\\trajectory.txt',data, delimiter=' ', fmt='%.6f')


def progress_indicator(result):
    global completed
    completed +=1
    #print(result.result())
    if (cfg.no_of_sampels-completed) % cfg.report_samples == 0:
        print(f'Remaining tasks: {cfg.no_of_sampels-completed}')

def execute_obj(system_type, system):
    sample = system_type(system())
    sample.integrate_fixed_step()
    if cfg.no_of_sampels == 1:
        return sample.y_arr
    else: 
        return sample.marker_list


def main_concurent(system_type, system):
    print(f'Remaining tasks: {cfg.no_of_sampels}')
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.workers) as executor:
        results = [executor.submit(execute_obj, system_type, system) for _ in range(cfg.no_of_sampels)]
        for result in results:
            result.add_done_callback(progress_indicator)

    if cfg.no_of_sampels == 1: 
        futures = [x.result() for x in concurrent.futures.as_completed(results)] 
        start1 = time.time()
        #print(futures[0])
        save_traj(futures[0])
        print(f'save time: {time.time()-start1}')
    else:
        start1 = time.time()
        save_marker([x.result() for x in concurrent.futures.as_completed(results)])
        print(f'save time: {time.time()-start1}')

def main_single( system_type, system):
    try:
        sample = system_type(system())
        sample.integrate_fixed_step()
    except TypeError as e:
        print('Not all required methods implemented!!!!. Details: '.upper())
        print(e)

    try:
        res = sample.y_arr
    except AttributeError as e:
        print(e)
        quit()

    save_traj(res)

    figure = plt.figure()
    ax = plt.axes()
    ax.plot(res[:,0], res[:,1])
    plt.show()

if __name__ == "__main__":

    
    fast_del =  not os.path.isfile('results\\trajectory.txt') or os.remove('results\\trajectory.txt')
    fast_del =  not os.path.isfile('results\\poincare_map.txt') or os.remove('results\\poincare_map.txt')
    fast_del =  not os.path.isfile('results\\one_period.txt') or os.remove('results\\one_period.txt')

        

    main_concurent( system_type = Periodic_NA, system = LiNonDim2023 )
    #main_single( system_type = Periodic_NA, system = LiNonDim2023 )
    #sample = Periodic_NA(LiNonDim2023())
    #print(sample.tf)
    #print(sample.__dict__)
    #sample.integrate_fixed_step()
    #print(sample.__dict__)
    
