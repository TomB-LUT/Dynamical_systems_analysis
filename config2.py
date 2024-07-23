from os import cpu_count

from functions.rk4 import rk4_np, rk_45, rk4_sd
# Zrób tak: każde sample to jest obiekt. Metody to: draw IC(), integrate()


no_of_sampels = 1
report_samples = 5_000
Maxmaps = 20


poincare_distance = 0.0005
workers = cpu_count()-12
poincare_t = 0.25
do_marker = True
sys_type = 1 # 2 - autonomiczny, 1 - nieautonomiczny 
dim = 4
save_poincare = True
error_check = False
integration_algorithm = rk4_np
poincare_iter = 2