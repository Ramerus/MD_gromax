import gromax_lib as gl
from MDAnalysis.coordinates.XTC import XTCReader
import time
import numpy as np
from numpy import save, zeros, array
from os import listdir, getcwd, chdir, mkdir
from os.path import isfile, join, isdir
import matplotlib.pyplot as plt
import math as m
#import drop_lib as dl
#from  draw_lib import graph_ift, surface_show_saddle, surface_show_drop
#'''
dx = 0.15
file_type = 'xtc'
num_files = 100
fig_name = 'test'
#constants for all experiments(sepends on mineral)
solid = 4.0 #hight of mineral
#solid = 0
shift = 10 #shift for liquid componets
#type of experiment
# 0 - saddle
# 1 - drop
texp = 0
#go to directory with all experiments
exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\drop_lib\\quartz_aro'
chdir(exp_path)
#lets find all experiments
path = './'
dirs = listdir(path)
#lets go throught all experiments
#print(dirs)
angles = []
for dir in dirs:
    files_xtc = gl.find_files(exp_path, dir, file_type)
    files_gro = gl.find_files(exp_path, dir, 'gro')
    system_data = gl.mol_desc(files_gro[0])
    exp_data = XTCReader(files_xtc[-1])
    exp = gl.SAC(exp_data, num_files, system_data, shift, solid)
    exp_result, exp_delta = exp.analyse()
    angles.append([dir, exp_result, exp_delta])
    #print(angles)
print(angles)
#'''
#angles_met = [['10', 28.41, 0.37], ['15', 33.86, 0.39], ['20', 33.16, 0.35], ['25', 34.87, 0.36], ['5', 30.81, 0.33]]
#angles_temp = [['283', 34.77, 0.35], ['303', 37.05, 0.34], ['323', 37.14, 0.33], ['343', 36.51, 0.32], ['363', 36.1, 0.32], ['383', 31.81, 0.36], ['403', 31.21, 0.5]]
#angles_press = [['200', 32.21, 0.45], ['280', 31.09, 0.53], ['360', 33.74, 0.46], ['440', 29.96, 0.54], ['520', 33.87, 0.36], ['600', 32.78, 0.36]]
#angles_salt = [['0', 27.27, 0.5], ['126', 25.9, 0.4], ['158', 24.9, 0.52], ['189', 24.57, 0.54], ['221', 25.67, 0.39], ['284', 23.38, 0.42], ['32', 34.78, 0.46], ['63', 31.77, 0.37], ['95', 28.08, 0.42]]
#pres_old_salt = np.array([float(i[0]) for i in angles])*(1000/(8753*18.015))
#ift_old = np.array([i[1] for i in angles])
exit()
print(pres_old)
a, b = np.polyfit(pres_old, ift_old, deg=1)
x = np.arange(0, 2, 0.1)
y_est = a*x + b
y_est2 = a*pres_old + b
y_err = m.sqrt(((y_est2 - ift_old)**2).mean())
    

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)
ax.set_ylim(0, 90)
ax.set_title('Salt', loc="left", fontstyle='italic', fontsize=20)
ax.plot(x, y_est, '-', color ='tab:green')
ax.set_xlabel(r'$m_{NaCl}$'+ '[mol/kg]', fontstyle='italic')
ax.set_ylabel('Angle', fontstyle='italic')
# Add error area
ax.fill_between(x, y_est - y_err, y_est + y_err, alpha=0.2)
ax.errorbar(pres_old, ift_old, yerr=1.7, marker='*', color ='tab:brown', 
                linestyle='none', ecolor='tab:brown', elinewidth=0.8, capsize=4, capthick=1)
#x.plot(np.array([float(i[0]) for i in angles]), np.array([i[1] for i in angles]), 'o', color='tab:brown')
plt.show()