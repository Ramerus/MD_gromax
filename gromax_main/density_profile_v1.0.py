import gromax_lib as gl
from MDAnalysis.coordinates.XTC import XTCReader

'''
To use this script you should have xtc file with trajectories
and at least 1 gro file of your system
'''
exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\drop_lib\\quartz_methane'
#exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\density_profile'
dx = 0.12
dir = '10'
file_type = 'xtc'
num_files = 2000
fig_name = 'met_10'
solid = 4.0 #hight of mineral(quartz) for type = 'angle'
#solid = 0
mass, mol_mas = gl.forcefield_loader()

files_xtc = gl.find_files(exp_path, dir, file_type)

files_gro = gl.find_files(exp_path, dir, 'gro')
system_data = gl.system_desc(files_gro[0], mass, mol_mas)

try:
    exp_data = XTCReader(files_xtc[-1])
except:
    print('No .xtc file! Check experimental directory!')
    exit()
den_w, den_ar, den_alk, den_s, den_met = gl.density_profiles(exp_data, system_data, dx, num_files, type = 'angle', solid = solid)
gl.dens_show_1d(den_w, den_ar, den_alk,  den_s, den_met, exp_data[0].dimensions[1]*0.1+dx, fig_name, dx = dx)
#exit()
'''calculating dissolewed components
w_in_oil = 0
alk_in_w = 0
met_in_w = 0
met_in_oil = 0

for i in range(len(den_w)):
    if den_w[i] < 0.0017:
        w_in_oil += den_w[i]
    elif den_w[i] > 0.98:
        alk_in_w += den_alk[i]
        met_in_w += den_met[i]
    if den_w[i] <= 0.98:
        met_in_oil += den_met[i]
print("Dissolved water in oil: ", w_in_oil/sum(den_alk)*100, " %")
print("Dissolved oil in water: ", alk_in_w/sum(den_w)*100, " %")
print("Dissolved methane in water: ", met_in_w/sum(den_w)*100, " %")
print("PPM of methane in oil: ", met_in_oil/(sum(den_alk)+sum(den_ar)))
'''
density = 0
mid = 0
for i in range(int(9//dx), int(12//dx)):
    density += den_alk[i]+den_ar[i]+den_met[i]
    mid += 1
density_m = density/mid
print(density_m)