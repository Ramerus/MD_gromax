import numpy as np
import math as m
import scipy as sc
import time
import numpy as np
from os import listdir, getcwd, chdir, mkdir
from os.path import isfile, join, isdir
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import difflib
from tqdm import tqdm
from numpy import save
from MDAnalysis.coordinates.XTC import XTCReader
from pylab import mean, sqrt, empty, newaxis

def density_profiles(exp_data, system_data, dx, num_files, type, solid = 0):
    '''
    type - 'angle', 'ift'
    '''
    number_frames =  exp_data.n_frames
    step_number = number_frames // num_files
    
    if step_number == 0:
        step_number = 1
        #num_files = number_frames
    num_files = number_frames // step_number
    xs, ys, zs  = exp_data[0].dimensions[0]*0.1, exp_data[0].dimensions[1]*0.1, exp_data[0].dimensions[2]*0.1
    if type == 'ift':
        zs += dx
        den_water = np.zeros(m.floor(zs/dx))
        den_aro = np.zeros((m.floor(zs/dx))) 
        den_alkane = np.zeros((m.floor(zs/dx))) 
        den_salt = np.zeros((m.floor(zs/dx)))
        den_methane = np.zeros((m.floor(zs/dx)))
    elif type == 'angle':
        xs += dx
        den_water = np.zeros(m.floor(xs/dx))
        den_aro = np.zeros((m.floor(xs/dx))) 
        den_alkane = np.zeros((m.floor(xs/dx))) 
        den_salt = np.zeros((m.floor(xs/dx)))
        den_methane = np.zeros((m.floor(xs/dx)))   

    if num_files > number_frames:
        step_number = 1
    print('COLLECTING DENSITIES')
    for i in tqdm(range(0, number_frames-1, step_number)):
        den_w, den_ar, den_alk, den_s, den_met= dens_map_1d(xs, ys, zs, exp_data[i]._pos, system_data, dx = dx, dz = dx, solid = solid, type = type)
        den_water += den_w/num_files
        den_aro += den_ar/num_files
        den_alkane += den_alk/num_files
        den_salt += den_s/num_files
        den_methane += den_met/num_files

    return den_water, den_aro, den_alkane, den_salt, den_methane

def forcefield_loader():
    '''
    Making dataframes with correct atom names, molecular names, masses
    mass -- equal to trappe forcefield file
    mol_mas -- masses of atoms
    '''
    exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\density_profile'
    mass = pd.read_excel(exp_path+'\\trappe.xlsx')
    mol_mas = pd.read_excel(exp_path+'\\molecular.xlsx')

    for i in range(len(mol_mas)):
        mol_mas.loc[i, 'atom'] = str(mol_mas.loc[i, 'atom']).replace(' ', '')
        mol_mas.loc[i, 'type'] = str(mol_mas.loc[i, 'type']).replace(' ', '')
    for i in range(len(mass)):
        mass.loc[i, 'atom'] = str(mass.loc[i, 'atom']).replace(' ', '')
    print('FORCEFIELD LOADED')
    return mass, mol_mas

def find_files(exp_path, dir, type):
    '''
    Find all files with some types(gro or xtc)
    in:
    exp_path -- path to experiment
    dir -- directory of experiment
    type -- type of file
    out:
    files -- all files with needed type in path
    '''
    
    chdir(join(exp_path,dir))
    path = join(exp_path,dir)
    files = listdir(join(path))
    files = [f for f in listdir(path) if (isfile(join(path, f))) and (('.'+type) in f)]
    
    return files

def composition(exp_path, dir):
    chdir(join(exp_path,dir))
    path = join(exp_path,dir)
    files = listdir(join(path))
    files = [f for f in listdir(path) if (isfile(join(path, f))) and (('system') in f)]
    name = []
    number = []
    with open(files[0]) as F:
        x = F.readline()
        while '[ molecules ]' not in x:
            x = F.readline()
        F.readline()
        x = F.readline()
        while '' != x:
            tup = x.split(' ')
            name.append(tup[0])
            number.append(int(tup[1].split('\n')[0]))
            x = F.readline()
        dic = {'name': name, 'number': number}
        compo = pd.DataFrame(dic)
    return compo

def atom_mass(atom, df, mol_mas, name):
    atom = atom.replace(' ', '')
    if atom == 'NA':
        atom = 'NA+'
    try:
        tr_atom = mol_mas[(mol_mas['name'] == name) & (mol_mas['atom'] == atom)]['type'].values[0]
    except:
        print('Warning!! Error in forcefield! ', atom, name)
        exit()
    mass = float(df[tr_atom == df['atom']]['mass'])                  
    return mass

def system_desc(file, mass_data, mol_mas):
    '''
    Make list of typles [(mol, atom_mass)]
    file -- filename of any .gro file
    '''
    system_data = []
    with open(file) as F:
        F.readline()
        n_lines = int(F.readline())
        for i in range(n_lines):
            x = F.readline()
            system_data.append([x[5: 10], float(atom_mass(x[10:15], 
                                                          mass_data, mol_mas, x[5: 10]))])
    print('System info loaded')
    return system_data

def mol_desc(file):
    '''
    Make list of  [mol_names]
    file -- filename of any .gro file
    '''
    system_data = []
    with open(file) as F:
        F.readline()
        n_lines = int(F.readline())
        for i in range(n_lines):
            x = F.readline()
            system_data.append(x[5: 10])
    print('System info loaded')
    return system_data


def dens_map_1d(xs, ys, zs, dd, system_data, solid=0, dx=0.4, dz=0.4, type = 'ift'):
    """Form 3 multy-arrays -- density profile of water, oil and salt

    Keyword arguments:
    xs, ys, zs -- size of sim.box 
    solid -- hight of mineral
    dd -- data with coords of atoms
    system_data -- data with moll names and atom masses
    dx -- x, y size of box for density calculation (default 0.4)
    dz -- z size of box for density calculation (default 0.4)
    type -- 'ift' or 'angle' because we need go on z or x
    """
    z0 = solid #zero for atoms
    zmax = zs - z0 # max z coord for atoms
    if type == 'ift':
        v_slice = dx*xs*ys #volume of one slice
        den_water = np.zeros(m.floor(zmax/dx))
        den_aro = np.zeros((m.floor(zmax/dx))) 
        den_alkane = np.zeros((m.floor(zmax/dx))) 
        den_salt = np.zeros((m.floor(zmax/dx)))
        den_methane = np.zeros((m.floor(zmax/dx)))
    
        for i in range(len(dd)):
            z = int((dd[i, 2]*0.1)//dx)
            if system_data[i][0] == 'WATER':
                den_water[z] += system_data[i][1]
            if system_data[i][0] in ['ASPHA', 'PYREN', 'NAPHT', 'PRBEN', 'BENZE']:
                den_aro[z] += system_data[i][1]
            if system_data[i][0] in ['SODIU', 'CHLOR']:
                den_salt[z] += system_data[i][1]
            if system_data[i][0] in ['TRCNT', 'NONCS', 'TETCS', 'NONDE', 'TETRD',
                                      'NONAN', 'OCTAN', 'HEPTN', 'BUTAN']:
                den_alkane[z] += system_data[i][1]
            if system_data[i][0] == 'METAN':
                den_methane[z] += system_data[i][1]
    elif type == 'angle':
        v_slice = dx*zmax*ys #volume of one slice
        den_water = np.zeros(m.floor(xs/dx))
    
        den_aro = np.zeros((m.floor(xs/dx))) 
        den_alkane = np.zeros((m.floor(xs/dx))) 
    
        den_salt = np.zeros((m.floor(xs/dx)))
        den_methane = np.zeros((m.floor(xs/dx)))
    
        for i in range(len(dd)):
            z = int((dd[i, 0]*0.1)//dx)
            if system_data[i][0] == 'WATER':
                den_water[z] += system_data[i][1]
            if system_data[i][0] in ['ASPHA', 'PYREN', 'NAPHT', 'PRBEN', 'BENZE']:
                den_aro[z] += system_data[i][1]
            if system_data[i][0] in ['SODIU', 'CHLOR']:
                den_salt[z] += system_data[i][1]
            if system_data[i][0] in ['TRCNT', 'NONCS', 'TETCS', 'NONDE', 'TETRD',
                                    'NONAN', 'OCTAN', 'HEPTN', 'BUTAN']:
                den_alkane[z] += system_data[i][1]
            if system_data[i][0] == 'METAN':
                den_methane[z] += system_data[i][1]

    gr_to_rho = 6.0221366516752e+2
    return (den_water)/(v_slice*gr_to_rho), (den_aro)/(v_slice*gr_to_rho), (den_alkane)/(v_slice*gr_to_rho), (den_salt)/(v_slice*gr_to_rho), (den_methane)/(v_slice*gr_to_rho)


def mas_remake(middle, mass):
    l = int(middle*0.5 - 0.5)
    r = int(middle*0.5 + 0.5)
    if (middle % 2 == 0):
        l = r = int(middle*0.5)
    n = len(mass)
    ans = []
    k = 1
    while (l >= 0) and (r<=n-1):
        ans.append(0.5*(mass[l]+mass[r]))
        l -= 1
        r += 1
    if (l < 0) and (k<2):
        k +=1
        middle = n-1-r
        l = int(middle*0.5 - 0.5)+r
        r = int(middle*0.5 + 0.5)+r
        if (middle % 2 == 0):
            l = r = int(middle*0.5)
        while (l >= 0) and (r<=n-1):
            ans.append(0.5*(mass[l]+mass[r]))
            l -= 1
            r += 1
    if (r > n-1) and (k<2):
        k +=1
        middle = l
        l = int(middle*0.5 - 0.5)
        r = int(middle*0.5 + 0.5)
        if (middle % 2 == 0):
            l = r = int(middle*0.5)
        while (l >= 0) and (r<=n-1):
            ans.append(0.5*(mass[l]+mass[r]))
            l -= 1
            r += 1
    return ans

def dens_show_1d(den_water, den_aro, den_alkane, den_salt, den_methane, zs, file, dx=0.4):
    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim(0, 1.4)
    ax.set_title(r"$CH_4 = 5\%$", loc="left", fontstyle='italic', fontsize=16)
    ax.set_xlabel('x[nm]')
    ax.set_ylabel(r'${\rho}[g/cm^3]$')
    x = np.array(list(range(len(den_water))))*dx
    ax.plot(x, den_water, '-', color ='tab:blue', label='Water')
    ax.plot(x, den_aro, '-', color ='y', label='Aromatic')
    ax.plot(x, den_alkane, '-', color ='tab:grey', label='Alkane')
    ax.plot(x, den_methane, '-', color ='red', label='Methane')
    ax.legend()
# change the color of the top and right spines to opaque gray
    ax.spines['right'].set_color((.8,.8,.8))
    ax.spines['top'].set_color((.8,.8,.8))
    # tweak the axis labels
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    xlab = ax.xaxis.get_label()
    ylab = ax.yaxis.get_label()
    xlab.set_style('italic')
    xlab.set_size(16)
    ylab.set_style('italic')
    ylab.set_size(16)
    # tweak the title
    ttl = ax.title
    ttl.set_weight('bold')
    ttl.set_size(16)
    fig.tight_layout()
    file_save = file + ".jpg"
    plt.savefig(file_save, dpi= 600)
    plt.close()    
    return

def aro_data_old(coords, system_data):
    x = []
    y = []
    z = []
    for i in range(len(coords)):
        if system_data[i][0] in ['ASPHA', 'PYREN', 'NAPHT', 'PRBEN', 'BENZE']:
            x.append(coords[i, 0])
            y.append(coords[i, 1])
            z.append(coords[i, 2])
    data = {'x': x, 'y': y, 'z': z}
    df = pd.DataFrame(data)
    return df

def aromatic_data(composition, coords, system_data):
    x = []
    y = []
    z = []
    exp_path = 'path to exp'
    mol_length = pd.read_excel(exp_path+'\\mol_length.xlsx')
    
    aromatic = ['perylene', 'trimetphen', 'pyrene', 'phenanthrene', 'dinaphtalene', 'propbenz', 'toluene']
    i = 0
    for mol in composition['name']:
        lenght = int(mol_length[mol_length['name'] == mol]['atom_numbers'])*int(composition[composition['name'] == mol]['number'])
        if mol in aromatic:
            x_sum, y_sum, z_sum, mol_mas = 0, 0, 0, 0
            for _ in range(int(composition[composition['name'] == mol]['number'])):
                for j in range(i, i+int(mol_length[mol_length['name'] == mol]['atom_numbers'])):
                    x_sum, y_sum, z_sum = x_sum+coords[i, 0]*system_data[i][1], y_sum+coords[i, 1]*system_data[i][1], z_sum+coords[i, 2]*system_data[i][1]
                    mol_mas += system_data[i][1]
                x.append(x_sum/mol_mas)
                y.append(y_sum/mol_mas)
                z.append(z_sum/mol_mas)
                i += int(mol_length[mol_length['name'] == mol]['atom_numbers'])
            
        else:
            i += lenght   
    
    data = {'x': x, 'y': y, 'z': z}
    df = pd.DataFrame(data)
    return df

class circle_fit():
    def __init__(self, xr, yr):
        self.x = xr
        self.y = yr
    
    def calc_R(self, xc, yc):
        """ calculate the distance of each 2D points from the center c=(xc, yc) """
        return sqrt((self.x-xc)**2 + (self.y-yc)**2)

    def f_2b(self, c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = self.calc_R(*c)
        return Ri - Ri.mean()

    def Df_2b(self, c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc     = c
        df2b_dc    = empty((len(c), self.x.size))

        Ri = self.calc_R(xc, yc)
        df2b_dc[ 0] = (xc - self.x)/Ri # dR/dxc
        df2b_dc[ 1] = (yc - self.y)/Ri                   # dR/dyc
        df2b_dc       = df2b_dc - df2b_dc.mean(axis=1)[:, newaxis]

        return df2b_dc
    def fitting(self):
        # coordinates of the barycenter
        x_m = mean(self.x)
        y_m = mean(self.y)
        center_estimate = x_m, y_m
        
        center_2b, ier = sc.optimize.leastsq(self.f_2b, center_estimate, Dfun=self.Df_2b, col_deriv=True)
        self.xc_2b, self.yc_2b = center_2b
        Ri_2b        = self.calc_R(self.xc_2b, self.yc_2b)
        self.R_2b         = Ri_2b.mean()
        delta_R = sum(np.fabs(Ri_2b - self.R_2b))/len(Ri_2b)
        
        residu_2b    = sum((Ri_2b - self.R_2b)**2)
        residu2_2b   = sum((Ri_2b**2-self.R_2b**2)**2)
        return self.xc_2b, self.yc_2b, self.R_2b, delta_R


class SAC():
    def __init__(self, exp_data, num_files, system_data, shift, solid):
        self.exp_data = exp_data
        self.num_files = num_files
        self.system_data = system_data
        self.shift = shift
        self.solid = solid
    def angle_calc(self, x, xc, yc):
        """Calculates the gradient, the height of tangent line and
        angle of contact angle

        Keyword arguments:
        solid -- the hight of the mineral surface
        x -- touch point coordinate on x axis
        xc -- circle center coordinate on x axis
        yc -- circle center coordinate on y axis
        """
        # определение координат нормального вектора(между касанием и центром окружности)
        xn = xc - x
        yn = yc - self.solid
        # свободный член уравнения касательной
        ck = x*xn/yn + self.solid
        # коэффициент наклона прямой
        ak = -xn/yn
        #считаем угол
        angle = m.degrees(m.atan(-xn/yn))
        self.xn = xn
        self.yn = yn
        return ck, ak, angle
    def surface_saddle(self, dp, op, f=0.4, mn=0.15, dx=0.4, dz=0.4):
        """Form 4 arrays -- 2 array of coords for water-oil surface points
        to the left and to the right 

        Keyword arguments:
        dp -- density profile of water
        op -- oil existanse
        f -- max % of water in box on surface (default 0.4)
        mn -- min % of water in box on surface (default 0.1)
        dx -- x, y size of box for density calculation (default 0.4)
        dz -- z size of box for density calculation (default 0.4)
        """
        md = dp.max() #maximum density of water
        x, y, z = dp.shape
        sxr = []
        syr = []
        sxl = []
        syl = []
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    # (low density of water) & (have min ammout of water) & (have oil)
                    if (dp[i, j, k] > mn*md) and (dp[i, j, k] < f*md) and \
                        (op[i, j, k] == 1):
                        if dx*i+dx/2>self.xs/2:
                            sxr.append(dx*i+dx/2)
                            syr.append(dz*k+dz/2+self.solid)
                        else:
                            sxl.append(dx*i+dx/2)
                            syl.append(dz*k+dz/2+self.solid)
        sxr = np.array(sxr)
        syr = np.array(syr)
        sxl = np.array(sxl)
        syl = np.array(syl)
        return sxr, syr, sxl, syl
    def dens_map_ci(self, data, dx=0.4, dz=0.4):
        """Form 2 arrays -- density profile of water and oil existanse

        Keyword arguments:
        xs, ys, zs -- size of sim.box 
        solid -- hight of mineral
        dd -- data with coords of atoms
        dx -- x, y size of box for density calculation (default 0.4)
        dz -- z size of box for density calculation (default 0.4)
        """
        z0 = self.solid #zero for atoms
        zmax = self.zs - z0 # max z coord for atoms
        xs = self.xs+dx 
        ys = self.ys+dx
        zmax = zmax+dz
        #density profile of simulation
        dmap = np.zeros((m.floor(xs/dx), m.floor(ys/dx), m.floor((zmax)/dz)))
        #oil boolean profile
        dprofile = np.zeros((m.floor(xs/dx), m.floor(ys/dx), m.floor((zmax)/dz))) 
        for i in range(len(data)):
            x = int(data[i, 0]*0.1//dx)
            y = int(data[i, 1]*0.1//dx)
            z = int((data[i, 2]*0.1-z0)//dz)
            if self.system_data[i] == 'WATER':
                dmap[x, y, z] += 1
            elif (self.system_data[i] != 'QUARZ') and (z < zmax // dz):
                dprofile[x, y, z] = 1
        return dmap, dprofile
    def contact_angle_saddle(self, xc, yc, rc):
        """Calculates the gradient, the height of tangent line and
        angle of contact angle for a saddle

        Keyword arguments:
        solid -- the hight of the mineral surface
        xc -- circle center coordinate on x axis
        yc -- circle center coordinate on y axis
        rc -- radius of the circle
        """
        # координата по оси х, пересечения
        y0 = self.solid-yc
        x1 = xc - m.sqrt(abs(rc**2-(y0)**2))
        x2 = xc + m.sqrt(abs(rc**2-(y0)**2))
        if abs(self.xs*0.5-x2) < abs(self.xs*0.5-x1):
            x1 = x2
        ck, ak, angle = self.angle_calc(x1, xc, yc)
        return ck, ak, angle 
    def ca_saddle(self, data,  f=0.4, mn=0.15, dx=0.4, dz=0.4):
        """Calculates CA, coords of surface and circle parametrs for a saddle 

        Keyword arguments:
        solid -- the hight of the mineral surface
        shift -- the middle of the oil part (10.0 for quarts)
        drop_path -- the path to the .gro file
        """
        dp, oi = self.dens_map_ci(data) #making density profiles
        xr, yr, xl, yl = self.surface_saddle(dp, oi) #find surface points
        #try:
        for i in range(1):
            angle_right = circle_fit(xr, yr) #fitting with cilinder 
            xcr , ycr, rcr, delta_R_r = angle_right.fitting()
            ckr, akr, ar = self.contact_angle_saddle(xcr, ycr, rcr) #counting CA
            delta_angle = (rcr*delta_R_r*self.yn)/(self.yn**2+self.xn**2)
            angle_left = circle_fit(xl, yl) #fitting with cilinder 
            xcl , ycl, rcl, delta_R_l = angle_left.fitting()
            ckl, akl, al = self.contact_angle_saddle(xcl, ycl, rcl) #counting CA
            delta_angle += (rcl*delta_R_l*self.yn)/(self.yn**2+self.xn**2)
            delta_angle /= 2
            angle = 0.5*((ar*ar)**0.5+al) #mean angle
            xd = np.hstack((xr, xl))
            yd = np.hstack((yr, yl))
            coords = np.vstack((xd, yd))
            return angle, delta_angle
        #except Exception:
        #    print('точек нет')
        #    return 0
    def analyse(self):
        number_frames =  self.exp_data.n_frames
        step_number = number_frames // self.num_files
    
        if step_number == 0:
            step_number = 1
            #num_files = number_frames
        self.num_files = number_frames // step_number
        
        #self.xs, self.ys, self.zs  = self.exp_data[0].dimensions[0]*0.1, self.exp_data[0].dimensions[1]*0.1, self.exp_data[0].dimensions[2]*0.1
        angles = []
        deltas = []
        for i in tqdm(range(0, number_frames-1, step_number)):
            self.xs, self.ys, self.zs  = self.exp_data[i].dimensions[0]*0.1, \
                self.exp_data[i].dimensions[1]*0.1, self.exp_data[i].dimensions[2]*0.1
            angle, delta  = self.ca_saddle(self.exp_data[i]._pos, f=0.6, mn=0.15, dx=0.3, dz=0.3)
            angles.append(angle) #contact angle
            deltas.append(delta)
        return round(np.mean(angles), 2), round(np.mean(deltas), 2)
