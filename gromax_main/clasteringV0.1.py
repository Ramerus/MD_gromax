import gromax_lib as gl
from MDAnalysis.coordinates.XTC import XTCReader

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\drop_lib\\quartz_methane'
#exp_path = 'C:\\Users\\peter\\Desktop\\files\\my_soft\\density_profile'
dx = 0.12
dir = '15'
file_type = 'xtc'
num_files = 2000
fig_name = 'met_5'
solid = 4.0 #hight of mineral(quartz) for type = 'angle'
#solid = 0

mass, mol_mas = gl.forcefield_loader()

files_xtc = gl.find_files(exp_path, dir, file_type)

files_gro = gl.find_files(exp_path, dir, 'gro')
system_data = gl.system_desc(files_gro[0], mass, mol_mas)

composition = gl.composition(exp_path, dir)


try:
    exp_data = XTCReader(files_xtc[-1])
except:
    print('No .xtc file! Check experimental directory!')
    exit()
coords = np.dot(np.array(exp_data[exp_data.n_frames-1]._pos), 0.1)
#exit()
#aro_data = gl.aromatic_data(composition, coords, system_data)
aro_data = gl.aro_data_old(coords, system_data)
#scaler = StandardScaler()
#x_scaled = scaler.fit_transform(aro_data)

#x_normal = normalize(x_scaled)
#x_normal = pd.DataFrame(x_normal)
#pca = PCA(n_components=2)
#x_principal = pca.fit_transform(x_normal)
x_principal = pd.DataFrame(aro_data)
#x_principal = x_normal
#x_principal.columns = ['V1', 'V2']
#print(x_principal.head())

dbscan = DBSCAN(eps=0.4, min_samples=20).fit(x_principal)
labels = dbscan.labels_
aro_data['cluster'] = dbscan.labels_
#print(aro_data.tail())



unique_labels = set(labels)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
core_samples_mask = np.zeros_like(labels, dtype=bool)
core_samples_mask[dbscan.core_sample_indices_] = True
#print(x_principal.loc[1, 'V1'])
#exit()
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = labels == k
    x = []
    y = []
    z = []
    for j in range(len(aro_data)):
        if class_member_mask[j] and core_samples_mask[j]:
            x.append(aro_data.loc[j, 'x'])
            y.append(aro_data.loc[j, 'y'])
            z.append(aro_data.loc[j, 'z'])

    #xy = x_principal[class_member_mask & core_samples_mask]
    ax.scatter(
        x,
        y,
        z,
        "o",
        color=tuple(col),
        edgecolor="k",
        s=14,
    )
    x = []
    y = []
    z = []
    for j in range(len(aro_data)):
        if (class_member_mask[j]) and (~core_samples_mask[j]):
             x.append(aro_data.loc[j, 'x'])
             y.append(aro_data.loc[j, 'y'])
             z.append(aro_data.loc[j, 'z'])
    ax.scatter(
        x,
        y,
        z,
        "o",
       color=tuple(col),
        edgecolor="k",
        s=14,
    )
ax.set_xlim(0, exp_data[0].dimensions[0]*0.1)
ax.set_ylim(0, exp_data[0].dimensions[1]*0.1)
ax.set_zlim(0, exp_data[0].dimensions[2]*0.1)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.title(f"Estimated number of clusters: {n_clusters_}")
plt.show()
#'''