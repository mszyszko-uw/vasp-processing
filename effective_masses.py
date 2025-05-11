# Autor Maciej Szyszko
# Skrypt do post-processingu - obróbki i rysowania obliczonych mas efektywnych z kodu EMC
# skrypt jest niezbyt ogólny i dosyć chaotyczny, jak będziemy chcieli to dodawać tą funkcjonalność to go poprawie


import os 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib


# wczytywanie mas efektywnych z emc
def get_matrix_from_file(file):

    with open(file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 25:
            for i in range(len(lines)):
                if lines[i][:2] == '->':
                    matrix = []
                    matrix.append(lines[i+2].split())
                    matrix.append(lines[i+3].split())
                    matrix.append(lines[i+4].split())
                    matrix = np.array(matrix, dtype=float)
                    #matrix = np.linalg.inv(matrix[:2, :2])
                    return matrix
        else:
            return None

k1 = np.array([0.105563180, 0.060946935, 0.000000000])
k1 = k1/np.linalg.norm(k1)
k2 = np.array([0.000000000, 0.121893870, 0.000000000])
k2 = k2/np.linalg.norm(k2)
P = np.array([k1[:2], k2[:2]])
#print(np.linalg.inv(P))

a1 = np.array([9.4659999999999993,    0.0000000000000000])
a2 = np.array([-4.7329999999999997,    8.1977960000000003])
main_path = os.getcwd()

subdirectories = os.listdir(main_path)
subdirectories = [os.path.join(main_path, item) for item in subdirectories if os.path.isdir(os.path.join(main_path, item))]
#print(subdirectories)

VB = {}
CB = {}

for subdir in subdirectories:
    configurations = os.listdir(subdir)
    for con in configurations:
        files = os.listdir(os.path.join(subdir, con))
        files = [os.path.join(subdir, con, item) for item in files if item.startswith('emcpy')]
        for file in files:
            matrix = get_matrix_from_file(file)
            if type(matrix) == np.ndarray:
                if matrix[0, 0] >= 0:
                    CB[subdir+"/"+con] = matrix
                else:
                    VB[subdir+"/"+con] = matrix


CB = {key: value for key, value in CB.items()}
VB = {key: value for key, value in VB.items()}

print(CB)
holes = {}
electrons = {}
#print(CB)

for key in CB.keys():
    eigenval, eigenvec = np.linalg.eig(CB[key])
    electrons[key] = sorted(eigenval)
    if key == '/home/maciej/Documents/2D_Materials/Masy efektywne/3x3_EMC/x-0.667/C00001':
        print(eigenval)
        print(eigenvec)
        print()
        m = np.array(eigenvec[0]).T
        m = np.linalg.inv(P)*m
        n = np.array(eigenvec[1])
        n = np.linalg.inv(P)*n
        max = m[0]*a1 + m[1]*a2
        max = max/np.linalg.norm(max)
        min = n[0]*a1 + n[1]*a2
        min = min/np.linalg.norm(min)
        #print("MAX", max)
        #print("MIN", min)
        plt.plot([0, eigenvec[0][0]], [0, eigenvec[0][1]])
        plt.plot([0, eigenvec[1][0]], [0, eigenvec[1][1]])
        plt.show()
    eigenval, eigenvec = np.linalg.eig(VB[key])
    holes[key] = sorted(eigenval)

#print(electrons.values())

e_max = []
e_min = []
h_max = []
h_min = []
e_anisotropy = []
h_anisotropy = []
x = []
for key in electrons.keys():
    x.append(float(key.split('/')[-2][-5:]))
    e_max.append(electrons[key][1])
    e_min.append(electrons[key][0])
    e_anisotropy.append(abs((electrons[key][1] - electrons[key][0])/electrons[key][1]))
    h_min.append(holes[key][0])
    h_max.append(holes[key][1])
    h_anisotropy.append(abs((holes[key][0] - holes[key][1])/holes[key][0]))
    
e_min = np.array(e_min)
e_max = np.array(e_max)
h_max = np.array(h_max)
h_min = np.array(h_min)

font = {'size': 14}
matplotlib.rc('font', **font)
plt.scatter(x, (e_max + e_min)/2, label='electrons', marker='^')
plt.scatter(x, -(h_max + h_min)/2, label='hole')
#plt.scatter(x, e_min, label='e min', marker='v')
#plt.scatter(x, h_max, label='h max', marker='*')
#plt.scatter(x, (h_min + h_max)/2, label='valence band')
plt.legend()
plt.xlabel('x')
plt.ylabel('Carrier mass [m$_{e}$]')
plt.show()

#plt.scatter(x, (h_max + h_min)/2, label='conduction band', marker='^')
#plt.scatter(x, e_min, label='e min', marker='v')
#plt.scatter(x, h_max, label='h max', marker='*')
#plt.scatter(x, (h_min + h_max)/2, label='valence band')
#plt.legend()
#plt.xlabel('x')
#plt.ylabel('Electron mass [m$_{e}$]')
#plt.show()

plt.scatter(x, e_anisotropy, label='electrons', marker='^')
plt.scatter(x, h_anisotropy, label='holes')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('Anisotropy')
plt.show()
