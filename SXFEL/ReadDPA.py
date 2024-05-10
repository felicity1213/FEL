import numpy as np
import matplotlib.pylab as plt

# Change parameter according to your .in file
nslice = 1  # set to 1 if run in steady-state mode
harmonic = 1
npart = 8192

# Creat empty matrix
z1 = np.zeros(npart*nslice)
p1 = np.zeros(npart*nslice)
x1 = np.zeros(npart*nslice)
xp1 = np.zeros(npart*nslice)
y1 = np.zeros(npart*nslice)
yp1 = np.zeros(npart*nslice)
zpp = np.zeros(npart*nslice)

# Read & Write
with open('particle_1.out.dpa', 'rb') as f:
    BB = np.fromfile(f, dtype=np.double)
    BB = BB.reshape([6,npart*nslice])
    
for i in range(1, nslice+1):
    for j in range(npart):
        zpp[j+(i-1)*npart] = BB[1,j]
        z1[j+(i-1)*npart] = BB[1,j] / harmonic + 2*np.pi*(i-1)
        p1[j+(i-1)*npart] = BB[0,j]
        x1[j+(i-1)*npart] = BB[2,j]
        xp1[j+(i-1)*npart] = BB[4,j] / BB[0,j]
        y1[j+(i-1)*npart] = BB[3,j]
        yp1[j+(i-1)*npart] = BB[5,j] / BB[0,j]

z1 = z1/2/np.pi
# print(z1)
# print(p1)

# Visualization
plt.figure()
plt.plot(z1, p1, 'r.', markersize=0.5)
plt.show()