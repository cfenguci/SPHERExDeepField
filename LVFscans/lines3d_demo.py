import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
#theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
#z = np.linspace(-2, 2, 100)
#r = z**2 + 1
#x = r * np.sin(theta)
#y = r * np.cos(theta)

data=np.loadtxt('/home/cmb/computing/spherex/VLFscans/output/traj.txt');
x=data[:,0];
y=data[:,1];
z=data[:,2];
ax.plot(x, y, z)
ax.legend()

plt.show()
