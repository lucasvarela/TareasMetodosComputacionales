from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

x = np.outer(np.ones(121), np.linspace(0,100,101))
t = np.outer(np.linspace(0,120,121), np.ones(101))
data = np.loadtxt(filename)
y = data;

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.plot_wireframe(t,x,y,rstride= 1, cstride =1,linewidth=0.2,  color ='b')
ax.set_xlabel('$t[s]$')
ax.set_ylabel('$x[m]$')
ax.set_zlabel('$y[m]$')

plt.savefig(filename[:-4] +'.pdf', format = 'pdf', transparent = True)
plt.close
