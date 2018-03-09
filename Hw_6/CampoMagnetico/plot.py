
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D

a=sys.argv[1]
b=a.split(".")
datos = np.loadtxt(""+a+"")
x = datos[:,0]
y = datos[:,1]
z = datos[:,2]
# se usa la libreria de matplot.lib para hacer la esfera y el circulo.
deta= np.linspace(0,2*np.pi,100)
fig = plt.figure(figsize=(12.0, 20.0))
ax = fig.add_subplot(2, 1, 1)
l=plt.plot( np.cos(deta), np.sin(deta) )

l = ax.plot(x, y)
plt.xlim(-2.5,2.5)
plt.ylim(-2.5,2.5)
ax.grid(True)
ax.set_xlabel('x')
ax.set_ylabel('y')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
ax = fig.add_subplot(2, 1, 2, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z=  np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(X, Y, Z,  rstride=4, cstride=4, color='r')
surf = ax.plot(x, y, z,label='trayectoria en x-y-z')

plt.savefig(""+b[0]+".pdf")

