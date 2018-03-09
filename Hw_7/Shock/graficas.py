import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

datos = np.loadtxt(filename)

x = datos[:,0]
u = datos[:,1]
p = datos[:,2]
rho_f = datos[:,3]

figura = plt.figure(figsize=(12.0, 17.0))

ax = figura.add_subplot(3, 1, 1)
ax.plot(x, u)
ax.set_title("Velocidad")
ax.set_xlabel("$x[m]$")
ax.set_ylabel("$u[m/s]$")

ax = figura.add_subplot(3, 1, 2)
ax.plot(x, p)
ax.set_title("Presion")
ax.set_xlabel("$x[m]$")
ax.set_ylabel("$p[Pa]$")

ax = figura.add_subplot(3, 1, 3)
ax.plot(x, rho_f)
ax.set_title("Densidad")
ax.set_xlabel("$x[m]$")
ax.set_ylabel(r"$\rho[kg/m^3]$")

plt.savefig(filename[:-4]+".pdf", format ='pdf', transparent=True)

