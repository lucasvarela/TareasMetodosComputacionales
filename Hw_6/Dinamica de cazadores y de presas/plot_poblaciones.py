
import numpy as np
import matplotlib.pyplot as plt
import sys

a=sys.argv[1]
b=a.split(".")
data = np.loadtxt(""+a+"")
plt.xlabel("Presa")
plt.ylabel("Depredador")
plt.title(b[0])
plt.plot(data[:,0], data[:,1])
plt.savefig(""+b[0]+".pdf")





