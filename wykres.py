import numpy as np
import matplotlib.pyplot as pt

M = np.loadtxt("dane.txt")
pt.plot(M[:,0],M[:,1], fmt="r-")
pt.show()