
from mpl_toolkits.mplot3d import Axes3D
import h5py
import matplotlib.pyplot as plt
import os
import numpy as np
import pickle as pl


h5 = h5py.File("out.h5")

plot_eta0 = 1
plot_eta1 = 1
plot_eta2 = 1

frame = 40

frame_str = str(frame).zfill(6)
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
    
if (plot_eta0):
        eta0 = np.transpose(h5["/eta0/"+frame_str])

if (plot_eta1):
        eta1 = np.transpose(h5["/eta1/"+frame_str])

if (plot_eta2):
        eta2 = np.transpose(h5["/eta2/"+frame_str])

w_sub = np.transpose(h5["/w_sub/" + frame_str])
a = w_sub.shape
x = a[0]
y = a[1]
X = np.arange(x)
Y = np.arange(y)
Y, X = np.meshgrid(Y,X)

phi = np.transpose(h5["/phi"])
colors = np.empty(phi.shape, dtype=str)
for x1 in range(x):
   for y1 in range(y):
         colors[x1,y1] = "w"
         if (eta0[x1,y1] >= 0.8) or (eta0[x1,y1] <= -0.8):
          colors[x1,y1] = "g"
         if (eta1[x1,y1] >= 0.8) or (eta1[x1,y1] <= -0.8):
          colors[x1,y1] = "b"
         if (eta2[x1,y1] >= 0.8) or (eta2[x1,y1] <= -0.8):
          colors[x1,y1] = "c"
         if (phi[x1,y1] >= 0.3) and (phi[x1,y1] <= 0.7):
          colors[x1,y1] = "k"

ax.set_zlim(0, 100)
ax.plot_surface(X, Y, w_sub, facecolors= colors, rcount = 100, ccount = 100)

plt.show()

