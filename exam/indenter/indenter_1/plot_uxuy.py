

import h5py
import matplotlib.pyplot as plt
import os
import numpy as np

h5 = h5py.File("out.h5")

plot_ux = 0
plot_uy = 1

for frame in range(1,41):

    frame_str = str(frame).zfill(6)
    fig = plt.figure(figsize=(6,6))
    phi = h5["/phi"]
    plt.contour(phi, levels=[0.2,0.4,0.6,0.8], colors="black")

    if (plot_ux):
        data0 = np.transpose(h5["/ux/"+frame_str])
        plt.contour(data0, levels=[0.5,0.9], colors="black")
        plt.contour(np.negative(data0), levels=[0.5,0.9], colors="black")
        plt.contourf(data0, levels=[0.001,0.005,0.01,0.02,0.05],cmap=plt.cm.Greens)
        plt.contourf(np.negative(data0), levels=[0.001,0.005,0.01,0.02,0.05],cmap=plt.cm.Blues)

    if (plot_uy):
        data1 = np.transpose(h5["/uy/"+frame_str])
        plt.contour(data1, levels=[0.5,0.9], colors="black")
        plt.contour(np.negative(data1), levels=[0.5,0.9], colors="black")
        plt.contourf(data1, levels=[0.001,0.005,0.01,0.02,0.05],cmap=plt.cm.Greens)
        plt.contourf(np.negative(data1), levels=[0.001,0.005,0.01,0.02,0.05],cmap=plt.cm.Blues)

    
    ax = fig.gca()
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])

    fig.savefig("displacement"+frame_str)
    os.system("convert -trim " + frame_str + ".png " + frame_str + ".png")

