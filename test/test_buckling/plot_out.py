

import h5py
import matplotlib.pyplot as plt
import os
import numpy as np

h5 = h5py.File("out.h5")

plot_eta0 = 1
plot_eta1 = 1
plot_eta2 = 1

for frame in range(1,41):

    frame_str = str(frame).zfill(6)
    fig = plt.figure(figsize=(6,6))
    phi = h5["/phi"]
    plt.contour(np.transpose(phi), levels=[0.2,0.4,0.6,0.8], colors="black")

    if (plot_eta0):
        data0 = np.transpose(h5["/eta0/"+frame_str])
        plt.contour(data0, levels=[0.5,0.9], colors="black")
        plt.contour(np.negative(data0), levels=[0.5,0.9], colors="black")
        plt.contourf(data0, levels=[0.5,0.75,1,2],cmap=plt.cm.Greys)
        plt.contourf(np.negative(data0), levels=[0.5,0.75,1,2],cmap=plt.cm.Greys)

    if (plot_eta1):
        data1 = np.transpose(h5["/eta1/"+frame_str])
        plt.contour(data1, levels=[0.5,0.9], colors="black")
        plt.contour(np.negative(data1), levels=[0.5,0.9], colors="black")
        plt.contourf(data1, levels=[0.5,0.75,1,2],cmap=plt.cm.Blues)
        plt.contourf(np.negative(data1), levels=[0.5,0.75,1,2],cmap=plt.cm.Blues)

    if (plot_eta2):
        data2 = np.transpose(h5["/eta2/"+frame_str])
        plt.contour(data2, levels=[0.5,0.9], colors="black")
        plt.contour(np.negative(data2), levels=[0.5,0.9], colors="black")
        plt.contourf(data2, levels=[0.5,0.75,1,2],cmap=plt.cm.Greens)
        plt.contourf(np.negative(data2), levels=[0.5,0.75,1,2],cmap=plt.cm.Greens)
    
    ax = fig.gca()
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])

    fig.savefig("displacement" + frame_str)
    os.system("convert -trim " + frame_str + ".png " + frame_str + ".png")

