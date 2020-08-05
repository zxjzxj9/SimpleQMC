#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def draw_energy_contour(file_name):
    data = np.loadtxt(file_name)
    # data: 11x11x4
    data = np.reshape(data, (11, 11, -1))
    c, alpha, mean, std = [data[..., i] for i in range(4)]
    #print(mean[3])
    #import sys; sys.exit()
    #norm = cm.colors.Normalize(vmax=-0.4, vmin=-0.5)
    #print(np.max(mean), np.min(mean))
    plt.contourf(c, alpha, mean)
    plt.xlabel("c parameter")
    plt.ylabel("alpha parameter")
    plt.savefig(file_name+".png")

if __name__ == "__main__":
    x = np.linspace(-1.0, 1.0, 11, endpoint=True)
    y = np.linspace(0.2, 1.8, 11, endpoint=True)
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    print(xv[3])
    # print(yv)
    draw_energy_contour("energy_surface.txt")    