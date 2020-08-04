#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def draw_energy_contour(file_name):
    data = np.loadtxt(file_name)
    
    data = np.reshape(data, (-1, 11, 11))
    c, alpha, mean, std = [data[i, ...] for i in range(4)]
    norm = cm.colors.Normalize(vmax=-0.4, vmin=-0.5)
    plt.contourf(mean)
    plt.xlabel("c parameter")
    plt.ylabel("alpha parameter")
    plt.savefig(file_name+".png")

if __name__ == "__main__":
    draw_energy_contour("energy_surface.txt")    