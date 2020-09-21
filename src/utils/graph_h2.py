#! /usr/bin/env python

import numpy as np
# import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

data = np.loadtxt("h2_pes.txt")
print(data)

xdata = data[:, 0]
ydata = data[:, 1]

def func(r, de, a):
    return -de + de*(1-np.exp(-a*(r-0.74)))**2

popt, pcov = curve_fit(func, xdata, ydata, maxfev=800)
print(popt) 