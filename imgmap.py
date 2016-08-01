#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import misc

rawimg = misc.imread('andromeda.jpg', flatten=True)
smimg = misc.imresize(rawimg, (181,361))

plt.imshow(smimg, cmap=plt.cm.gray)
plt.show()

data = np.genfromtxt("tfunc.csv", delimiter=",")[:,4:]

timg = np.copy(smimg)
for point in data:
    timg[int(point[0]), int(point[1])] = smimg[int(point[2]), int(point[3])]

plt.imshow(timg, cmap=plt.cm.gray)
plt.show()
