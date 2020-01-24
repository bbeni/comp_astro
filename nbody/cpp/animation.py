#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob 

import matplotlib.pyplot as plt



def plot_snapshot(filename):
	df = pd.read_csv(filename)

	pos = df.loc[:,['x', 'y', 'z']]

	x = df.loc[ : , 'z']
	y = df.loc[ : , 'y']


	plt.plot(x, y, '.')


plot_snapshot('initial.csv')
plot_snapshot('snapshot0.csv')
plot_snapshot('snapshot1.csv')
plot_snapshot('snapshot2.csv')

plt.show()

def load_all_snapshots():
	files = glob.glob('snapshot*.csv')
	files.insert(0, 'initial.csv')

	positions = []
	for filename in files:
		df = pd.read_csv(filename)
		pos = df.loc[:,['x', 'y', 'z']].values
		positions.append(pos)

	positions = np.array(positions)

	return positions

positions = load_all_snapshots()
print(positions.shape)

# plot

n = 1
xs = positions[:, 0, n]
ys = positions[:, 1, n]
zs = positions[:, 2, n]

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(xs, ys, zs, marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()