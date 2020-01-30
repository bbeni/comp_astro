#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob 

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D



def load_all_snapshots(every=1, end=-1):
	files = glob.glob('snapshot*.csv')
	files = sorted(files)[:end:every]
	#files.insert(0, 'initial.csv')

	positions = []
	dfs = []
	for filename in files:
		df = pd.read_csv(filename)
		dfs.append(df.copy())
		pos = df.loc[:,['x', 'y', 'z']].values
		positions.append(pos.copy())

	positions = np.array(positions)
	return dfs, positions

dfs, positions = load_all_snapshots(every=10, end=10000)
print(positions.shape)


# energy
def plot_energies():
	e_kins = []
	e_pots = []

	for df in dfs:
		e_kin = 0.5*np.dot(df.loc[:,'m'], np.square(df.loc[:,['vx', 'vy', 'vz']]).sum(axis=1))
		e_kins.append(e_kin)

		pos = df.loc[:,['x', 'y', 'z']].values
		m = df.loc[:, ['m']].values
		epsilon = df.loc[:, 'softening'].values
		
		e_pot = 0
		for i in range(pos.shape[0]):
			for j in range(pos.shape[0]):
				if i == j: 
					continue
				xij = pos[i] - pos[j]
				underij = np.sum(np.square(xij)) + np.square(epsilon[i])
				e_pot -= m[i]*m[j]/underij	# G = 1

		e_pots.append(e_pot)

	e_kins = np.array(e_kins)
	e_pots = np.array(e_pots).reshape((len(e_pots),))/2
	e = e_kins + e_pots

	print(e_pots)
	print(e_kins)
		
	print(e_kins.shape, e_pots.shape)

	plt.plot(e/e[0])
	plt.ylabel("dE/E")
	plt.xlabel("t")

	plt.show()

	plt.plot(e_kins, label='E_kin')
	plt.plot(e_pots, label='E_pot')
	plt.plot(e_pots + e_kins, label='E_tot')
	plt.legend()
	plt.xlabel('#time step')
	plt.ylabel('Energy')
	plt.show()

plot_energies()


# plot

for particle_nr in range(20):
	x = positions[:, particle_nr, 0]
	y = positions[:, particle_nr, 1]

	plt.plot(x, y, '--')
	#plt.scatter(x, y)

plt.show()


# animation

from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot(positions[0, :, 0], positions[0, :, 1], 'r.')

def init():
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    return ln,

def update(frame):
    xdata = positions[frame, :, 0]
    ydata = positions[frame, :, 1]
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=range(positions.shape[0]),
                    init_func=init, blit=True)
plt.show()


# 3d plot not working

def plot_particle_3d(particle_nr):
	xs = positions[:, particle_nr, 0]
	ys = positions[:, particle_nr, 1]
	zs = positions[:, particle_nr, 2]



	ax.scatter(xs, ys, zs, marker='o')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(10):
	plot_particle_3d(i)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()