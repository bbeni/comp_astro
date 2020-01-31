#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob 

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


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

def plot_energies():
	tracking_df = pd.read_csv("simulation_track.csv")
	print(tracking_df)

	etot = tracking_df["e_pot"] + tracking_df["e_kin"]

	plt.plot( tracking_df["e_kin"], label="ekin")
	plt.plot( tracking_df["e_pot"], label="epot")
	plt.plot( etot, label="e")
	plt.xlabel('#time step')
	plt.ylabel('Energy')
	plt.legend()
	plt.show()

	plt.plot(etot/etot[0])
	plt.ylabel("dE/E")
	plt.xlabel("#time step")
	plt.show()

def plot_trajectories(positions):
	for particle_nr in range(20):
		x = positions[:, particle_nr, 0]
		y = positions[:, particle_nr, 1]

		plt.plot(x, y, '--')
		plt.scatter(x[-1], y[-1])

	plt.show()

def plot_animation(positions, gif=None):

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

	if gif:
		ani.save(gif, writer='imagemagick')

	plt.show()


def plot_trajectories_3d(positions):
	
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


if __name__ == "__main__":

	plot_energies()

	dfs, positions = load_all_snapshots(every=2, end=200)	
	print(positions.shape)

	plot_trajectories(positions)
	plot_animation(positions, "asdf.gif")
	plot_trajectories_3d(positions)







