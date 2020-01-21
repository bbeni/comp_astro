#!/usr/bin/env python3

'''
The data file is made of an header and sequential arrays(each array as a single column data entry), as follows: Header
Number of particles(N), Number of gas particles(N*=0 there is actually no gas particle), Number of star particles

Arrays(running in i=1,...,N)
Masses[i]
x[i]
y[i]
z[i]
Vx[i]
Vy[i]
Vz[i]
softening[i]
potential[i]

STEP 1
Preliminarly, verify the form of the density function rho(r) by constructing it from the particles and compare it with the analytical density function described in the original paper by Hernquist (from 1990 on Astrophysical Journal available on the web). Add Poissonian error bars to the numerical density profile.
Note that we use a system of units in which G=1. Start by assuming units of length and mass for your calculations (units of velocity and time follow automatically from the assumption G=1).
Poissonian error: The exercise requires to compare the Ernquist density profile with the density profile of the data. The simplest way to do that is to divide the space in spherical bins (shells) and count how many particles you have in all the bins, as you would do when building a histogram. Then you can compare these values to the expected values, i.e. the average amount of particles you would expect in each bin given the Ernquist density profile. In doing that you should consider that the number of particles you count in a given shell is a random variable that follows a Poisson distribution with lambda = expected number of particles (average value). When you compare the 2 values you would need some error estimate to evaluate how similar they are, and the standard error you have is the standard deviation of the Poisson distribution, i.e. sqrt(lambda). Please be careful about how you choose your bins, in order to have reasonable results.

STEP 2
Compute the direct N-Body forces between particles (note the array potential[i] is not needed for this purpose). Start by assuming a softening of the order of the mean interparticle separation in the system, then repeat the force calculation by experimenting with different values of the softening.
[analytical vs numerical graph]
To check the direct force calculation result, and its dependence on softening choice, at a given radius r compute the analytical force expected based on the application of Newton's second theorem for spherical potentials and compare it to the average force felt by particles in a shell at the same radius, plotting the result (use the book "Galactic Dynamics" by Binney and Tremaine as main reference for the theoretical notions, in particular sec. 2.2 (most recent version of the book) or 2.1 (1987 version).
Compute the relaxation timescale of the numerical model given the number of particles and the physical crossing timescale (use the half-mass radius R_hm and the circular velocity computed at the half-mass radius, vc = sqrt(GM(<R_hm)/R_hm). Do you expect the value of the gravitational softening to change the relaxation timescale? In particular, do you expect it to increase or decrease if the softening is increased above the interparticle separation? Can you expain why?

'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def parse_data(filename='data.ascii'):
	def next_line_as_float(file):
		return float(next(file).strip())

	with open(filename) as f:
		# header
		N, _, _ = next(f).split()
		N = int(N)

		masses = [next_line_as_float(f) for _ in range(N)]
		x = [next_line_as_float(f) for _ in range(N)]
		y = [next_line_as_float(f) for _ in range(N)]
		z = [next_line_as_float(f) for _ in range(N)]
		vx = [next_line_as_float(f) for _ in range(N)]
		vy = [next_line_as_float(f) for _ in range(N)]
		vz = [next_line_as_float(f) for _ in range(N)]
		softening = [next_line_as_float(f) for _ in range(N)]
		potential = [next_line_as_float(f) for _ in range(N)]

	# np.matrix(masses, x, y, z, vx, vy, vz, softening, potential)

	df = pd.DataFrame(data={"m":masses, 'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'softening':softening, 'potential':potential})
	return df


star_df = parse_data()

center = np.array([0,0,0])

pos = star_df.loc[:, ['x', 'y', 'z']]


# spherical shell bins
radii = np.apply_along_axis(np.linalg.norm, 1, pos.values)
star_df.insert(0, 'dist_to_center', radii)

n_bins = 4000
max_r = star_df.dist_to_center.max()
space = np.linspace(0, max_r/100, n_bins)
step = max_r/n_bins/100

mid_space = space + step/2


average_mass_per_bin = []
for r in space:
    ms = star_df.m[(star_df.dist_to_center > r) & (star_df.dist_to_center < r+step)]
    if len(ms) == 0:
        m = 0
    else:
        m = np.average(ms)
    average_mass_per_bin.append(m)

average_mass_per_bin = np.array(average_mass_per_bin)

print(len(radii[radii < 1]))
bins = [len(radii[(r < radii) & (radii <= r+step)]) for r in space]
volumes = [4/3*np.pi*((r+step)**3 - r**3) for r in space]
errors = np.sqrt(bins)



rho = np.array(bins)/np.array(volumes) * average_mass_per_bin
errors = errors/np.array(volumes) * average_mass_per_bin

print(rho[:30])

M = np.sum(star_df.m)

def hernquist_density(r, a):
    return M/2*np.pi * a/(r*(r+a)**3)

from scipy.optimize import curve_fit

popt, pcov = curve_fit(hernquist_density, mid_space, rho)
a = popt[0]
print("a = ", a)


high_res_space = np.linspace(step/2, 0.09, 300)


plt.yscale('log')
plt.errorbar(mid_space[:50], rho[:50], yerr=errors[:50], fmt=".")
plt.plot(high_res_space, hernquist_density(high_res_space, a))
plt.show()








## animation
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = domain, []
ln, = plt.plot([], [], "r.")

def init():
    ax.set_xlim(0, 1)
    ax.set_ylim(-1, 3)
    return ln,

def update(frame):
    ydata = f[frame]
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=range(0, n_steps),
                    init_func=init, blit=True)
plt.show()