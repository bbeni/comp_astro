#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


dx = 0.01
v0 = 1

#CFL condition for Delta t
dt = 0.5 * dx/v0


domain = np.linspace(0,1,100)
#f0 = np.heaviside(domain - 0.5, 0.5)
f0 = 1 + np.exp(-(domain-0.5)**2/0.1**2)

def step(u):
	N = len(u)
	vp = np.max(v0, 0) # velocity +
	vm = np.min(v0, 0) # velocity -

	u_next = []

	for i in range(0, N):
		uxm = (u[i] - u[i-1])/dx
		uxp = (u[(i+1) % N ] - u[i])/dx
		u_next.append(u[i] - dt*(vp*uxm + vm*uxp))

	return np.array(u_next)


## integrate 100 timesteps
n_steps = 100
f = [f0]
for i in range(n_steps):
	f.append(step(f[-1]))



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