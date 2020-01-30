#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sys import argv

import matplotlib.pyplot as plt

if len(argv) != 2:
	print("Usage: {} csvfile.csv".format(argv[0]))
	exit()

_, csvfilename = argv
df = pd.read_csv(csvfilename)

acc = df.loc[:,['ax', 'ay', 'az']]
df['acc'] = np.linalg.norm(acc, axis=1)


df = df.sort_values(by=['r'])
print(df)

# half mass stuff

N = df.shape[0]
mid = int(df.shape[0]/2)
R_hm = df.loc[mid].r

M_hm = np.sum(df[df.r < R_hm].m)

t_cross = np.sqrt(R_hm**3/M_hm)
t_rlx = N/(8*np.log(N)) * np.sqrt(R_hm**3/M_hm)

print("t_cross", t_cross)
print("t_rlx", t_rlx)




# ernquist stuff

cutoff_r = 5

r_space = np.linspace(0.1,cutoff_r,1000)
def analytical_solution(r_space):
	# ernquist density
	y = []
	G = 1
	for r in r_space:
		m = np.sum(df[df.r < r].m)
		y.append(G*m/r**2)

	return np.array(y)


cutmask1 = df.r < cutoff_r
cutmask2 = r_space < cutoff_r


plt.plot(df.r[cutmask1], df.acc[cutmask1], '.', label='numerical solution')
plt.plot(r_space[cutmask2], analytical_solution(r_space)[cutmask2], label='analytical solution (erquist)')
plt.xlabel("Radius")
plt.ylabel("Acceleration |a|")
plt.legend()
plt.show()
#print(df)
#print(df.ax)



