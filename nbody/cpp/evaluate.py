import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


df = pd.read_csv('test.csv')

acc = df.loc[:,['ax', 'ay', 'az']]
df['acc'] = -np.linalg.norm(acc, axis=1)


df = df.sort_values(by=['r'])
print(df)


cutoff_r = 10

r_space = np.linspace(0.1,cutoff_r,1000)
def analytical_solution(r_space):
	y = []
	G = 1
	for r in r_space:
		m = np.sum(df[df.r < r].m)
		y.append(-G*m/r**2)

	return np.array(y)


cutmask1 = df.r < cutoff_r
cutmask2 = r_space < cutoff_r


plt.plot(df.r[cutmask1], df.acc[cutmask1], '.')
plt.plot(r_space[cutmask2], analytical_solution(r_space)[cutmask2])
plt.show()
#print(df)
#print(df.ax)