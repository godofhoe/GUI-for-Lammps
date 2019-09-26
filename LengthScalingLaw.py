# Check length scaling law (Lognormal)
import numpy as np
import matplotlib.pyplot as plt

f = open('L92total.txt','r')
data = []
for line in f:
    data.append([float(cell) for cell in line.split()])
    
f.close()

L = [row[0] for row in data]

h = np.histogram(L)
b = h[1]
#b = np.arange(10, 110, 10)

plt.hist(L, b)
plt.show()
