import matplotlib.pyplot as plt
import numpy as np
import sys

l = sys.argv[1]
num = sys.argv[2]

with open("movie_Cl/frame.dat") as f:
    data = f.read()
data = data.split('\n')
del data[-1]
x = [row.split(' ')[0] for row in data]
y1 = [row.split(' ')[1] for row in data]
plot1 = plt.plot(x,y1, label = r'simplified Cl')

plt.legend()
plt.title("l = " + l)
axes = plt.gca()
axes.set_ylim([-10e7,10e7])
plt.savefig("movie_Cl/frame_" + num)
