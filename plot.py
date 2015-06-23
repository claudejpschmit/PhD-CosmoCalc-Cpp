import matplotlib.pyplot as plt
import numpy as np
import sys

k1 = sys.argv[1]
k2 = sys.argv[2]
num = sys.argv[3]

with open("movie/frame.dat") as f:
    data = f.read()
data = data.split('\n')
del data[-1]
x = [row.split(' ')[0] for row in data]
y1 = [row.split(' ')[1] for row in data]
y2 = [row.split(' ')[2] for row in data]
plot1 = plt.plot(x,y1, label = r'long')
plot2 = plt.plot(x,y2, label = r'simple')

plt.legend()
plt.title("k1 = " + k1 + ", k2 = " + k2)

plt.savefig("movie/frame_" + num)
