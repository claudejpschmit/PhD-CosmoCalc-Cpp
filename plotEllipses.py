from matplotlib import pyplot as plt
from matplotlib import patches as ptc
from matplotlib.ticker import MaxNLocator
import sys
from math import *

filename = sys.argv[1]
param_filename = sys.argv[2]
with open(filename) as f:
    content = [float(x) for x in f if x != '\n']
with open(param_filename) as f:
    params = [line.rstrip('\n') for line in f]

num_params = int(content[0])
ellipse_number = 0
fig, ax = plt.subplots()
for i in range(0,num_params - 1):
    for j in range(i+1, num_params):

        w = 2*sqrt(content[ellipse_number*5 + 1])
        h = 2*sqrt(content[ellipse_number*5 + 2])
        theta = content[ellipse_number*5 + 3]
        x = content[ellipse_number*5 + 4]
        y = content[ellipse_number*5 + 5]
        frame_index = (num_params-1)*i + j
        ax1 = plt.subplot(num_params-1,num_params-1, frame_index)
        plt.subplots_adjust(hspace = 0.001, wspace = 0.001)
        if (j-i) > 1:
            ax1.xaxis.set_ticklabels([])
            ax1.yaxis.set_ticklabels([])
        if (j-i) == 1:
            plt.ylabel(params[i],rotation='horizontal') 
            ax1.yaxis.labelpad = 50
            labels = ax1.yaxis.get_major_ticks()
            labels[-1].label1.set_visible(False)
            #ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
            if (i == num_params - 2):
                plt.xlabel(params[i+1])
                ax1.xaxis.labelpad = 10
        ellipse = ptc.Arc(xy = (x,y), width=w, height=h, angle=theta)
        ellipse2 = ptc.Ellipse(xy = (x,y), width=w, height=h, angle=theta)
        ax1.add_artist(ellipse)
        ax1.add_artist(ellipse2)
       
        m = max(w,h)
        limx = max(abs(w*cos(theta)), abs(h*sin(theta)))
        limy = max(abs(w*sin(theta)), abs(h*cos(theta)))
        #print("ellipse ", ellipse_number, w,h, limx, limy , theta)
        #limx = max(abs(cos(theta)*w), abs(cos(theta)*h))
        #limy = max(abs(sin(theta)*h), abs(sin(theta)*w))
        ax1.set_xlim([ x - limx, x + limx])
        ax1.set_ylim([ y - limy, y + limy])

        #p_err = 1 
        #ax1.set_xlim([(1-p_err) * x, (1+p_err) * x])
        #ax1.set_ylim([(1-p_err) * y, (1 + p_err) * y])

        ellipse_number += 1

plt.show()
