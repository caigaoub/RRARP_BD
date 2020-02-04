import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import re
import math
import time

def set_plot_attributes(plt, ax):
    plt.grid(alpha=.5)
    # plt.xticks(np.arange(-10,1000,1))
    # plt.yticks(np.arange(-10,1000,1))
    ax.set(xlim=(0,1000), ylim=(0,1000))
    # ax.set_aspect('equal', adjustable='box')


def plot_instance(instancefile):
    file_ = open(instancefile, 'r')
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)
    nb_tars_ = int(list_[0])
    fig, ax  = plt.subplots()
    # set_plot_attributes(plt, ax)
    ax.set_title(instancefile)
    minx = 100000
    maxx = -100
    miny = 100000
    maxy = -100
    #departure
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    departure = [float(list_[0])]
    departure.append(float(list_[1]))
    minx = min(minx, departure[0])
    maxx = max(maxx, departure[0])
    miny = min(miny, departure[1])
    maxy = max(maxy, departure[1])
    rect = pch.Rectangle((departure[0], departure[1]), 5, 5, edgecolor='r', facecolor='r', alpha=1)
    ax.add_patch(rect)
    #arrival
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    arrival = [float(list_[0])]
    arrival.append(float(list_[1]))
    minx = min(minx, arrival[0])
    maxx = max(maxx, arrival[0])
    miny = min(miny, arrival[1])
    maxy = max(maxy, arrival[1])
    rect = pch.Rectangle((arrival[0], arrival[1]), 5, 5, edgecolor='r', facecolor='r', alpha=1)
    ax.add_patch(rect)
    #target circles
    itr = 1
    line_ = file_.readline()
    while itr <= nb_tars_:
        list_ = re.split("\t|\n", line_)
        tar_loc_x = float(list_[0])
        tar_loc_y = float(list_[1])
        tar_rad = float(list_[2])
        # print(tar_loc_x,tar_loc_y,tar_rad)
        minx = min(minx, tar_loc_x-tar_rad)
        maxx = max(maxx, tar_loc_x+tar_rad)
        miny = min(miny, tar_loc_y-tar_rad)
        maxy = max(maxy, tar_loc_y+tar_rad)
        circle = plt.Circle((tar_loc_x, tar_loc_y), radius=tar_rad, edgecolor='k',facecolor='k',alpha=0.6)
        ax.add_artist(circle)
        itr += 1   
        line_ = file_.readline()
    file_.close()
    ax.set(xlim=(minx,maxx), ylim=(miny,maxy))
    ax.set_aspect('equal', adjustable='box')

    plt.show()

################################################################################################

#  transformation files
# instancefile = '../ret/inst_n_20/n_20_h_1.txt'
instancefile = '../ret/cluster_n_30/n_30_c_5_1.txt'

plot_instance(instancefile)