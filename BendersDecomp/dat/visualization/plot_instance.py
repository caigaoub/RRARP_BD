import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import matplotlib
import re
import math
import time
from sys import argv
from operator import itemgetter 


# def plot_instance(instancefile, instsize, count):
def plot_instance(instancefile):
    file_ = open(instancefile, 'r')
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)
    nb_tars_ = int(list_[0])
    fig, ax  = plt.subplots()
    # ax.set_title(instancefile)
    minx = 100000
    maxx = -100
    miny = 100000
    maxy = -100
    #departure depot
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    departure = [float(list_[0])]
    departure.append(float(list_[1]))
    minx = min(minx, departure[0]-3)
    maxx = max(maxx, departure[0]+3)
    miny = min(miny, departure[1]-3)
    maxy = max(maxy, departure[1]+3)

    rect = pch.Rectangle((departure[0]-0.1, departure[1]-0.1), 0.2, 0.2,edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(rect)
    #arrival depot
    line_ = file_.readline()
    list_ = re.split("\t|\n", line_)    
    arrival = [float(list_[0])]
    arrival.append(float(list_[1]))
    minx = min(minx, arrival[0]-3)
    maxx = max(maxx, arrival[0]+3)
    miny = min(miny, arrival[1]-3)
    maxy = max(maxy, arrival[1]+3)
    rect = pch.Rectangle((arrival[0]-0.1, arrival[1]-0.1), 0.2, 0.2,edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(rect)
    # circle = pch.Circle((arrival[0], arrival[1]), radius=0.1, edgecolor='r', facecolor='r', alpha=1)
    # ax.add_artist(circle)
    # ax.annotate("("+str(nb_tars_+1)+": "+str(arrival[0])+","+str(arrival[1])+")", xy=(arrival[0], arrival[1]), xytext=(arrival[0]+0.3, arrival[1]+0.3))
    itr = 1
    line_ = file_.readline()
    tar_locs = []
    while itr <= nb_tars_:
        list_ = re.split("\t|\n", line_)
        tar_loc_x = float(list_[0])
        tar_loc_y = float(list_[1])
        tar_locs.append([tar_loc_x, tar_loc_y])
        tar_rad = float(list_[2])
        # print(tar_loc_x,tar_loc_y,tar_rad)
        minx = min(minx, tar_loc_x-tar_rad - 0.2)
        maxx = max(maxx, tar_loc_x+tar_rad + 0.2)
        miny = min(miny, tar_loc_y-tar_rad - 0.2)
        maxy = max(maxy, tar_loc_y+tar_rad + 0.2)
        circle = plt.Circle((tar_loc_x, tar_loc_y), radius=tar_rad, edgecolor='k', facecolor='darkgrey', fill=True,alpha=0.5,lw=1)
        ax.add_artist(circle)
        # circle = plt.Circle((tar_loc_x, tar_loc_y), radius=0.01, edgecolor='r',facecolor='r',alpha=1)
        ax.add_artist(circle)
        ax.annotate(str(itr), xy=(tar_loc_x, tar_loc_y), xytext=(tar_loc_x-0.1, tar_loc_y-0.1),size=11)
        # ax.annotate("("+str(itr)+": "+str(tar_loc_x)+","+str(tar_loc_y)+")", xy=(tar_loc_x, tar_loc_y), xytext=(tar_loc_x+0.3, tar_loc_y+0.3))
        itr += 1   
        line_ = file_.readline()

    if False:
        '''plot inner path'''
        filename = '/home/latte/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/visualization/optpath.txt'
        fileopath_ = open(filename, 'r')
        pathx = []
        pathy = []
        seq = [0]
        for i in range(2*nb_tars_+2):
            line_ = fileopath_.readline()
            if i % 2 == 1:
                seq.append(int(line_.split('\t')[0]))
            pathx.append(float(line_.split('\t')[1]))
            pathy.append(float(line_.split('\t')[2]))
        plt.plot(pathx, pathy,color='b', alpha=1,linewidth=1.3,label="Opt-Sol")
        fileopath_.close()
        filename = '/home/latte/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/visualization/optpath2.txt'
        fileopath_ = open(filename, 'r')
        pathx = []
        pathy = []
        seq = [0]
        for i in range(2*nb_tars_+2):
            line_ = fileopath_.readline()
            if i % 2 == 1:
                seq.append(int(line_.split('\t')[0]))
            pathx.append(float(line_.split('\t')[1]))
            pathy.append(float(line_.split('\t')[2]))
        plt.plot(pathx, pathy,color='k',linestyle='dashed', alpha=1,linewidth=1,label="H-Sol")
        # plt.legend(bbox_to_anchor=(0.76, 0.98), loc='upper left', borderaxespad=0.,prop={'size': 9})
        fileopath_.close()

    file_.close()
    ax.set(xlim=(minx,maxx), ylim=(miny,maxy))
    ax.set_aspect('equal', adjustable='box')
    # plt.xlabel('x', fontsize=10)
    # plt.ylabel('y', fontsize=10)
    # plt.xticks(np.arange(minx,maxx,1))
    # plt.yticks(np.arange(miny,maxy,1))
    # ax.axis('off')
    # plt.grid(alpha=.3)
    # plt.savefig('{}.pdf'.format('/home/latte/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/visualization/example1'), bpi=500)
    plt.show()




if __name__ == "__main__":
    instancefile = argv[1]
    plot_instance(instancefile)

    # instsize = int(argv[2])
    # minidx = int(argv[3])
    # maxidx = int(argv[4])
    # for i in range(minidx, maxidx+1):
        # plot_instance(instancefile, instsize, i)
