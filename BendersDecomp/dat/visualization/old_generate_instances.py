import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from decimal import *


def generate_instance(nb_targets, radius, densitycoef):
    precision = 1000.0
    # densitycoef = 0.5
    area = (densitycoef**2) * nb_targets * radius**2 * math.pi
    width = math.sqrt(area)
    height = width
    tar_locs = []
    rewards_ratio = []
    #depot 1 and 2
    while True:
        x1 =  round(random.uniform(0, width)*precision)/precision
        y1 =  round(random.uniform(0, height)*precision)/precision
        x2 =  round(random.uniform(0, width)*precision)/precision
        y2 =  round(random.uniform(0, height)*precision)/precision
        if x1 != x2 or y1 != y2:
            tar_locs.append((x1,y1,0))
            tar_locs.append((x2,y2,0))
            break
    # target location 
    tidx = 1
    itr = 0
    while True:
        x =  round(random.uniform(0, width)*precision)/precision
        y =  round(random.uniform(0, height)*precision)/precision
        # print(x, y)
        itr += 1
        # if itr > 10000:
        #     print('endless loop')
        #     exit()
        covered = False
        for e in tar_locs:
            dist = math.sqrt((e[0]-x)**2 + (e[1]-y)**2) - e[2] - radius - 0.5
            if dist < 0:
                covered = True
                break
        if covered == False:
            tar_locs.append((x,y,radius))
            tidx += 1
        if tidx > nb_targets:
            break
    # print(tar_locs)
    rewards_ratio = []
    for i in range(0, nb_targets):
        rewards_ratio.append(round(random.uniform(0.7, 0.9)*precision)/precision)
    return tar_locs, rewards_ratio, width, height


def plot_instance(tar_locs, width, height):
    fig, ax  = plt.subplots()
    circle = pch.Circle((tar_locs[0][0], tar_locs[0][1]), radius=0.05, edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(circle)
    ax.annotate("("+str(0)+": " + str(tar_locs[0][0])+","+str(tar_locs[0][1])+")", xy=(tar_locs[0][0], tar_locs[0][1]), xytext=(tar_locs[0][0], tar_locs[0][1]))

    circle = pch.Circle((tar_locs[1][0], tar_locs[1][1]), radius=0.05, edgecolor='r', facecolor='r', alpha=1)
    ax.add_artist(circle)
    ax.annotate("("+str(1)+": " + str(tar_locs[1][0])+","+str(tar_locs[1][1])+")", xy=(tar_locs[1][0], tar_locs[1][1]), xytext=(tar_locs[1][0], tar_locs[1][1]))
    
    for i in range(2,len(tar_locs)):
        circle = plt.Circle((tar_locs[i][0], tar_locs[i][1]), radius=tar_locs[i][2], edgecolor='k',facecolor='k',alpha=0.2)
        ax.add_artist(circle)
        circle = plt.Circle((tar_locs[i][0], tar_locs[i][1]), radius=0.01, edgecolor='r',facecolor='r',alpha=0.8)
        ax.add_artist(circle)
        ax.annotate("("+str(i-1)+": "+str(tar_locs[i][0])+","+str(tar_locs[i][1])+")", xy=(tar_locs[i][0], tar_locs[i][1]), xytext=(tar_locs[i][0], tar_locs[i][1]))
    

    ax.set(xlim=(-1,width+1), ylim=(-1,height+1))
    ax.set_aspect('equal', adjustable='box')
    plt.show()
def write_instance(tar_locs, rewards_ratio, filename):
    file = open(filename, "w")
    file.write(str(len(tar_locs)-2) + '\n')
    file.write(str(tar_locs[0][0]) + '\t' + str(tar_locs[0][1]) + '\n')
    file.write(str(tar_locs[1][0]) + '\t' + str(tar_locs[1][1]) + '\n')
    for i in range(2,len(tar_locs)):
        file.write(str(tar_locs[i][0]) + '\t' + str(tar_locs[i][1]) + '\t' + str(tar_locs[i][2]) + '\n')
    # print(rewards_ratio)
    for i in range(0,len(rewards_ratio)-1):
        file.write(str(rewards_ratio[i]) + '\t')
    file.write(str(rewards_ratio[-1]))
    file.close()
if __name__ == "__main__":
    ''' basic setting '''
    radius = 1.0
    size_range = [6, 20]
    
    ''' instance generation '''
    # densitycoefs = [4, 2, 1.5]
    # for nb_targets in np.arange(size_range[0],size_range[1]+1):
    # tar_locs, rewards_ratio, width, height = generate_instance(16, radius, 0.3)
    # plot_instance(tar_locs, width, height)
    # write_instance(tar_locs, rewards_ratio, '../n_16_m_1.dat')
   path = '/home/latte/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/dat/'
    # path = '/projects/academic/josewalt/caigao/RRARP_BD/BendersDecomp/dat/'
    for nb_targets in range(size_range[0], size_range[1]+1):
        densitycoef = 7
        print("currently generate instance with size ", nb_targets)
        for i in range(1,11):
            print('e',i)
            tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
            write_instance(tar_locs, rewards_ratio, path +'n_'+ str(nb_targets) +'_' + 'e_' + str(i) + '.dat')

        densitycoef = 6
        for i in range(1,11):
            print('m',i)

            tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
            write_instance(tar_locs, rewards_ratio,  path +'n_' + str(nb_targets) +'_' + 'm_' + str(i) + '.dat')

        densitycoef = 5
        for i in range(1,11):
            print('h',i)
            tar_locs, rewards_ratio, width, height = generate_instance(nb_targets, radius, densitycoef)
            write_instance(tar_locs, rewards_ratio,  path +'n_'+ str(nb_targets) +'_' + 'h_' + str(i) + '.dat')
