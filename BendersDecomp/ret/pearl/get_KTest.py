import numpy as np
import os
import matplotlib.pyplot as plt




path_model_outs = '../model_outs/'
print("Path to model results exists? ",os.path.exists(path_model_outs))

# print([name for name in os.listdir(path_model_outs) if os.path.isfile(os.path.join(path_model_outs, name))])




nb_targets = [8, 8]
algos = [3]
levels = ['e']


for n in range(nb_targets[0], nb_targets[1]+1):
    for l in levels:
        fig, ax  = plt.subplots()
        plt.xticks(np.arange(1, 51, 1))
        for i in range(1,2):
            instname = 'n_'+str(n)+ '_' +l + '_' + str(i)
            OBJ = []
            TIME = []
            BDCUTS = []
            STCUTS = []
            NODES = []
            for k in range(8,51,1):
                out_algo = instname + '_algo_3' + '.KTest_'+ str(k)
                # print(out_algo)
                if os.path.isfile(os.path.join(path_model_outs, out_algo)):
                # print(instname)
                    fileOO = open(path_model_outs+out_algo, "r")                    
                    name = fileOO.readline().split(":")[1]
                    algidx = int(fileOO.readline().split(":")[1])
                    objval = float(fileOO.readline().split(":")[1])
                    objval = float("{0:.3f}".format(objval))
                    time_1 = float(fileOO.readline().split(":")[1])
                    time_GRB = float(fileOO.readline().split(":")[1])
                    time_GRB = float("{0:.3f}".format(time_GRB))
                    BDcuts = int(fileOO.readline().split(":")[1])
                    STcuts = int(fileOO.readline().split(":")[1])
                    Uscuts = int(fileOO.readline().split(":")[1])
                    status_GRB = int(fileOO.readline().split(":")[1])
                    gap = float(fileOO.readline().split(":")[1])
                    gap = float("{0:.3f}".format(gap))
                    nodes_GRB = int(fileOO.readline().split(":")[1])
                    OBJ.append(objval)
                    TIME.append(time_1)
                    BDCUTS.append(BDcuts)
                    STCUTS.append(STcuts)
                    NODES.append(nodes_GRB)
            property = 'OBJ '
            plt.title('Ten instances with size: ' + ': ' + str(n) + ', level: ' + l)
            plt.plot(range(8,len(OBJ)+8,1), OBJ,'o-')
            plt.xlabel('k')
            plt.ylabel(property)
            # plt.plot(range(1,len(TIME_GRB)+1), TIME_GRB,linewidth=1)
            # plt.plot(range(1,len(BDCUTS)+1), BDCUTS,linewidth=1)
            # plt.plot(range(1,len(STCUTS)+1), STCUTS, 'g',linewidth=3)
            # plt.plot(range(1,len(NODES)+1), NODES, 'g',linewidth=3)

        plt.show()
        # pause()






