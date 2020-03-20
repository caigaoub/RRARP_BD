import numpy as np
import os
import matplotlib.pyplot as plt



## use dictionary

path_model_outs = '../model_outs/KTest/'
print("Path to model results exists? ", os.path.exists(path_model_outs))

# print([name for name in os.listdir(path_model_outs) if os.path.isfile(os.path.join(path_model_outs, name))])


''' define a class for storing the data '''
class Metric:

    def __init__(self):
        self._name = ''
        self._OBJs = {}
        self._TIMEs = {}
        self._GAPs = {}
        self._BDCUTs = {}
        self._STCUTs = {}
        self._NODEs = {}

''' Basic setting'''
nb_targets = [6, 7, 8]
levels = ['e', 'm', 'h']
K = range(8,21,1)

MTC = {} #'nli': Metric
for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            instname = 'n_'+str(n)+ '_' +l + '_' + str(i)
            x = Metric()
            x._name = instname
            for k in K:
                out_algo = instname + '_algo_2' + '.KTest_'+ str(k)
                # print(out_algo)
                if os.path.isfile(os.path.join(path_model_outs, out_algo)):
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
                    x._OBJs[k] = objval
                    x._TIMEs[k] = time_1
                    x._BDCUTs[k] = BDcuts
                    x._STCUTs[k] = STcuts
                    x._GAPs[k] = gap
                    x._NODEs[k] = nodes_GRB
            MTC[str(n)+l+str(i)] = x
            # print(str(n)+l+str(i), MTC[str(n)+l+str(i)]._OBJs)

for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            print(str(n)+l+str(i), MTC[str(n)+l+str(i)]._OBJs)


print('--------------------------------------------------------------')

''' take average '''
AVG_MTC = {}
for n in nb_targets:
    for l in levels:
        y = Metric()
        y._name = str(n)+'_'+l
        for k in K:
            avg_obj = 0.0
            avg_time = 0.0
            avg_BDcuts = 0.0
            avg_STcuts = 0.0
            avg_gap = 0.0
            avg_nodes = 0.0
            for i in range(1,11):
                instname = str(n)+l+str(i)
                print(instname, MTC[instname]._OBJs)
                avg_obj += MTC[instname]._OBJs[k]/10.0
                avg_time += MTC[instname]._TIMEs[k]/10.0
                avg_BDcuts += float(MTC[instname]._BDCUTs[k])/10.0
                avg_STcuts += float(MTC[instname]._STCUTs[k])/10.0
                avg_gap += MTC[instname]._GAPs[k]/10.0
                avg_nodes += float(MTC[instname]._NODEs[k])/10.0
            # print(avg_obj)
            y._OBJs[k] = float("{0:.3f}".format(avg_obj))
            y._TIMEs[k] = float("{0:.3f}".format(avg_time))
            y._BDCUTs[k] = float("{0:.3f}".format(avg_BDcuts))
            y._STCUTs[k] = float("{0:.3f}".format(avg_STcuts))
            y._GAPs[k] = float("{0:.3f}".format(avg_gap))
            y._NODEs[k] = float("{0:.3f}".format(avg_nodes))
        AVG_MTC[str(n)+l] = y
        # print(AVG_MTC[str(n)+l]._OBJs)


subsetK = [8,12,16,20]
fileWWW = open('k-test.txt', "w")
for n in nb_targets:
    for l in levels:
        # print('n-'+str(n)+'-'+l, '&', (e for e in list(AVG_MTC[str(n)+l]._OBJs.values())[::4]))
        fileWWW.write('n-'+str(n)+'-'+l + ' & ')
        for k in subsetK:
            fileWWW.write(str(AVG_MTC[str(n)+l]._OBJs[k])+ ' & ')
        fileWWW.write(' & ')
        for k in subsetK:
            fileWWW.write(str(AVG_MTC[str(n)+l]._TIMEs[k])+ ' & ')
        fileWWW.write(' & ')
        for k in subsetK:
            fileWWW.write(str(AVG_MTC[str(n)+l]._BDCUTs[k])+ ' & ')
        fileWWW.write(' & ')

        for k in subsetK:
            fileWWW.write(str(AVG_MTC[str(n)+l]._STCUTs[k])+ ' & ')
        fileWWW.write(' & ')

        for k in subsetK:
            fileWWW.write(str(AVG_MTC[str(n)+l]._NODEs[k])+ ' & ')
        fileWWW.write(' & ')
        fileWWW.write('\n')
fileWWW.close()







        # fig, ax  = plt.subplots()

        # plt.xticks(np.arange(1, 51, 1))

            # property = 'OBJ '
            # plt.title('Ten instances with size: ' + ': ' + str(n) + ', level: ' + l)
            # plt.plot(range(8,len(OBJ)+8,1), OBJ,'o-')
            # plt.xlabel('k')
            # plt.ylabel(property)
            # plt.plot(range(1,len(TIME_GRB)+1), TIME_GRB,linewidth=1)
            # plt.plot(range(1,len(BDCUTS)+1), BDCUTS,linewidth=1)
            # plt.plot(range(1,len(STCUTS)+1), STCUTS, 'g',linewidth=3)
            # plt.plot(range(1,len(NODES)+1), NODES, 'g',linewidth=3)

        # plt.show()
        # pause()
