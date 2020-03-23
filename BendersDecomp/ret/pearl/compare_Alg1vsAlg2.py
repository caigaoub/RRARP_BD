import numpy as np
import os
import matplotlib.pyplot as plt



## use dictionary

path_model_outs = '/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/model_outs/Alg1vsAlg2/'
print("Path to model results exists? ", os.path.exists(path_model_outs))

# print([name for name in os.listdir(path_model_outs) if os.path.isfile(os.path.join(path_model_outs, name))])


''' define a class for storing the data '''
class Metric:

    def __init__(self):
        self._name = ''
        self._OBJs = -1
        self._TIMEs = -1
        self._GAPs = -1
        self._BDCUTs = -1
        self._STCUTs = -1
        self._NODEs = -1

''' Basic setting'''
nb_targets = [6, 7, 8, 9]
levels = ['e', 'm', 'h']
algs = [1, 2]
MTC = {} #'nli': Metric
for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            instname = 'n_'+str(n)+ '_' +l + '_' + str(i)
            for al in algs:
                x = Metric()
                out_algo = instname + '_algo_' + str(al) + '.out'
                x._name = out_algo
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
                    gap = float(fileOO.readline().split(":")[1]) * 100.0
                    gap = float("{0:.3f}".format(gap))
                    nodes_GRB = int(fileOO.readline().split(":")[1])
                    x._OBJs = objval
                    x._TIMEs = time_1
                    x._BDCUTs = BDcuts
                    x._STCUTs = STcuts
                    x._GAPs = gap
                    x._NODEs = nodes_GRB
                MTC[str(n)+l+str(i)+str(al)] = x

for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            for al in algs:
                print(str(n)+l+str(i)+str(al), MTC[str(n)+l+str(i)+str(al)]._OBJs)




fileWWW = open('/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/pearl/alg1_alg2_instance_8.txt', "w")
for n in [8]:
    for l in levels:
        for i in range(1,11):
            fileWWW.write('n-'+str(n)+'-'+l+'-'+str(i))
            instname = str(n)+l+str(i)+str(1)
            fileWWW.write(' & ' + str(MTC[instname]._OBJs))
            fileWWW.write(' & ' + str(MTC[instname]._TIMEs))
            fileWWW.write(' & ' + str(MTC[instname]._BDCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._STCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._GAPs))
            fileWWW.write(' & ' + str(MTC[instname]._NODEs))
            instname = str(n)+l+str(i)+str(2)
            fileWWW.write(' & ')
            fileWWW.write(' & ' + str(MTC[instname]._OBJs))
            fileWWW.write(' & ' + str(MTC[instname]._TIMEs))
            fileWWW.write(' & ' + str(MTC[instname]._BDCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._STCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._GAPs))
            fileWWW.write(' & ' + str(MTC[instname]._NODEs))
            fileWWW.write('\\\\'+'\n')
        fileWWW.write('\\midrule'+'\n')
fileWWW.close()



# print('--------------------------------------------------------------')

''' take average '''
AVG_MTC = {}
for n in nb_targets:
    for l in levels:
        for al in algs:
            y = Metric()
            y._name = str(n)+'_'+l
            avg_obj = 0.0
            avg_time = 0.0
            avg_BDcuts = 0.0
            avg_STcuts = 0.0
            avg_gap = 0.0
            avg_nodes = 0.0
            for i in range(1,11):
                out_algo = str(n)+l+str(i)+str(al)
                avg_obj += MTC[out_algo]._OBJs/10.0
                avg_time += MTC[out_algo]._TIMEs/10.0
                avg_BDcuts += float(MTC[out_algo]._BDCUTs)/10.0
                avg_STcuts += float(MTC[out_algo]._STCUTs)/10.0
                avg_gap += MTC[out_algo]._GAPs/10.0
                avg_nodes += float(MTC[out_algo]._NODEs)/10.0
            y._OBJs = float("{0:.3f}".format(avg_obj))
            y._TIMEs = float("{0:.3f}".format(avg_time))
            y._BDCUTs = float("{0:.3f}".format(avg_BDcuts))
            y._STCUTs = float("{0:.3f}".format(avg_STcuts))
            y._GAPs = float("{0:.3f}".format(avg_gap*100.0))
            y._NODEs = float("{0:.3f}".format(avg_nodes))
            AVG_MTC[str(n)+l+str(al)] = y



# fileWWW = open('/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/pearl/alg1_alg2.txt', "w")
# for n in nb_targets:
#     for l in levels:
#         fileWWW.write('n-'+str(n)+'-'+l)
#         instname = str(n)+l+str(1)
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._OBJs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._TIMEs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._BDCUTs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._STCUTs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._GAPs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._NODEs))
#         instname = str(n)+l+str(2)
#         fileWWW.write(' & ')
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._OBJs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._TIMEs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._BDCUTs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._STCUTs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._GAPs))
#         fileWWW.write(' & ' + str(AVG_MTC[instname]._NODEs))
#         fileWWW.write('\\\\'+'\n')
#     fileWWW.write('\\midrule'+'\n')
# fileWWW.close()

