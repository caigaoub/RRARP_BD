import numpy as np
import os
import matplotlib.pyplot as plt


path_model_outs = '/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/model_outs/FischettiTest/'
print("Path to model results exists? ", os.path.exists(path_model_outs))


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
IdxSet = [0, 1, 2]
MTC = {} #'nli': Metric
for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            instname = 'n_'+str(n)+ '_' +l + '_' + str(i)
            for idx in IdxSet:
                x = Metric()
                out_algo = instname + '_algo_' + str(3) + '.FisTest_'+str(idx)
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
                    gap = float(fileOO.readline().split(":")[1])
                    gap = float("{0:.3f}".format(gap))
                    nodes_GRB = int(fileOO.readline().split(":")[1])
                    x._OBJs = objval
                    x._TIMEs = time_1
                    x._BDCUTs = BDcuts
                    x._STCUTs = STcuts
                    x._GAPs = gap
                    x._NODEs = nodes_GRB
                MTC[str(n)+l+str(i)+str(idx)] = x

for n in nb_targets:
    for l in levels:
        for i in range(1,11):
            for idx in IdxSet:
                print(str(n)+l+str(i)+str(idx), MTC[str(n)+l+str(i)+str(idx)]._OBJs)




fileWWW = open('/home/cai/Dropbox/Box_Research/Github/RRARP_BD/BendersDecomp/ret/pearl/fis-test_txt', "w")
for n in [9]:
    for l in levels:
        for i in range(1,11):
            fileWWW.write('n-'+str(n)+'-'+l+'-'+str(i))
            instname = str(n)+l+str(i)+str(0)
            fileWWW.write(' & ' + str(MTC[instname]._OBJs))
            fileWWW.write(' & ' + str(MTC[instname]._TIMEs))
            fileWWW.write(' & ' + str(MTC[instname]._BDCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._STCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._GAPs))
            fileWWW.write(' & ' + str(MTC[instname]._NODEs))
            instname = str(n)+l+str(i)+str(1)
            fileWWW.write(' & ')
            fileWWW.write(' & ' + str(MTC[instname]._OBJs))
            fileWWW.write(' & ' + str(MTC[instname]._TIMEs))
            fileWWW.write(' & ' + str(MTC[instname]._BDCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._STCUTs))
            fileWWW.write(' & ' + str(MTC[instname]._GAPs))
            fileWWW.write(' & ' + str(MTC[instname]._NODEs))
            # instname = str(n)+l+str(i)+str(2)
            # fileWWW.write(' & ')
            # fileWWW.write(' & ' + str(MTC[instname]._OBJs))
            # fileWWW.write(' & ' + str(MTC[instname]._TIMEs))
            # fileWWW.write(' & ' + str(MTC[instname]._BDCUTs))
            # fileWWW.write(' & ' + str(MTC[instname]._STCUTs))
            # fileWWW.write(' & ' + str(MTC[instname]._GAPs))
            # fileWWW.write(' & ' + str(MTC[instname]._NODEs))
            fileWWW.write('\\\\'+'\n')
        fileWWW.write('\\midrule'+'\n')
fileWWW.close()


