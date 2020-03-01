import numpy as np
import os





path_model_outs = '../model_outs/'
print("Path to model results exists? ",os.path.exists(path_model_outs))

# print([name for name in os.listdir(path_model_outs) if os.path.isfile(os.path.join(path_model_outs, name))])




nb_targets = [6, 35]
algos = [2,3,4]
levels = ['e','m','h']


fileII = open("results.txt", "w")
fileII.write("\t\t\t\t\t\t"+"algo_2"+"\t\t\t\t\t\t\t\t\t\t\t\t\t\t"+"algo_3"+"\t\t\t\t\t\t\t\t\t\t\t\t\t\t"+"algo_4"+'\n')
fileII.write('{:10s}'.format("Name")+'\t'+'{:7s}'.format("Obj")+'\t'+'{:7s}'.format("Time")+'\t'+'{:7s}'.format("BDCs")+'\t'+ '{:7s}'.format("STCs") +'\t'+'{:7s}'.format("Gap")+'\t'+'{:7s}'.format("Nodes")+'\t' + '{:7s}'.format("~~"))
fileII.write('\t'+'{:7s}'.format("Obj")+'\t'+'{:7s}'.format("Time")+'\t'+'{:7s}'.format("BDCs")+'\t'+ '{:7s}'.format("STCs") +'\t'+'{:7s}'.format("Gap")+'\t'+'{:7s}'.format("Nodes")+'\t' + '{:7s}'.format("~~"))
fileII.write('\t'+'{:7s}'.format("Obj")+'\t'+'{:7s}'.format("Time")+'\t'+'{:7s}'.format("BDCs")+'\t'+ '{:7s}'.format("STCs") +'\t'+'{:7s}'.format("Gap")+'\t'+'{:7s}'.format("Nodes")+'\t' + '{:7s}'.format("~~"))
fileII.write('\n')
for n in range(nb_targets[0], nb_targets[1]+1):
	for l in levels:
		for i in range(1,11):
			instname = 'n_'+str(n)+ '_' +l + '_' + str(i)
			fileII.write('{:10s}'.format(instname) + '\t')
			for al in algos:
				out_algo = instname + '_algo_' + str(al) + '.out'
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
					gap = float(fileOO.readline().split(":")[1]) * 100.0
					gap = float("{0:.2f}".format(gap))
					nodes_GRB = fileOO.readline().split(":")[1].split("\n")[0]
					# print(name, algidx, objval, time_GRB, BDcuts, STcuts, Uscuts, status_GRB, gap, nodes_GRB)
					fileII.write('{:7s}'.format(str(objval))+'\t'+ '{:7s}'.format(str(time_GRB))+'\t' + '{:7s}'.format(str(BDcuts))+'\t' + '{:7s}'.format(str(STcuts)) + '\t' + '{:7s}'.format(str(gap)) +'\t' + '{:7s}'.format(nodes_GRB)+'\t'+'{:7s}'.format("~~"))				
					fileOO.close()
				else:
					fileII.write('{:7s}'.format("-")+'\t'+'{:7s}'.format("-")+'\t'+'{:7s}'.format("-")+'\t'+'{:7s}'.format("-")+'\t'+'{:7s}'.format("-")+'\t'+'{:7s}'.format("-")+'\t' +'{:7s}'.format("~~"))
					
			fileII.write('\n')

fileII.close()







