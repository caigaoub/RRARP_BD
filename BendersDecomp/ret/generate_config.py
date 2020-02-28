import numpy as np




"""  generate configuration files """

# algo_idx = 4
# nb_dstzn = 8




instsize_min = 6
instsize_max = 35
idx_config = 1
diff_levels = ['e','m','h']
subidx_min = 1 
subidx_max = 10
for instsize in range(instsize_min, instsize_max+1):
	for level in diff_levels:
		for sidx in range(subidx_min, subidx_max+1):
			print(idx_config)
			file = open("./configs/config_" + str(idx_config), "w")
			file.write("n_"+str(instsize) + '_' + level + '_' + str(sidx)+".dat")
			file.close()
			idx_config += 1
	#		print(idx_config)
	print("done with instance with", instsize, "targets")
