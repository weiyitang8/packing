import os,sys,shutil
from collections import Counter, defaultdict
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pylab as plt
import numpy as np



input_dir= sys.argv[1]
result_w_detail=open(os.path.join(input_dir,'result_w_detail.txt'))
result=result_w_detail.read()
design_set=result.split('##')

def zero():
	return 0
# ojbs=pd.DataFrame()
result={}
for design in design_set:
	resultdict=defaultdict(zero)
	l=design.rsplit('\n')
	# print(l)
	# sys.exit()
	if len(l) <4:continue
	# if l
	for line in l:
		if len(line)<2:continue
		if line[0]=='*':continue

		if line[0]=='#':
			index=line[-12:-4]
			# print(index)
			continue
		else:
			pair=line.rsplit('  [')
			# print(pair[0])
			pair[1]= pair[1].replace(']','')
			knobs=pair[1].rsplit(',')
			# print(knobs)

			if len(knobs)>1:
				for i in range(len(knobs)):
					new_pair=pair[0]+knobs[i].replace(' ','')
					new_pair=new_pair.replace("'k",' k')
					new_pair=new_pair.replace("'",'')
					# print(new_pair)
					resultdict[new_pair]=1
			
			else:
				new_pair=pair[0]+knobs[0]
				new_pair=new_pair.replace("'k",' k')
				new_pair=new_pair.replace("'",'')

				resultdict[new_pair]=1
				# print(new_pair)
	result[index]=resultdict
	# resultdict=pd.Series(resultdict,columns=index)
	# ojbs=pd.concat([ojbs,resultdict],axis=1)
	# print(ojbs)
ojbs=pd.DataFrame(result)
ojbs=ojbs.fillna(0)
ojbs=ojbs.sort_index()
ojbs=ojbs.T
ojbs=ojbs.sort_index()
# ojbs=ojbs.sort_columns()

print(ojbs)

disMat=sch.distance.pdist(ojbs,metric='euclidean')
print(disMat)
print(len(disMat))
print(disMat.mean())
Z=sch.linkage(disMat,method='complete',metric='distance')


P=sch.dendrogram(Z,labels=ojbs.index)
plt.show()

h_clust=sch.fcluster(Z,t=1.8,criterion='distance')
# print(len(set(h_clust)))
print(h_clust)
design_names_sorted=ojbs.index.to_list()
# print(design_names_sorted)
design_clusters=defaultdict(list)
# design_cluster_key=defaultdict(list)

for n in np.arange(len(design_names_sorted)):
	design_clusters[h_clust[n]].append(design_names_sorted[n])
	# design_cluster_key[h_clust[n]].append(design_names_sorted[n])

for k,v in sorted(design_clusters.items()):
	print('this is cluster Number', k)
	for p in v:
		print(result[p])



singletons_dir=os.path.join(input_dir,'singletons')
if not os.path.exists(singletons_dir):
	os.mkdir(singletons_dir)

for cluster_index,v in sorted(design_clusters.items()):
	# print(cluster_index,v)
	cluster_size=len(v)
	if cluster_size>1:
		clust_dir=os.path.join(input_dir,'cluster_%s_%d')%(cluster_index,cluster_size)
		if not os.path.exists(clust_dir):
			os.mkdir(clust_dir)
		for design in v:
			design_path = os.path.join(input_dir,'design_input_struct_%s.pdb')%design
			shutil.copy(design_path,clust_dir)
	else:
		design_path = os.path.join(input_dir,'design_input_struct_%s.pdb')%v[0]
		shutil.copy(design_path,singletons_dir)






	# print(resultdict)
			# print(pair[1])
			# sys.exit()
			# pair[1]=list(pair[1])


			# print('----')
			# print(pair)




# result_line=result_w_detail.readlines()
# score_sum=[]
# packcounter=Counter(result_line)

# top30=packcounter.most_common(30)
# for pack in top30:
# 	print(pack)
# del packcounter['\n']
# del packcounter['###    design_input_struct_030_0021.pdb\n']
# temp=packcounter.items()
# for pack in temp:
# 	if len(pack[0]) <3: 
# 		del pack
# 	elif '#' in pack[0]: 
# 		del pack
# 	elif '*' in pack[0]:
# 		score_sum.append(pack)
# 		del pack

