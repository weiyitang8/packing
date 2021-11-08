import os, sys
import pandas as pd


# input_dir= sys.argv[1]
# input_dir = '/Users/tangweiyi/desktop/design/output'
input_dir = '/home/sting-gpu/weiyi/design/setup5_output2/good_packing'
dir_files = os.listdir(input_dir)
pdb_files = [x for x in dir_files if x[-3:]=='pdb']

data = list(range(len(pdb_files)))


for i, pdb in enumerate(pdb_files):
	path = os.path.join(input_dir,pdb)

	inF = open(path,'r')
	pdb_text = inF.readlines()

	for l in pdb_text[-500:]:
		if len(l) <1: continue
		line= l.rsplit()
		if len(line) <2:continue

		if 'holes-ABX_interface' == line[0]:
			holes_interface 		= float( line[1] )
			data[i]=holes_interface
df= pd.DataFrame(data, columns=['holes_interface'],index=pdb_files)
df= df.sort_values(by='holes_interface')
print(df)
			

