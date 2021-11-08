import os, sys
import numpy as np, pandas as pd
import Bio.PDB 
from Bio.PDB.PDBParser import PDBParser
import pack_analysisV4_sumonly
import shutil

AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

input_dir = sys.argv[1]
# output_dir = sys.argv[2]
# input_dir='/Users/tangweiyi/desktop/design/output'
dir_files = os.listdir(input_dir)
pdb_files = [x for x in dir_files if x[-3:]=='pdb']
data       ={}

parser = PDBParser(QUIET=True)

result_sum_score = open(os.path.join(input_dir,'result_sum_score.txt'),'r')
result_line =result_sum_score.readlines()

result_sum_score_blast = open(os.path.join(input_dir,'result_sum_score_blast.txt'),'w')
txt= '\n # Blast of good packing designs\n\n'

for pdb in pdb_files:
	path=os.path.join(input_dir,pdb)
	structure=parser.get_structure('struct',path)
	structure=structure[0]
	chainX=structure['X']
	seq=[]
	
	for res in chainX:
		seq.append(AA[res.get_resname()])
	seq=''.join(seq)

	

	pdb_text=open(path,'r').readlines()
	for l in pdb_text[-500:]:
		if len(l) <1: continue
		line= l.rsplit()
		if len(line) <2:continue
		if 'holes-ABX_interface' == line[0]:
			holes_interface 		= float( line[1] )

	
	for l in result_line:
		if len(l)<1 : continue
		line=l.rsplit()
		if len(line) <2 :continue
		if pdb == line[0]:
			txt+= line[0]+'     '+ line[1]+'     ' + line[2]+ '     '+seq+'    '+'%.3f'%holes_interface+ ' \n'
			data[pdb] = [line[1],line[2],seq,holes_interface]




result_df=pd.DataFrame(data,index=['finalsc','#sumknobs','Seq','DAball']).T
result_df=result_df.sort_values(by=['Seq','DAball','finalsc'])
result_df.to_csv(os.path.join(input_dir,'result_sum_score_blast.csv'))
print(result_df)
# print(txt)



result_sum_score_blast.write(txt)
result_sum_score_blast.close()	



# result_sum_score = open(os.path.join(input_dir,'result_sum_score.txt'),'r')
result_sum_score_blast = open(os.path.join(input_dir,'result_sum_score_blast.txt'),'r')
lines= result_sum_score_blast.readlines()

seq_set=[]
seq_cluster={}


for line in lines:
	if len(line) <1 :continue
	l=line.rsplit()
	if len(l)<2 : continue
	if '#' in l[0]: continue
	if l[3] not in seq_set:
		# print(l[3])
		seq_set.append(l[3])

for seq in seq_set:
	seq_cluster[seq]=[]
	for line in lines:
		if len(line)<1:continue
		l=line.rsplit()
		if len(l)<2:continue
		if l[3] == seq:
			seq_cluster[seq].append([l[0],l[2],l[4]])

# print(seq_cluster)
slim_seq_cluster=seq_cluster.copy()
for seq in slim_seq_cluster.items():
	for design in seq[1][:]:
		if float(design[2]) > 0.1:
			seq[1].remove(design) 

# for seq in slim_seq_cluster.items():
# 	print('--------------')
# 	if len(seq[1]) >1:
# 		output_file_dir=os.path.join(input_dir,'{}'.format(seq[0]))
# 		if not os.path.exists(output_file_dir):
# 			os.mkdir(output_file_dir)
		
# 		for design in seq[1]:				
# 			outfile=os.path.join(output_file_dir,'{}'.format(design[0]))
		
# 			if not os.path.exists(outfile):
# 				shutil.copyfile(os.path.join(input_dir,design[0]),outfile)
# 				print(design)


for seq in slim_seq_cluster.items():
	print('--------------')		
	for design in seq[1]:				
			print(design)





		 

