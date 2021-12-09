import os, sys,numpy as np, pandas as pd,shutil
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pylab as plt
from collections import defaultdict

input_dir= sys.argv[1]
X_1st_resnum=int(sys.argv[2])
dir_files=os.listdir(input_dir)
pdb_files= [x for x in dir_files if x[-3:]=='pdb']

AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


L=np.array([0,0,0,0,0])
I=np.array([0,0,0,0,0])
M=np.array([0,0,0,0,1])
V=np.array([0,0,0,0,2])

C=np.array([0,0,0,1,0])

A=np.array([0,2,0,0,0])
G=np.array([0,2,0,1,0])

S=np.array([0,2,1,0,0])
T=np.array([0,2,1,0,1])

P=np.array([0,2,2,0,1])

F=np.array([0,4,0,0,0])
Y=np.array([0,4,0,0,1])
W=np.array([0,4,0,0,2])

E=np.array([-2,-2,-1,0,0])
D=np.array([-2,-2,-1,0,0])
N=np.array([-2,-2,-1,0,-1])
Q=np.array([-2,-2,-1,0,-1])

K=np.array([-2,-2,-1,0,0])
R=np.array([-2,-2,-1,0,-1])
H=np.array([-2,-2,-2,0,0])

AA_index={'L':L,'I':I,'V':V,'M':M,'C':C,'A':A,'G':G,'S':S,'T':T,'P':P,'F':F,'Y':Y,'W':W,'E':E,'D':D,'N':N,'Q':Q,'K':K,'R':R,'H':H}


parser = PDBParser(QUIET=True)


seq_set=list(np.arange(len(pdb_files)))
name_set=list(np.arange(len(pdb_files)))
sequence_set=defaultdict(list)


for i, pdb in enumerate(pdb_files):
	
	pdb_name=pdb[-12:-4]
	name_set[i]=pdb_name
	path=os.path.join(input_dir,pdb)
	
	structure=parser.get_structure('struct',path)
	structure=structure[0]
	chainX=structure['X']
	chainX_resnum=len(list(chainX.get_residues()))
	seq=list(np.arange(chainX_resnum))
	sequence=''
	
	step=0
	for res in chainX:
		# print(type(AA_index[AA[res.get_resname()]]))

		# print(type(AA[res.get_resname()]))
		# print(seq[step])
		sequence+=(AA[res.get_resname()])
		seq[step]=AA_index[AA[res.get_resname()]]
		step+=1
	seq_set[i]=seq
	sequence_set[pdb_name].append(sequence)
print(len(seq_set))

seq_set=np.array(seq_set)

print(seq_set.shape)

# print(seq_df)
# disMat=ssd.pdist(seq_set,metric='euclidean')
l,m,n=seq_set.shape
print(l,m,n)
D=np.zeros([l,l])
for i in range(l):
	for j in range (i+1,l):
		D[i,j]=np.linalg.norm(seq_set[i]-seq_set[j])
		D[j,i]=D[i,j]


# print(D)
# print(len(D))
# print(D.mean())
D=ssd.squareform(D)
Z=sch.linkage(D,method='complete',metric='distance')
P=sch.dendrogram(Z,labels=name_set)
plt.show()


h_clust=sch.fcluster(Z,t=2.5,criterion='distance')
# print(len(set(h_clust)))
print(h_clust)


design_clusters=defaultdict(list)
# design_cluster_key=defaultdict(list)

for n in np.arange(len(name_set)):
	design_clusters[h_clust[n]].append(name_set[n])
	# design_cluster_key[h_clust[n]].append(design_names_sorted[n])


for k,v in sorted(design_clusters.items()):
	print('this is cluster Number', k)
	for p in v:

		print(p,sequence_set[p])



singletons_dir=os.path.join(input_dir,'singletons')
if not os.path.exists(singletons_dir):
	os.mkdir(singletons_dir)

for cluster_index,v in sorted(design_clusters.items()):
	# print(cluster_index,v)
	cluster_size=len(v)
	if cluster_size>1:
		clust_dir=os.path.join(input_dir,'seq_cluster_%s_%d')%(cluster_index,cluster_size)
		if not os.path.exists(clust_dir):
			os.mkdir(clust_dir)
		for design in v:
			design_path = os.path.join(input_dir,'design_input_struct_%s.pdb')%design
			shutil.copy(design_path,clust_dir)
	else:
		design_path = os.path.join(input_dir,'design_input_struct_%s.pdb')%v[0]
		shutil.copy(design_path,singletons_dir)




















# input_MD_dir= sys.argv[1]




# structure=parsePDB(os.path.join(input_MD_dir,'input.pdb'))
# repr(structure)
# ensemble=parseDCD


