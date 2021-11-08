import sys,os
# sys.path.append('/Users/tangweiyi/Desktop/lib_py')
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
import pack_analysisV4_4traj as pack
# import matplotlib.pyplot as plt
# import itertools
# import seaborn as sns
# import shutil
from collections import defaultdict

# class Traj_HKpair:

#     def __init__(self,pair,exist,time,contactparm):
#         self.pair=pair
#         self.name=str(pair)
#         # self.exist=np.zeros(8,dytpe=int)
#         self.exist=0
#         self.time=0
#         self.contactparm=0

# #####################INPUT SECTION for large output set#############################
# input_dir=sys.argv[1]
# chainABC_str = 'FB'
# chainX_str='A'
# scorecutoff=0
# thirdfilter=0



# dir_files = os.listdir(input_dir)
# traj_pdb_files = [x for x in dir_files if x[-3:] == 'pdb' and x[:5] =='frame']
# ori_pdb_file = [x for x in dir_files if x[-3:] == 'pdb' and d[:6] == 'design'][0]

# traj_pdb_files.sort()
# HK_set=[]

# # pdb_files = [x for x in dir_files if x[-3:] == 'pdb']


# # result_sum_score= open(os.path.join(output_dir,'result_sum_score.txt'),'w')
# # txt = '\n#  Good design finalscore is above ' +str(scorecutoff)+'\n'
# # txt+= '#  Design_name          finalscore   #packed_knobs\n'

# # result_w_detail=open(os.path.join(output_dir,'result_w_detail.txt'),'w')
# # txt2 = '\n#  Good design finalscore is above ' +str(scorecutoff)+'\n'
    
# # print('\n -----------------------\n the output directory is built.')

# ori_protein  = parser.get_structure(ori_pdb_file,os.path.join(input_dir,ori_pdb_file))[0]
# for contactparm2 in np.linspace(2.8,4,4):
#         sum5,sum6,finalscore,packingnames   = pack.excute(protein,"AB","X",scorecutoff,thirdfilter,contactparm2)  
#         for pair in packingnames:
#             pair=Traj_HKpair(pair,1,-100,contactparm2)
#             HK_set.append(pair)


# # result=defaultdict()


# for num,input_PDB_name in enumerate(traj_pdb_files): 
    
#     time=input_PDB_name[-7:-3]
#     protein  = parser.get_structure(input_PDB_name,os.path.join(input_dir,input_PDB_name))[0]
#     residues=list(protein.get_residues())

#     for res in residues:
#         resid=res.get_id()
#         resname=res.get_resname()
#         if resname == "ACE" or resname == "NME":
#             res.get_parent().detach_child(resid)

#     for contactparm2 in np.linspace(2.8,4,4):
#         sum5,sum6,finalscore,packingnames   = pack.excute(protein,chainABC_str,chainX_str,scorecutoff,thirdfilter,contactparm2)
#         for pair in packingnames:
#             pair=Traj_HKpair(pair,1,time,contactparm2)
#             HK_set.append(pair)
           
#             # print(pair[0],pair[1])
#             # txt2+='{}  {}\n'.format(pair[0],pair[1])

#     #   print(sum3,sum4,sum5,sum6,finalscore)
#         #Huge matrix for further plotting


        
#         # if finalscore >scorecutoff:
            
#         #     #print out the result 
#         #     txt2+= '###    ' +str(protein.get_parent().id)+'\n'
#         #     txt2+='*** The finalscore is %.1f,  there are %.1f knobs packing into holes\n' %(finalscore, sum6+sum5) 
#         #     print('*** The finalscore is ', finalscore, 'there are',sum6+sum5, 'knobs packing into holes')
#         #     for pair in packingnames:
#         #         print(pair[0],pair[1])
#         #         txt2+='{}  {}\n'.format(pair[0],pair[1])
            
#         #     #writing all the print-out to a result txt
#         #     txt+= str(input_PDB_name)+'         %.1f        '%finalscore + '%.1f'%(sum5+sum6)+' \n'
#     #         print(txt)
        

# for pair in HK_set:



# # result_sum_score.write(txt)
# # result_w_detail.write(txt2)
# # result_sum_score.close()
# # result_w_detail.close()


# class Traj_HKpair:

#     def __init__(self,pair,exist,time,contactparm):
#         self.pair=pair
#         self.name=str(pair)
#         # self.exist=np.zeros(8,dytpe=int)
#         self.exist=0
#         self.time=0
#         self.contactparm=0

#####################INPUT SECTION for large output set#############################
input_dir=sys.argv[1]
chainABC_str = 'FB'
chainX_str='A'
scorecutoff=0
thirdfilter=0


parser = PDBParser(QUIET=True)
dir_files = os.listdir(input_dir)
traj_pdb_files = [x for x in dir_files if x[-3:] == 'pdb' and x[:5] =='frame']
ori_pdb_file = [x for x in dir_files if x[-3:] == 'pdb' and x[:6] == 'design'][0]

traj_pdb_files.sort()
print(traj_pdb_files)

# pdb_files = [x for x in dir_files if x[-3:] == 'pdb']

# result_w_detail=open(os.path.join(output_dir,'result_w_detail.txt'),'a+')



traj_pack_w_detail=open(os.path.join(input_dir,'traj_pack_w_detail.txt'),'w')
txt = 'This is the detailed packing information'
txt+= '###    ' +ori_pdb_file+'\n'


hk_set=defaultdict(list)

for contactparm2 in np.linspace(2.8,4,4):
    

    txt += '\n\n\nThe contact distance cutoff is {}\n'.format(contactparm2)

    
    ori_protein  = parser.get_structure(ori_pdb_file,os.path.join(input_dir,ori_pdb_file))[0]
    sum5,sum6,finalscore,packingnames   = pack.excute(ori_protein,"AB","X",scorecutoff,thirdfilter,contactparm2)  
    for pair in packingnames:
        # txt2+='{}  {}'.format(pair[0],pair[1]) + "  -100  \n"
        pair=str(pair).replace("A","F")
        pair=pair.replace("X","A")
        hk_set[str(pair)].append(-100)

    for num,input_PDB_name in enumerate(traj_pdb_files): 
        
        time=input_PDB_name[-7:-3]
        protein  = parser.get_structure(input_PDB_name,os.path.join(input_dir,input_PDB_name))[0]
        residues=list(protein.get_residues())
        # print(input_PDB_name)
        # print(residues)
        for res in residues:

            resid=res.get_id()
            resname=res.get_resname()
            # print(resid)
            if resname == "ACE" or resname == "NME":
                res.get_parent().detach_child(resid)
                # print(resid)


        sum5,sum6,finalscore,packingnames   = pack.excute(protein,chainABC_str,chainX_str,scorecutoff,thirdfilter,contactparm2)
        for pair in packingnames:
            hk_set[str(pair)].append(time)

    for key, value in hk_set.items():
        txt+= key +"   "
        for v in value:
            txt += str(v).replace(".","")+"   "

        txt +='\n'
    hk_set=defaultdict(list)


traj_pack_w_detail.write(txt)
traj_pack_w_detail.close()
           

        




print(txt)



























