import sys,os
# sys.path.append('/Users/tangweiyi/Desktop/lib_py')
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
import pack_analysisV3 as pack
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
import shutil


###Run multiple Rosetta Design output and make the plot
### input 1) Path to directory of the rosetta output PDB files
### input 2) the cutoff of finalscore to decide what hole-knob-pair details you what to print out 
### input 3) the cutoff of filter to deicde whether this design worth going through 
### input 4) output directory name for the GOOD/pass the thirld filter proteins and result txt files.
### input chains(no spacing) of the target protein interface eg. AB by default
### input chain(single letter) of the "design TM antibody" eg. X by default


###The plotting part was meant to use on the jupyter notebook though

#################INPUT SECTION for Jupyter notebook#############################3
# input_dir=sys.argv[3]
# chainABC_str = sys.argv[4]
# chainX_str=sys.argv[5]
# scorecutoff=float(sys.argv[6])



#####################INPUT SECTION for large output set#############################
input_dir=sys.argv[1]
chainABC_str = 'AB'
chainX_str='X'
scorecutoff=int(sys.argv[2])
thirdfilter=int(sys.argv[3])
output_dir=sys.argv[4]
######################INPUT SECTION for single run##############################
# input_dir='/home/sting-gpu/weiyi/design/Rosetta/setup5_output2/good_design'  
# chainABC_str = 'AB'
# chainX_str='X'
# scorecutoff=7
# thirdfilter=7
# output_dir='/home/sting-gpu/weiyi/design/Rosetta/setup5_output2/good_packing/dalphaball'
###########################################################

dir_files = os.listdir(input_dir)
pdb_files = [x for x in dir_files if x[-3:] == 'pdb']

# print(pdb_files)

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

result_sum_score= open(os.path.join(output_dir,'result_sum_score.txt'),'w')
txt = '\n#  Good design finalscore is above ' +str(scorecutoff)+'\n'
txt+= '#  Design_name          finalscore   #packed_knobs\n'

result_w_detail=open(os.path.join(output_dir,'result_w_detail.txt'),'w')
txt2 = '\n#  Good design finalscore is above ' +str(scorecutoff)+'\n'
    
print('\n -----------------------\n the output directory is built.')
    
parser   =PDBParser(QUIET=1)
num_knobs=np.zeros((len(pdb_files),7))
total_num=len(pdb_files)

        
    
for num,input_PDB_name in enumerate(pdb_files): 
    if num%10==0:
        print('now it is ',num, 'of',total_num) #just for impatient weiyi
#   print(input_PDB_name)
    protein  = parser.get_structure(input_PDB_name,os.path.join(input_dir,input_PDB_name))[0]
    sum1,sum2,sum3,sum4,sum5,sum6,finalscore,packingnames   =pack.excute(protein,chainABC_str,chainX_str,scorecutoff,thirdfilter)
# #   print(sum3,sum4,sum5,sum6,finalscore)
#     #Huge matrix for further plotting
#     num_knobs[num][0]=sum1  #holeABC pass the first filter, vector method
#     num_knobs[num][1]=sum2  #holeX pass the first filter, vector method
#     num_knobs[num][2]=sum3  #holeABC pass the second filter, sasa/chain
#     num_knobs[num][3]=sum4  #holeX pass the second filter, sasa/chain
#     num_knobs[num][4]=sum5  #holeABC pass the second filter,sasa/knob
#     num_knobs[num][5]=sum6  #holeX pass the second filter, sasa/knob
#     num_knobs[num][6]=finalscore #final score, fully packing +1, partially packing +0.5
    
    if finalscore >scorecutoff:
        
        #print out the result 
        txt2+= '###    ' +str(protein.get_parent().id)+'\n'
        txt2+='*** The finalscore is %.1f ,there are %.1f knobs packing into holes\n' %(finalscore, sum6+sum5) 
        print('*** good design,finalscore is above', scorecutoff,  '   ###', protein.get_parent().id)
        print('*** The finalscore is ', finalscore, 'there are',sum6+sum5, 'knobs packing into holes')
        for pair in packingnames:
            print(pair[0],pair[1])
            txt2+='{}  {}\n'.format(pair[0],pair[1])
        
        #writing all the print-out to a result txt
        txt+= str(input_PDB_name)+'         %.1f        '%finalscore + '%.1f'%(sum5+sum6)+' \n'
#         print(txt)
        
        
        #copy the good packing designs to a new directory
        outfile=os.path.join(output_dir,input_PDB_name)
        if not os.path.exists(outfile):
            print(outfile, 'is copying to the good_packing folding')
            shutil.copyfile(os.path.join(input_dir,input_PDB_name),outfile)


result_sum_score.write(txt)
result_w_detail.write(txt2)
result_sum_score.close()
result_w_detail.close()
