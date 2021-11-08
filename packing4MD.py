from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import *
import sys, os
import pack_analysisV3 as pack


#####################INPUTSECTION##############################
input_dir=sys.argv[2]
input_mdoutput=sys.argv[3]



# input_dir='/home/sting-gpu/weiyi/design/MD_input_files/design_051_0007/wkr/lastFrame_run_nvt.pdb'
# input_dir='/home/sting-gpu/weiyi/design/MD_input_files/design_051_0007/wkr/lastFrame_run_em.pdb'

# input_dir='/home/sting-gpu/weiyi/design/MD_input_files/design_051_0007/wkr/lastFrame_run.pdb'
# input_dir='/home/sting-gpu/weiyi/design/MD_input_files/design_059-0002/design_setup5_059-0002.pdb'

chainABC_str = 'DE'
chainX_str='A'
scorecutoff=0
thirdfilter=0

###########################################################3


# print(pdb_files)

input_PDB_name = input_dir.split('/')[-1]
print(input_PDB_name)
parser   = PDBParser(QUIET=1)
protein  = parser.get_structure(input_PDB_name,input_dir)[0]

residues=list(protein.get_residues())
# print(residues)
for res in residues:
    resid=res.get_id()
#     print(resid)
    resname=res.get_resname()
    if resname == "ACE" or resname == "NME":
#         print(resname)
        res.get_parent().detach_child(resid)
#         print(resid)



sum1,sum2,sum3,sum4,sum5,sum6,finalscore   =pack.excute(protein,chainABC_str,chainX_str,scorecutoff,thirdfilter)
print(sum1,sum2,sum3,sum4,sum5,sum6,finalscore )