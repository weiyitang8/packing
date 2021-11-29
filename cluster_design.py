from prody import *
from pylab import *
import sys, os
input_base='/home/sting-gpu/weiyi/design/Rosetta/setup11/goodones/'
# test1= parsePDB('/home/weiyi/design/Rosetta/setup11/goodones/design_input_struct_003_0002')


test1= parsePDB(os.path.join(input_base,'design_input_struct_003_0002.pdb'))
test2= parsePDB(os.path.join(input_base,'design_input_struct_003_0003.pdb'))

showProtein(test1,test2)
matches=matchChains(test1,test2)
for match in matches:
	print(match[0],"\n",match[1])
	print(calcRMSD(match[0],match[1]))
