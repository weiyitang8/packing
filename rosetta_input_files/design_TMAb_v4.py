# Marco Mravic, Scripps 2021

## Input 1: target structure	path of this will be used to find span file w/ same name
## Input 2: path to rosetta main
## input 3: path to Rosetta scripts XML with protocol 
## input 4: path to resfile
## input 5: path to constraints

## Example command line
# python design_TMAb_v4.py design_input_struct-Full-target.pdb /home/sting-gpu/rosetta/main design_TMAb_protocol6.xml  resfile_setup5.txt fullTarget_constraints.txt setup5_output2/

# 
# python ../../../bin/design_TMAb_v1.py design_input_struct.pdb ~/rosetta/ ../../../bin/design_TMAb_protocol1.xml ./resfile.txt ./target_region_constraints.txt setup1_output
#
# python ../../../bin/design_TMAb_v3.py design_input_struct-Full-target.pdb ~/rosetta/ ../../../bin/design_TMAb_protocol3.xml ./resfile_setup3.txt ./fullTarget_constraints.txt setup3_output
#		162 seconds per run
# seq 1 3 | xargs -n 1 -P 3 python ../../../bin/design_TMAb_v3.py design_input_struct-Full-target.pdb ~/rosetta/ ../../../bin/design_TMAb_protocol5.xml ./resfile_setup5.txt ./fullTarget_constraints.txt setup5_output


import sys, os, subprocess as sp, numpy as np
from prody import *

inPDB		= sys.argv[1]
rosiBase 	= sys.argv[2]
protocolPth = sys.argv[3]
resfile_path= sys.argv[4]
cstfile_path= sys.argv[5]
outputDir 	= sys.argv[6]
if len(sys.argv) == 8:
	offset  = sys.argv[7]
else:
	offset  = 0

rosiScrps	= os.path.join( rosiBase, 'source/bin/rosetta_scripts.linuxgccrelease' )
rosiDB 		= os.path.join( rosiBase, 'database/')
rosiScoreF 	= os.path.join( rosiBase, 'database/scoring/score_functions/aa_composition/')
rosiHoles 	= os.path.join( rosiBase, 'source/external/DAlpahBall/DAlphaBall.gcc')

spanF 		= inPDB[:-4] + '.span'
if not os.path.exists( spanF ):
	print('ERROR no span file at: %s \n\n' % spanF )
	sys.exit()


#check the last file's suffix numbering.  Next files should be increasing index
if len( os.listdir(outputDir) ) == 0:
	pdbs= []
else:
	pdbs = [ x for x in os.listdir(outputDir) if 'pdb' in x ]

if len(pdbs) < 1:
	suffix = "_%s" % str(1).zfill(3)
else:
	suffices = [0]
	for x in pdbs:
		print (x, x.split('_') )
		suf = int(x.split('_')[-2])
		suffices.append(suf)
	suffix = "_%s" % str( np.max(suffices) + 1 + int(offset) ).zfill(3)



cmd = [  rosiScrps, 
'-parser:protocol', protocolPth, 			# Path to Rosetta script (see above)
'-in:file:s', inPDB,							# Input PDB structure   
'-nstruct', '20', 							# Generate 1 model
'-mp:setup:spanfiles', spanF,				# Input spanfile
'-relax:jump_move', 'true', 				# Allow jumps to move during relax
'-out:overwrite',
'-packing:resfile', resfile_path , 
'-packing:pack_missing_sidechains', '1',
'-ignore_zero_occupancy', 'false',
'-mp:lipids:composition', 'POPC',
'-out:path:all', outputDir,
'-parser:script_vars', 'cst_file=%s' % cstfile_path, 'AAcomp_file=./TMab_AA_composition.comp'
'-aa_composition_setup_file', os.path.join( rosiScoreF, 'TMab_AA_comp_less2-MWY.comp' ),
'-holes:dalphaball', rosiHoles, 
'-out:suffix', suffix
 ]

sp.call( cmd )

 