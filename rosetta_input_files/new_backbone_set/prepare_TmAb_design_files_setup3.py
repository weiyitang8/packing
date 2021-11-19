import sys, os, numpy as np, shutil, random, re
from PDButilv3 import *
from Ez_potential import *
from collections import Counter

# example usage:
# >> python ../../bin/prepare_TmAb_design_files.py ./ ../../bin/9aa_helix_ideal_ACE-CT3-termini.pdb ./ design_set1_input.txt 


## inputs
target_region_files 	= sys.argv[1]
ideal_helix_path		= sys.argv[2]
design_base_dir			= sys.argv[3]
build_file				= sys.argv[4]


TM_AAs		= [ 'A', 'G', 'I', 'L','F', 'S', 'T', 'V' ]
weightAA	= [0.16, 0.07, 0.14, 0.24, 0.11, 0.07, 0.07, 0.14]
#TM_AAs		= [ 'A', 'G', 'I', 'L','F', 'S', 'T', 'V' ]
#weightAA 	= [0.125, 0.075, 0.125, 0.35, 0.05, 0.075, 0.075, 0.125]



# I/O of target file inputs.  Must rename your files to this format; check if exists & quit if not
# target_region_*; e.g. target_region_resfile.txt, target_region.pdb, target_region_constraints.txt
ideal_helix 		= parsePDB(ideal_helix_path)
target_files = [ os.path.join(target_region_files , x) for x in ['fullTarget_resfile.txt', 'fullTarget_relaxed.pdb', 'fullTarget_constraints.txt', 'target_region.pdb'] ]
for i in target_files:
	if not os.path.exists(i):
		print('\nPROBLEM w/ inputs, required file does not exist in location: %s\n quitting... get it together user' % i )
		sys.exit()
target_resfile, full_target_PDBpath, target_constraint, target_regionPDBpath = tuple(target_files)

target_regionPDB 	= parsePDB( target_regionPDBpath ) 
full_targetPDB 		= parsePDB( full_target_PDBpath )


#  extract residue ranges & chains used from the fragment of the target used for PDB binder matching. 
#  assuming continuous residue ranges...
target_chains 	=  sorted( list(set(target_regionPDB.getChids())))
byChain_selString, sel_str, step	= {}, '', 0
for ch in target_chains:
	resi_set = target_regionPDB.select( 'ca chain %s ' % ch ).getResnums()
	byChain_selString[ch] = 'chain %s resnum %s' % ( ch, ' '.join( [ str(x) for x in resi_set] ) )
	sel_str += '( chain %s resnum %s ) ' % ( ch, ' '.join( [ str(x) for x in resi_set] ) ) 
	step +=1
	if step >= len (target_chains): 
		break
	else:
		sel_str += ' or '

# calculate the transformation to match the target region to the full target, relaxed/oriented in rosetta
design_target_frag_selection = target_regionPDB.select( sel_str )
part_of_full_target_to_align = full_targetPDB.select( sel_str )
empty, matrix2target = superpose( design_target_frag_selection, part_of_full_target_to_align )

# heavy atoms within target for distance matrix with interacting docked TM 
target_forDistances = target_regionPDB.select('name N C CA O CB CG')


# read input resfile contents - template to add design to
target_resfile_txt = ''
for x in open(target_resfile):
	target_resfile_txt +=x


### read input file with pdb candidate path, noting if helix residues (poly-ala) added to N or C term 
for i in open(build_file):
	# skip header or blank end of file
	if len(i) < 10 or i[0] == '#':	
		continue		
	# write inital pdb, helix-extended PDB, resfile, & light coordinate constriant file for first design round
	line = tuple( i.split() )
	frag_path, dir_prefix, Naa_2_extend, Caa_2_extend = line	
	Naa_2_extend, Caa_2_extend = int(Naa_2_extend), int(Caa_2_extend)
	binder_pdb_path 	= frag_path
	binder_workingDir  	= os.path.join( design_base_dir	 , dir_prefix + '_' + os.path.basename(frag_path)[:-4] ) 

	extendPDB_path 		= os.path.join( binder_workingDir ,  'ext_n%d_c%d-%s' % ( Naa_2_extend, Caa_2_extend,  os.path.basename(frag_path) ) )
	merged_target_TMAb  = os.path.join( binder_workingDir ,  'design_input_struct.pdb')
	res_path 	= os.path.join( binder_workingDir ,  'resfile.txt')
	resHB_path	= os.path.join( binder_workingDir ,  'resfile_HB.txt')

	if not os.path.exists( binder_workingDir ): 
		os.mkdir( binder_workingDir )

	## do not overwrite directories if they have already been prepared...
	if os.path.exists( res_path ):
		print ('\t already wrote: ', res_path, '... quiting')
		continue

	print (	binder_workingDir)
	shutil.copy( binder_pdb_path, binder_workingDir )

	# add short helical extension to existing PDB, as requested in build file
	if Naa_2_extend > 6 or Caa_2_extend > 6: 
		print('\nError, termini extension too long. Use only 6 or less\n quitting... '  )
		sys.exit()
	pdb 					 	= parsePDB( frag_path )
	step = 110
	pdb.setChids(['X' for x in pdb.getChids()])
	pdb.setSegnames(['X' for x in pdb.getSegnames()])
	for r in pdb.iterResidues():
		r.setResnum(step)
		step+=1
	last_3_Nalign 			 	= pdb.select('ca')[:3]
	last_3_Calign 			 	= pdb.select('ca')[-3:]
	start_resnum, end_resnum 	= last_3_Nalign[0].getResnum(), last_3_Calign[-1].getResnum()

	helix_ext 					= ideal_helix.copy()
	ext_last_3_N 				= ideal_helix.select('ca')[ Naa_2_extend  : 3 + Naa_2_extend ]
	ext_last_3_C 				= ideal_helix.select('ca')[- 3 - Caa_2_extend :  - Caa_2_extend ]

	nTransformation 			= calcTransformation( ext_last_3_N , last_3_Nalign )
	ext_last_3_N_ALL			= applyTransformation( nTransformation, ideal_helix.select( 'resnum 0 to %d' % (Naa_2_extend)  ).copy() ).copy()
	N_start_resnum 				= start_resnum - Naa_2_extend - 2
	for r in ext_last_3_N_ALL.iterResidues():
		r.setResnum(N_start_resnum + 1)
		N_start_resnum += 1	

	cTransformation 			= calcTransformation( ext_last_3_C , last_3_Calign )
	ext_last_3_C_ALL			= applyTransformation( cTransformation, ideal_helix.select( 'resnum %s to 10' % (10-Caa_2_extend) ).copy()   ).copy()
	C_start_resnum 				= end_resnum +1
	for r in ext_last_3_C_ALL.iterResidues():
		r.setResnum(C_start_resnum)
		C_start_resnum += 1	

	## transform this TM fragment to the relative conformation / geometry relative to the full target PDB 
	tmAb_template 	= ext_last_3_N_ALL + pdb + ext_last_3_C_ALL
	applyTransformation( matrix2target , tmAb_template )
	tmAb_merged		= full_targetPDB.select('protein').copy()  + tmAb_template 
	writePDB( merged_target_TMAb, tmAb_merged )

	#  auto-detect which residues to design for interaction interface
	#   if residue is gly, compute distance to closest backbone or CB/CG atom on target
	# (this uses target's modeled sidechain rotamers, which may or may not change during protein design trials) 
	design_flg_list = []
	seq 			= [ '-' for x in np.arange(len(tmAb_template.select('ca')) + 2) ]
	last_six_indx 	= [ len(tmAb_template.select('ca')) - x for x in [0,1,2,3,4,5]  ]
	step 			= 1
	n4aa_EZ, c4aa_EZ= {}, {}

	for r in tmAb_template.iterResidues():
		print(step)
		design_flg 	= 0
		CA_CB 		= r.select( 'name CA CB ')

		#  ignore N & C caps
		if not CA_CB:
			design_flg_list.append( design_flg )
			
		#  handle amino acids, assign as interfacial (designable) or lipid-facing (randomized)
		#  for first and last 4 amino acids, use Z value to 'best' position potential ARG, LYS, TRP, TYR
		else:
			notGly_flag	= len(CA_CB) - 1
			dMat = buildDistMatrix( CA_CB , target_forDistances )
			closest_SC_dist = np.min(dMat)

			if notGly_flag:
				CB_vs_CA_dmat 	= np.mean( dMat[1] - dMat[0] )		# is the CB of sidechain closer to interface than CA?
				design_flg 		= 0

				if closest_SC_dist < 8:
					design_flg = 1
				elif closest_SC_dist < 10 and CB_vs_CA_dmat < -0.2:
					design_flg = 1
				elif closest_SC_dist < 8.5 and CB_vs_CA_dmat < 0:
					design_flg = 1
				else:
					pass

			else:			# glycines
				if closest_SC_dist < 8:
					design_flg = 1

			design_flg_list.append( design_flg )
			seq[step] = '_'

			## EZ-potential
			if not design_flg:

				Z 		= np.fabs(  calcCenter(CA_CB)[2]  )
				print(Z)
				probs 	= [ round( compute_Ez(x, Z), 4 ) for x in ['TRP', 'TYR', 'ARG', 'LYS', 'PHE'] ]
				probs.append(Z)

				if Z > 13.5 :					# if outside hydrophobic core or headgroup region, default to Lys 
					seq[step] = 'K'
				#	print(r, step, '-', round(Z, 2) )

				elif Z < 10 :						# if well inside hydrophobic core, default to semi-randomized AA
				#	print(r, step, '-', round(Z, 2) )
					pass

				elif step in [0,1,2,3,4,5,6]:		# if at headgroup region, calculate EZ potential for F, Y, W
					n4aa_EZ[step] =  probs
				#	print( r, step, '-',round(Z, 1), probs)

				elif step in last_six_indx:			# if at headgroup region, calculate EZ potential for F, Y, W
					n4aa_EZ[step] =  probs
				#	print( r, step, '-',round(Z, 1), probs)
				else:
				#	print(r, step, '-', round(Z, 2) )
					pass
			else:
			#	print(r, step, 'X')
				seq[step] = 'X'

			
			step += 1

	#print( [ round( compute_Ez(x, 20), 4 ) for x in ['TRP', 'TYR', 'ARG', 'LYS'] ] )
	step = 0
	for n, r in enumerate( tmAb_template.iterResidues() ):
	#	print (n, r, design_flg_list[n])
		pass

	print(n4aa_EZ)
	W_top = sorted( n4aa_EZ.items(), key=lambda x: x[1][0] )[0][0]
	print (W_top)
	seq[W_top] = 'W'
	del n4aa_EZ[W_top]

	# Y_top = sorted( n4aa_EZ.items(), key=lambda x: x[1][1] )[1][0]
	# print (Y_top)
	# seq[Y_top] = 'Y'
	# del n4aa_EZ[Y_top]

	

	numRandom_AA = len( [ x for x in seq if x == '_'] )


	## randomization, quit only upon returning a success flag
	# no G/S/AxxxG/S/A motifs or 'G/A/SxxxxxxG/S/A', G/S/TxxxT/G/S,   
	# LxxI(V/T/C/I)xxL is also rejected
	known_motif_tegex = []
	randomTMaa_success_flg = 0
	#while randomTMaa_success_flg == 0:
	while randomTMaa_success_flg == 0:
		step = 0
		random_AAs = random.choices(TM_AAs	,  weights=tuple(weightAA), k=numRandom_AA)

		seq2 = []
		for i in seq :
			if i=='_':
				seq2.append( random_AAs[step] )
				step +=1
			else:
				seq2.append( i )

		seq2_str = ''.join( seq2 ) 

		## sequence counter filters
		## no more than 1 of S,T,G; no more than 4 L
		counts = Counter( seq2 )
		if counts['S'] > 2 or counts['G'] > 1 or counts['T'] > 2 or counts['L'] > 4 or counts['I'] > 3 or counts['F'] < 2:
			#print (counts.most_common(), '--------------BAD SEQ, count filter broken!!! \n')
			continue

		## sequence pattern filters
		match = re.search(r"[GAS]\w\w\w[GAS]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 1 filter broken!!! \n')
			continue

		match = re.search(r"[GST]\w\w\w[GST]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 2 filter broken!!! \n')
			continue

		match = re.search(r"[GSA]\w\w\w\w\w\w[GSA]", seq2_str )
		if match:
			#print ( '--------------BAD SEQ, seq pattern 3 filter broken!!! \n')
			continue
		
		randomTMaa_success_flg = 1



	print ('\t', seq2_str)
	# write resfile & constraint file


	# copy constraints for target.  write light constraint file for core of TM domain matched
	shutil.copy( target_constraint, binder_workingDir)

	### TBD section .... write constraints for the for core of TM domain matched
	#for r in pdb.copy().iterResidues ():
	#cstStr += 'CoordinateConstraint CA %s N %s %f %f %f HARMONIC 0.0 0.3\n' % ( r[0], ca[0][0], r[1][0], r[1][1], r[1][2] )



	# resfile
	header_design = '### TMab - DE NOVO design, lipid-facing residues fixed from semi-random computer selection'
	txt, txt_HB = '\n%s\n' % header_design, '\n%s\n' % header_design
	cstStr = ''
	residue_CA 		= tmAb_template.select('ca')
	seq2_str 		= seq2_str[1:-1]
	start_resnum 	= residue_CA[0].getResnum()
	for n in np.arange( len(residue_CA) ):
		#print( n , residue_CA[n], residue_CA[n].getResname(), seq2_str[n] )
		assignment = seq2_str[n]
		resnum = residue_CA[n].getResnum()
		if  assignment == 'X':
			txt 	+= '%d %s PIKAA AFGILMYWSTCV \n' % ( resnum, 'X' )
			txt_HB 	+= '%d %s PIKAA AFGILMYWSTCVNQH \n' % ( resnum, 'X' )
		else: 
			txt 	+= '%d %s PIKAA %s \n' % ( resnum, 'X', assignment )
			txt_HB 	+= '%d %s PIKAA %s \n' % ( resnum, 'X', assignment )

	resfile 	= target_resfile_txt + txt
	resfile_HB 	= target_resfile_txt + txt_HB

	for p, t in zip( [res_path, resHB_path], [resfile, resfile_HB] ):
		f = open(p, 'w')
		f.write(t)
		f.close()

sys.exit()


# write resfile 
txt = 'start\n'
for k,v in chResi.items():
			if k in ['X', 'Y', 'Z']:
				for r in v:
					txt += '%d %s ALLAAxc\n' % ( r, k )
			else:
				for r, c in zip( v, IntF ):
					if c: txt += '%d %s NATAA\n' % ( r, k )
					else: txt += '%d %s PIKAA G\n' % ( r, k )

# write constraint file 
if ca[0] == r: # skip the first residue, since it is the reference coordinate for constraints
			cstStr += 'CoordinateConstraint CA %d CA %d %f %f %f HARMONIC 0.0 0.3\n' % ( r[0], ca[0][0], r[1][0], r[1][1], r[1][2] )
		
