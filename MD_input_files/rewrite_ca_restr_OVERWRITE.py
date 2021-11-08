import os, sys, numpy as np


## generate position restraints for CA atoms in each chain
#  and append as "#include" note to topology .itp for each chain
 
# collect these topol_chain_x.itp files 
filenames	 = [ x for x in os.listdir( './' ) if 'topol_' in x and '.itp' in x ]

text_4_files, write_or_not_flag = [], []
for f in filenames:

	text_already_there_flag = 0
	# text to append to .itp file for turning on CA restraints
	posre_CA_text = '\n\n; Include position restraints for C-alpha of Protein\n'
	posre_CA_text += '#ifdef POSRES_CA\n\n[ position_restraints ]\n'
	posre_CA_text += ';  i funct       fcx        fcy        fcz\n' 

# open  file to find the atom indecies of CA atoms
	oF = open( f, 'r')
	for line in oF.readlines():
		if '#ifdef POSRES_CA' in line: text_already_there_flag = 1
		l = line.split()

		if len(l) == 8:

			
			if l[4] == 'CA':
				print (l)
				posre_CA_text += "{:>4}   1       2000       2000       2000\n".format( l[0] )


	posre_CA_text += '\n#endif\n'

	text_4_files.append(posre_CA_text)
	write_or_not_flag.append(text_already_there_flag)
	oF.close()


### edit the .itp files, appending Ca restraint text
for i,  txt , f in zip( filenames, text_4_files, write_or_not_flag):
	# print(i)

	if  '#' in i: continue		# skip back-up files

	if f: continue				# posre_CA_text already present in this .ipt file, don't double write

	oF = open( i, 'a')
	oF.write(txt)
	oF.close()
	print ('added 2kcal restr of CA to %s' % i )

