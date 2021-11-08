
import os
topol_names	 = [ x for x in os.listdir( './' ) if 'topol_' in x and '.itp' in x ]
posre_names  = [ x for x in os.listdir( './' ) if 'posre_' in x and '.itp' in x ]

# text_4_files, write_or_not_flag=[],[]
for topol in topol_names:
	# text_already_there_flag=0

	for posre in posre_names:
		chain_name=posre[-11:-4]
		if chain_name == topol[-11:-4]:
			posre_text='\n\n; Include positions restraints for %s \n' %chain_name
			posre_text+='\n#ifdef POSRES \n'
			posre_text+='#include "{}"'.format(posre)
			posre_text+='\n#endif\n'

			oF = open(topol,'a')
			# if not '#ifdef POSRES' in oF: continue
			oF.write(posre_text)
			oF.close()

			print('inclue POSRES for %s' %chain_name)
			







