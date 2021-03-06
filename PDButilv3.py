## PDB file utilities ##
##################### Written By: Marco Mravic  ###  Degrado Lab UCSF Biophysics Aug 2015

import numpy as np, subprocess as sp, re
from prody import * 
from collections import Counter
from random import randint


###################### Globally defined helpers  ####################

divalent_metals = ['BA','BE', 'CD', 'CA', 'CR', 'CO', 'CU', 
'EU', 'GD', 'GE', 'FE', 'LA', 'PD', 'MG', 'MN', 'HG', 'NI', 
'OS', 'PT', 'RU', 'SR', 'SN','U', 'V', 'Y', 'ZN' ]

natAA = {}
natAA["ALA"] = 'A'; natAA["CYS"] = 'C'; natAA["ASP"] = 'D'; natAA["GLU"] = 'E'; natAA["PHE"] = 'F';
natAA["GLY"] = 'G'; natAA["HIS"] = 'H'; natAA["ILE"] = 'I'; natAA["LYS"] = 'K';
natAA["LEU"] = 'L'; natAA["MET"] = 'M'; natAA["ASN"] = 'N'; natAA["PRO"] = 'P'; natAA["GLN"] = 'N';
natAA["ARG"] = 'R'; natAA["SER"] = 'S'; natAA["THR"] = 'T'; natAA["VAL"] = 'V'; natAA["TRP"] = 'W'; natAA["TYR"] = 'Y';

# Arbitrary amino acid distribution from a database of PDBs
aaProp = {}
aaProp["ALA"] = 7.73; aaProp["CYS"] = 1.84; aaProp["ASP"] = 5.82; aaProp["GLU"] = 6.61; aaProp["PHE"] = 4.05;
aaProp["GLY"] = 7.11; aaProp["HIS"] = 2.35; aaProp["HSD"] = 2.35; aaProp["ILE"] = 5.66; aaProp["LYS"] = 6.27;
aaProp["LEU"] = 8.83; aaProp["MET"] = 2.08; aaProp["ASN"] = 4.50; aaProp["PRO"] = 4.52; aaProp["GLN"] = 3.94;
aaProp["ARG"] = 5.03; aaProp["SER"] = 6.13; aaProp["THR"] = 5.53; aaProp["VAL"] = 6.91; aaProp["TRP"] = 1.51; aaProp["TYR"] = 3.54;


TM_aaFrequency = {}
TM_aaFrequency['A'] = 0.1186; TM_aaFrequency['R'] = 0.0216; 
TM_aaFrequency['N'] = 0.0189; TM_aaFrequency['D'] = 0.0115;
TM_aaFrequency['C'] = 0.0101; TM_aaFrequency['Q'] = 0.0177;
TM_aaFrequency['E'] = 0.0162; TM_aaFrequency['G'] = 0.0819;
TM_aaFrequency['H'] = 0.0121; TM_aaFrequency['I'] = 0.1057;
TM_aaFrequency['L'] = 0.1514; TM_aaFrequency['K'] = 0.0173;
TM_aaFrequency['M'] = 0.0351; TM_aaFrequency['P'] = 0.0258; 
TM_aaFrequency['F'] = 0.0814; TM_aaFrequency['S'] = 0.0528;
TM_aaFrequency['T'] = 0.0544; TM_aaFrequency['W'] = 0.0283; 
TM_aaFrequency['Y'] = 0.0360; TM_aaFrequency['V'] = 0.1031;

## Unnatural/PTM Amino acid converter
#unNatAA = { 'ABA':'A', 
#'CSO':'C' , 'CSD':'C', 'CME':'C', 'OCS':'C', 
#"HSD":'H',
#'KCX':'K', 'LLP':'K', 'MLY':'K', 'M3L':'K', 
#'MSE':'M', 
#'PCA':'P', 'HYP':'P',
#'SEP':'S', 'TPO':'T', 'PTR':'Y'
#  }

# Same Three to One res ID converter but includes/converts unnatural/modified amino acids
UnNatAA={}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'N';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'


NatAA_3 = {}
NatAA_3['A'] = 'ALA'; NatAA_3['C'] = 'CYS'; NatAA_3['D'] = 'ASP'; 
NatAA_3['E'] = 'GLU'; NatAA_3['F'] = 'PHE'; NatAA_3['G'] = 'GLY'; NatAA_3['H'] = 'HIS'; NatAA_3['I'] = 'ILE'; NatAA_3['K'] = 'LYS'; 
NatAA_3['L'] = 'LEU'; NatAA_3['M'] = 'MET'; NatAA_3['N'] = 'ASN'; NatAA_3['P'] = 'PRO'; NatAA_3['Q'] = 'GLN'; 
NatAA_3['R'] = 'ARG'; NatAA_3['S'] = 'SER'; NatAA_3['T'] = 'THR'; NatAA_3['V'] = 'VAL'; NatAA_3['W'] = 'TRP'; NatAA_3['Y'] = 'TYR'; 


rename_Rosetta_termini = {'AcetylatedNtermProteinFull':'ACE-', 'CtermProteinFull':'-CO2', 'MethylatedCtermProteinFull':'-NME' }



#################### End of Global Definitions ####################



############### Class def's ##############


# mapped PDB chain containing which pdb residue corresponds to each uniprot entry 
#Has optiono to initialize or update a chain with a DBREF line (pdbfile) or list of segments from EMBL SIFTS mapping .tsv file 
class PdbChain():
	
	# does mapping at initialization given a DBREF line from a PDB file
	def __init__(self, init_flg, init_info):
		self.name 		= ''
		self.uniprot 	= ''
		self.pdbRes 	= np.array( [] )
		self.uniRes  	= np.array( [] )
		self.exceptFlg 	= False				# This flag is True if the   

		if init_flg == 'pdb':
			self.pdbLine_init(init_info)	# initialize with DBREF lines
		elif init_flg == 'SIFTS':
			self.SIFTS_init(init_info)	# initlize with list of segments in SIFTS .tsv file format
		else:						#Stay empty
			pass


	def __repr__(self):
		return '%s, %s' % (self.name, self.uniprot) 

	# For each chain in pdb entry, should have uniprot mapping at chain level (except large structures or synthesized constructs)
	def pdbLine_init(self, pdbLine):
		self.name = pdbLine[7:11] + pdbLine[12].lower() 
		self.uniprot = pdbLine[33:39]
		self.pdbRes = np.arange( int( pdbLine[14:18].strip() ) , int ( pdbLine[19:24].strip() ) + 1 )
		self.uniRes  = np.arange( int( pdbLine[55:60].strip() ) , int ( pdbLine[62:68].strip() ) + 1 )
		return

	# reads in format from EMBL-EBI mapping https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html, pdb to uniprot
	def SIFTS_init(self, lineList):
		if len(lineList) == 0: print ('WHY', self.name, self.uniprot)
		for i in lineList:
			ln 			= i.rstrip().split(',')
			chain 		= ''.join( [ln[0].upper(), ln[1].lower() ] )
			uniprot 	= ln[2]
			try:
				p_res 	= np.arange(  int( ln[5] ) , int( ln[6] ) +1 )
			except ValueError:
				p_res 	= np.array( [] )
				self.exceptFlg = True
			try:	
				u_res 		= np.arange(  int( ln[7] ) , int( ln[8] ) +1 )
			except ValueError:
				u_res 	= np.array( [] )
				self.exceptFlg = True

			self.pdbRes = np.append( self.pdbRes, p_res )
			self.uniRes = np.append( self.uniRes, u_res )

		self.name		= chain
		self.uniprot 	= uniprot
		return


	def pdb2uni(self, resNum):
		val = 'X'
		try:
			val = list( self.pdbRes ).index( resNum )
		except ValueError:
			print ("Requested pdb residue %d not in chain %s" % (resNum, self.name) )
			pass

		try:
			val = self.uniRes[val]
		except IndexError:
			print  ("Requested pdb residue %d not in uniprot mapping of chain %s" % (val, self.name) )
			pass
		
		return val


	def uni2pdb(self, resNum):
		val = 'X'
		try:
			val = list( self.uniRes ).index( resNum )
		except ValueError:
			print ( "Requested uniprot residue %d not in uniprot entry %s" % (resNum, self.uniprot) )
			pass

		try:
			val = self.pdbRes[val]
		except IndexError:
			print  ( "Requested uniprot residue %d not in pdb mapping of chain %s" % (val, self.name) )
			pass
		
		return val


########## Metal pair class  ##########
#Initialize with prody atom objects for metals and list of prody residue objects contacting
# Tihs class slims down info each atom carries to only store essentials (reduce filesize/RAM used)
class biMsite():
	def __init__(self, atom1, atom2, residue_matesList, pdb):
		self.name = "%s_%s+%s_%s" % ( atom1.getElement(), str( atom1.getSerial() ), atom2.getElement(), str( atom2.getSerial() ) )
		self.pdb = pdb
		self.contacts = []
		#For coordinating residues, loss of metal-res link, but bidentate should be store as 'A_ASP180-OD2+OD1'
		#self.contacts = [ '%s_%s%s-%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName() ) for mate in residue_matesList ] 
		for mate in residue_matesList:
			biFlg = 0
			fullID 	= '%s_%s%s-%s=%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ), mate.getName(), str(int(mate.getBeta()) ) )	 
			if len( self.contacts ) == 0:
				self.contacts.append( fullID )
			else:
				rBaseId = '%s_%s%s' % ( mate.getChid(), mate.getResname(), str( mate.getResnum() ) ) 

				for c in self.contacts:
					if c.split('-')[0] == rBaseId:
						new = c.split('=')[0] + '+' + mate.getName() + '=' + str(int(mate.getBeta()) ) + '=' + c.split('=')[-1]
						old = c
						biFlg +=1
				if biFlg > 0 :
					self.contacts.remove( old )
					self.contacts.append( new )
				else:
					self.contacts.append( fullID )


				

	def __repr__(self):
		return self.name




###### class end ######################



####################### Class definitions end #####################


############### Functions  ########################################

## Return the count of each amino acid found (includes non-natural) in a PDB model
def cntAA( pdbPath ):
	aaList = []
	print ('*')
	pdbF = parsePDB( pdbPath )
	print ('*')

	for res in pdbF.iterResidues():
		try:
			name = unNatAA[ res.getResname() ]
		except KeyError:
			name = res.getResname()
		if name in natAA.keys():
			aaList.append( name )

	final = Counter ( aaList )	

	return final

# determine the frequency of each amino acids in percent (e.g. 15.6) for a pdb database, given a list file with a path to each file
def freqAA (pathListFile):
	freq = Counter([])

	with open( pathListFile ) as file:
		for i in file:
			freq.update(   cntAA( i.rstrip() ) )
			print (freq.items(), '\n')

	freq2 = {}		
	for k,v in freq.items():
		val = 100 * float(v) / sum( freq.values() )
		freq2[k] = round(val, 2)

	print (freq2.items(), '\n', sum( freq2.values() )/100.0 )
	return freq2

#freqAA( 'localZNdbFiles.txt' )


# Structural clustering using pre-user-defined symmetric distance matrix  
# Input a threshold level (presumably RMSD-cut off) such that all clusters are within threshold value
def heirarchy_cluster( matrix, threshold_RMSD = 2.0):
	import scipy.cluster.hierarchy as sp_clust

	# Complete linkage clustering
	link_matrix = sp_clust.complete( matrix )
	# Make flat clusters at the tree point where distance threshold is met: default RMSD < 2.0
	clusters 	= sp_clust.fcluster( link_matrix, threshold_RMSD, criterion = 'distance' )

	return clusters


def kmedoid_clustering( matrix, threshold_RMSD, lookupHash):
	from Bio.Cluster.cluster import kmedoids
	from collections import Counter
	import time

	cost = 4.0
	c = 135
	while cost > threshold_RMSD:
		trials	= np.size( matrix )
		start 	= time.time()
		clust 	= kmedoids( matrix , c, trials)
		print (time.time()-start, 'elapsed for ', c, 'clusters', clust[-1], 'identical trajectory(s) from total',trials,)
		# for this clustering scheme, calculate mean RMSD of each cluster
		clusterStats = {}
		index = 0
		# Look for each RMSD to the cluster centroid, store max within the cluster
		for p in clust[0]:
			#print p, index, matrix[p][index]
			try:
				clusterStats[p] = max( matrix[p][index], clusterStats[p] )
			except KeyError:
				clusterStats[p] = matrix[p][index]

			index += 1

		cost = max( clusterStats.values() )*4.0
		print ('maximum rmsd: ', cost)
	print (Counter( clust[0] ))
	return clust


## extend an alpha helix both C-terminally and N terminally, given the 4 pro/pre-ceding residue
## Uses GCN4 average internal coordinates of a 10 residue ideal helix
## input prody atom group of object to extend and atom group of ideal GCN4 mimic helix
## Some of the ends of the helix with non-ideal helical dihedrals should get minimized in rosetta later if necessary
def ext_aHelixBOTH ( obj, idealH ):

	# N term
	chainObj 	= obj.getChids()[0]
	first4		= obj.select( 'calpha' ).getResnums()[:4]
	selStrTarg 	= 'chain %s resnum %s' % ( chainObj, ' '.join( [str(x) for x in first4] ) )
	selStrMob	= 'chain A resnum 6 7 8 9'
	objNtarg 	= obj.select( selStrTarg )
	idealHmobN	= idealH.select( selStrMob ).copy()
	transN		= superpose( idealHmobN, objNtarg )[1]
	# fit chain and resnums to match those 
	mobileN 	= idealH.select('chain A resnum 2 3 4 5').copy()
	mobileN.setChids( [ 'X' for x in mobileN.getChids() ] )
	preResNums	= [ x - 4  for x in first4 ]
	for r, new in zip(  mobileN.iterResidues(), preResNums ):
		r.setResnum( new )
	applyTransformation( transN, mobileN )
	obj = mobileN + obj

	# C term
	last4		= obj.select( 'calpha' ).getResnums()[-4:]
	selStrTarg 	= 'chain %s resnum %s' % ( chainObj, ' '.join( [str(x) for x in last4] ) )
	selStrMob	= 'chain A resnum 2 3 4 5'
	objCtarg 	= obj.select( selStrTarg )
	idealHmobC	= idealH.select( selStrMob ).copy()
	transC		= superpose( idealHmobC, objCtarg )[1]
	# fit chain and resnums to match those 
	mobileC 	= idealH.select('chain A resnum 6 7 8 9').copy()
	mobileC.setChids( [ 'X' for x in mobileC.getChids() ] )
	postResNums	= [ x + 4  for x in last4 ]
	for r, new in zip(  mobileC.iterResidues(), postResNums ):
		r.setResnum( new )
	applyTransformation( transC, mobileC )
	obj = obj + mobileC


	return obj


def random_membrane_res():
	return ['A', 'F', 'V', 'I', 'L', 'L', 'L', 'L', 'L', 'L'][randint(0,9)]


### ONLY CAN RUN THIS FROM INSIDE FATCAT DIRECTORY
def calcAlign( path1, path2 ):

	cmd = [ 'bash', 'runCE.sh', '-file1', path1, '-file2', path2, '-printCE' ]

	output = sp.Popen( cmd, stdout=sp.PIPE, stderr=sp.PIPE )

	match = re.search(r'Rmsd = (\d+\.\d+)A', output.communicate()[0] )
	if match:
		return match.group(1)
	else:
		print ('ERROR IN FINDING ALIGNMENT MATCH... exiting')
		sys.exit()



