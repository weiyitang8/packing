from Bio.PDB import * 
import numpy as np,sys,os,pandas as pd
import itertools
# import matplotlib.pyplot as plt
import sasa_sourcecode as sasa_src


def bbfull(atomlist):

	bbfullatoms=[atom for atom in atomlist if atom.name=="CA" or atom.name=="C"\
	  or atom.name=="N"  or atom.name=="O"  or atom.name=="H"  or atom.name=="HA" ]

	return bbfullatoms

#################################find holes###########################

def findholes(chainset):

	chainholeset=[]

	for chain in chainset:
		holenum =len(chain)-7
		holeset =list(range(holenum))
		resset = chain.get_list()
		
		for i in np.arange(holenum):
			holeset[i]=[[resset[i+3],resset[i+4]],[resset[i],resset[i+7]]]
		holeset.insert(0,[[resset[1],resset[2]],[resset[5],0]])
		holeset.insert(1,[[resset[2],resset[3]],[resset[6],0]])
		holeset.insert(2,[[resset[3],resset[4]],[resset[7],0]])

		holeset.append([[resset[-4],resset[-3]],[0,resset[-7]]])
		holeset.append([[resset[-3],resset[-2]],[0,resset[-6]]])
		holeset.append([[resset[-2],resset[-1]],[0,resset[-5]]])

		chainholeset.append(holeset)
		

	return chainholeset


def getholedir(hole):

	holevec34=hole[0][1]['CA'].get_coord()-hole[0][0]['CA'].get_coord()
	cen=(hole[0][1]['CA'].get_coord()+hole[0][0]['CA'].get_coord())/2
	if hole[1][1] == 0:
		holevec07=hole[1][0]['CA'].get_coord()-cen
	elif hole[1][0] ==0:
		holevec07=cen-hole[1][1]['CA'].get_coord()
	else:
		holevec07=hole[1][1]['CA'].get_coord()-hole[1][0]['CA'].get_coord()	
	holevec=np.cross(holevec34,holevec07)
	unitholevec=holevec/np.linalg.norm(holevec)

	return unitholevec


def norm_smash(hole,vec):
	vec=vec/np.linalg.norm(vec)
	vec[2]=0
	hole.dir=vec
	return vec

def sasahole(hole):


	bbatoms=getholebbfull(hole)
	sasa_hole=[]
	for res in bbatoms:
		bbatoms       =bbfull(res)
		sasa4hole     =(sum([atom.sasa.round(3) for atom in res])).round(3)
		sasa_hole.append(sasa4hole)
	return sasa_hole

def getknobname(knob):

	namelist=['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in knob for knob in chain]

	return namelist

def getholebbfull(hole):
	if hole[1][1] == 0:
		holebbfull=[bbfull(res) for res in hole[0]]+[bbfull(hole[1][0])]
	elif hole[1][0] ==0:
		holebbfull=[bbfull(res) for res in hole[0]]+[bbfull(hole[1][1])]
	else:
		holebbfull = [bbfull(res) for res in hole[0]]+[bbfull(res) for res in hole[1]]

	return holebbfull

#################################bbfull#################################
class Hole:

	def __init__(self,holes):
		self.hole = holes
		self.eli  = 1
		self.knob = []
		self.knobsc = 0
		self.score=0

		self.id   = 'hole'+holes[0][0].get_full_id()[2]+str(holes[0][0].get_id()[1])+'-'+str(holes[0][1].get_id()[1])
		# print(self.id)
		# self.elisc  = 0
		# self.prepack_knobsc =0
		# self.d=[]
		# self.eknobsc =0
		# self.prepack_knob=[]
		# self.distance = []
		# self.knob=[]
		# self.knobname=[]
		# self.dsasa_knob=[]

	def __rper__(self):
		return self.id

	
	def npsasasc(self):

		sasa_hole=sasahole(self.hole)

		self.npsasals  = sasa_hole
		self.npsasa    = sum(sasa_hole).round(3)
	
	def psasasc(self):

		sasa_hole=sasahole(self.hole)

		self.psasals  = sasa_hole
		self.psasa    = sum(sasa_hole).round(3)

	def dsasasc(self):
		
		self.dsasa    = (self.npsasa-self.psasa).round(3)
		self.dsasals  = (np.array(self.npsasals)-np.array(self.psasals)).round(3)
 	

	def patchsasasc(self):
		
		sasa_hole=sasahole(self.hole)

		self.patchsasals  = sasa_hole
		self.patchsasa    = sum(sasa_hole).round(3)

	def patch1(self):
		dsasa34=np.array([dsasa for dsasa in self.dsasals[1:3]])
		# print(dsasa34)

		self.patchdsasa= dsasa34

		if (dsasa34 == 0).any():
			self.knob=[]
			self.knobsc=0




	#########################vector method###############


	def facingholes_AB(self,unit_packvec,unit_prepackvec=[]):

		
		holevec    = getholedir(self.hole)
		holevec[2] = 0
		facingproduct = np.dot(holevec,unit_packvec)

		#self.prepack = 1
		#self.facingx = 1

		if not facingproduct > 0:
			self.eli = 0
			#self.facingx -= 1

		if unit_prepackvec != []:       #prepack if prepack can happend
			prepackproduct = np.array([np.dot(holevec,p) for p in unit_prepackvec])

			if (prepackproduct>0.75).any():
				self.eli = 0
				# self.elisc=0
				#self.prepack -=1



	def findknobs(self,chainX,contactparm):



		holeCA34       = [res['CA'] for res in self.hole[0]]
		holeCA34coord  = np.array([ca.get_coord() for ca in holeCA34])
		# holebbfull   = [bbfull(res) for res in self.hole]
		knoblist     = []
		
		for knob in chainX: 
			knobCAcoord  = knob['CA'].get_coord()
			distance     = np.array([holeCA34coord-knobCAcoord]).round(3)
			distance     = np.linalg.norm(distance[0],axis=1)


			if (distance<12).all():
				knoblist.append(knob)


		if len(knoblist) > 0:
			for knob in knoblist[:]:
				knob_sc  = knob.get_list()[4:]
				knob_sch = [atom  for atom in knob_sc if atom.element != 'H']
				if knob_sch !=[]:
					knobCA        = knob['CA']			
					dis_knobatom  = np.array([knobscatom-knobCA for knobscatom in knob_sch])
					faratom       = knob_sch[np.argmax(dis_knobatom)]
					knobvec       = faratom.get_coord()-knobCA.get_coord()
					unitknobvec   = knobvec/np.linalg.norm(knobvec)
					unitholevec   = self.dir
					dotproduct    = np.dot(unitholevec,unitknobvec)

					if dotproduct > 0:
						knoblist.remove(knob)

				else:
					knoblist.remove(knob)



		if len(knoblist) >0:
			holebbfull=getholebbfull(self.hole)
			for knob in knoblist[:]:
				knob_sc=knob.get_list()[4:]
				a=0
				for res in holebbfull:
					for bbatom in res:
						distance = np.array([bbatom - knob_scatom for knob_scatom in knob_sc])
						distance = np.min(distance)
						if distance < contactparm:
							a+=1
				if a == 0:
					knoblist.remove(knob)

		self.knob.append(knoblist)
		self.knobsc+=len(knoblist)




##########################################################################################################



def excute(protein,chainABC_str,chainX_str,x=0,y=0,contactparm2=3): #chainABC_str has to be 2 letters with NO spacing in between, chainX_str has to be one letter, x is the cut-off for the final score.

	chainABC  = [protein[a] for a in chainABC_str]
	chainA    = chainABC[0]
	chainB    = chainABC[1]
	chainX    = protein[chainX_str]
	

	A_atoms = list(chainA.get_atoms())
	B_atoms = list(chainB.get_atoms())
	X_atoms = list(protein[chainX_str].get_atoms())


	chainset=[chainA,chainB,chainX]
	chain_twinset=list(itertools.combinations(chainset,2))

	geometry={}
	for chain_twin in chain_twinset:
		dis=chain_twin[1].get_list()[0]['CA']-chain_twin[0].get_list()[0]['CA']
		# print(chain_twin,dis,chain_twin[1].get_list()[0],chain_twin[1].get_list()[0].get_parent(),chain_twin[0].get_list()[0])
		# print(chainA.get_list())
		if dis >20:

			geometry[chain_twin]='antiparellel'
			# print('anti')
		else:
			geometry[chain_twin]='parellel'
			# print('pp')



	ABgeo=geometry[chain_twinset[0]]
	AXgeo=geometry[chain_twinset[-1]]

	if ABgeo =='parellel' and AXgeo =='antiparellel':
		topA=np.array([res['CA'].coord for res in chainA.get_list()[:4]])
		topB=np.array([res['CA'].coord for res in chainB.get_list()[:4]])
		topX=np.array([res['CA'].coord for res in chainX.get_list()[-4:]])

		botA=np.array([res['CA'].coord for res in chainA.get_list()[-4:]])
		botB=np.array([res['CA'].coord for res in chainB.get_list()[-4:]])
		botX=np.array([res['CA'].coord for res in chainX.get_list()[:4]])


	elif ABgeo =='antiparellel' and AXgeo =='parellel':
		topA=np.array([res['CA'].coord for res in chainA.get_list()[:4]])
		topB=np.array([res['CA'].coord for res in chainB.get_list()[-4:]])
		topX=np.array([res['CA'].coord for res in chainX.get_list()[:4]])

		botA=np.array([res['CA'].coord for res in chainA.get_list[-4:]])
		botB=np.array([res['CA'].coord for res in chainB.get_list[:4]])
		botX=np.array([res['CA'].coord for res in chainX.get_list[-4:]])


	top_center  = np.array([topA,topB,topX]).mean(axis=0).mean(axis=0)
	bot_center  = np.array([botA,botB,botX]).mean(axis=0).mean(axis=0)
	center      = (top_center+bot_center)/2
	dis         = np.linalg.norm(top_center-center)

	p,q=vectors.Vector(0,0,-dis),vectors.Vector(top_center-center)

	protein.transform(rotmat(p,q),np.array(-center,'f'))


	step_size     = [len(chain)//3 for chain in chainABC]
	chainlistset  = [chain.get_list() for chain in chainABC]

	chain_3parts  = list(range(len(chainlistset)))
	centerset     = list(range(len(chainlistset)))

	for i,chain in enumerate(chainlistset):
		top=chain[:step_size[i]]     
		mid=chain[step_size[i]-1:2*step_size[i]+1]
		bot=chain[2*step_size[i]:]

		# print([res.id[1] for res in top])
		# print([res.id[1] for res in mid])
		# print([res.id[1] for res in bot])

		# top=chain[:step_size[i]]     
		# mid=chain[step_size[i]:2*step_size[i]]
		# bot=chain[2*step_size[i]:]

		chain_3parts[i] =[top,mid,bot]
		centerset[i]    =[0,0,0]
		for j,part in enumerate(chain_3parts[i]):
			coordset        = np.array([res['CA'].coord for res in part])
			centerset[i][j] = coordset.mean(axis=0).round(3)
	centerset=np.array(centerset)
	# print(centerset)
	           


	if ABgeo == 'parellel':
		first_plane  =[]
		second_plane =[]
		for chain in centerset:
			first_plane.append((chain[0][-1]+chain[1][-1])/2)
			second_plane.append((chain[2][-1]+chain[1][-1])/2)
			
		first_plane=np.array(first_plane).mean().round(3)
		second_plane=np.array(second_plane).mean().round(3)

		# print(first_plane)
		# print(second_plane)

		Xtop=[]
		Xmid=[]
		Xbot=[]
		# chainXlist=chainX.get_list()
		for res in chainX.get_list():
			coord=res['CA'].get_coord()[2]
			# print(res.get_id(),coord)
			if coord < first_plane:
				Xtop.append(res)
			elif coord > second_plane:
				Xbot.append(res)
			else:
				Xmid.append(res)

		# Xtop.append(chainX[Xtop[0].get_id()[1]-2])     ####weird logic on list sort
		# Xtop.append(chainX[Xtop[0].get_id()[1]-1])	###better not touch anything here
		Xmid.append(chainX[Xmid[-1].get_id()[1]+1])
		Xmid.append(chainX[Xmid[0].get_id()[1]-1])
		# Xbot.append(chainX[Xbot[-1].get_id()[1]+1])
		# Xbot.append(chainX[Xbot[-1].get_id()[1]+1])
		# Xtop.sort()
		Xmid.sort()

		# print([res.id[1] for res in Xtop])
		# print([res.id[1] for res in Xmid])
		# print([res.id[1] for res in Xbot])


		Xchain_3parts=[Xtop,Xmid,Xbot]
		Xcenterset=[0,0,0]
		# print("Xtopmidbot",Xchain_3parts)
		for i,part in enumerate(Xchain_3parts):
			coordset=np.array([res['CA'].get_coord().round(3) for res in part])
			Xcenterset[i]=coordset.mean(axis=0).round(3)

	Xcenterset=np.array(Xcenterset)


	packvecset=centerset-Xcenterset
	unit_packvecset=packvecset

	for i,chain in enumerate(packvecset):
		for j,part in enumerate(chain):
			unit_packvecset[i][j]=(part/np.linalg.norm(part)).round(3)
			unit_packvecset[i][j][2]=0


	centerset_twinset=list(itertools.permutations(centerset,2))

	prepackvecset=list(range(len(centerset_twinset)))
	for i,centerset_twin in enumerate(centerset_twinset):
		prepackvecset[i]=centerset_twin[1]-centerset_twin[0]

	unit_prepackvecset=list(range(len(centerset_twinset)))
	for i,twin in enumerate(prepackvecset):
		unit_prepackvec=[0,0,0]
		for j,part in enumerate(twin):
			unit_prepackvec[j]=(part/np.linalg.norm(part)).round(3)
			unit_prepackvec[j][2]=0
		unit_prepackvecset[i]=unit_prepackvec

#########################################################
	

	ABCholeset = findholes(chainABC)
	Xholeset   = findholes([chainX])
	aholes     = ABCholeset
	xholes     = Xholeset

	for i,chain in enumerate(ABCholeset):		
		for j,hole in enumerate(chain):
			aholes[i][j] = Hole(hole)

	for i, chain in enumerate(Xholeset):
		for j,hole in enumerate(chain):
			xholes[i][j] = Hole(hole)

	topid=[res.get_id()[1] for res in top]
	midid=[res.get_id()[1] for res in mid]
	botid=[res.get_id()[1] for res in bot]

	Xtopid=[res.get_id()[1] for res in Xtop]
	Xmidid=[res.get_id()[1] for res in Xmid]
	Xbotid=[res.get_id()[1] for res in Xbot]

	# print('############################')


	for i,chain in enumerate(aholes):
		for j,hole in enumerate(chain):
			holevec=getholedir(hole.hole)
			holevec[2]=0
			hole.dir=holevec.round(3)
			holeid=[res.get_id()[1] for res in hole.hole[0]]

			if holeid[0] in topid and holeid[1] in topid:
				facingproduct = np.dot(holevec,unit_packvecset[i][0]).round(3)
				prepackproduct = np.dot(holevec,unit_prepackvecset[i][0]).round(3)
				# print(facingproduct)

			elif holeid[0] in midid and holeid[1] in midid:
				facingproduct = np.dot(holevec,unit_packvecset[i][1]).round(3)
				prepackproduct = np.dot(holevec,unit_prepackvecset[i][1]).round(3)
				# print(facingproduct)


			elif holeid[0] in botid and holeid[1] in botid:
				facingproduct = np.dot(holevec,unit_packvecset[i][2]).round(3)
				prepackproduct = np.dot(holevec,unit_prepackvecset[i][2]).round(3)
				# print(facingproduct)


			
			if not facingproduct <0:
				hole.eli=0
			if prepackproduct >0.75:
				hole.eli=0
			# if hole.eli ==1 :
			# 	print(hole.id,hole.eli,facingproduct,prepackproduct)





	for i, chain in enumerate(xholes):
		for j,hole in enumerate(chain):
			holevec=getholedir(hole.hole)
			holevec[2]=0
			hole.dir=holevec
			holeid=[res.get_id()[1] for res in hole.hole[0]]
			# print(holeid,holevec)
			# print(Xtop)

			# holeCAcoord=np.array([res['CA'].coord for res in hole.hole])
			# holemeancoord= holeCAcoord.mean(axis=0)
			if holeid[0] in Xtopid and holeid[1] in Xtopid:		
				# holevec=holemeancoord-Xcenterset[0]
				# holevec=norm_smash(hole,holevec)
				facingproductA = np.dot(holevec,unit_packvecset[0][0]).round(3)
				facingproductB = np.dot(holevec,unit_packvecset[1][0]).round(3)
				# print(facingproductA)
				# print(facingproductB)


			elif holeid[0] in Xmidid and holeid[1] in Xmidid:
				facingproductA = np.dot(holevec,unit_packvecset[0][1]).round(3)
				facingproductB = np.dot(holevec,unit_packvecset[1][1]).round(3)
				# print(facingproductA)
				# print(facingproductB)


			elif holeid[0] in Xbotid and holeid[1] in Xbotid:

				facingproductA = np.dot(holevec,unit_packvecset[0][2]).round(3)
				facingproductB = np.dot(holevec,unit_packvecset[1][2]).round(3)
				# print(facingproductA)				
				# print(facingproductB)


			if not (facingproductA >0 and facingproductB >0):
				hole.eli=0

			# if hole.eli ==1 :
			# 	print(hole.id,hole.eli,facingproductA,facingproductB)

	###########################################knobfinder##########################

	sum1=0
	knoblist=[]
	for chain in aholes:
		for hole in chain:
			if hole.eli == 1:
				hole.findknobs(chainX,contactparm2)
				if hole.knobsc>0:
					# print(hole.id,hole.eli,prepackproduct, hole.knobsc,['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in hole.knob for knob in chain])
					sum1 = hole.knobsc+sum1
					# knoblist.append(hole.knob)


	sum2=0
	for chain in xholes:
		for hole in chain:
			if hole.eli== 1:
				hole.findknobs(chainA,contactparm2)
				hole.findknobs(chainB,contactparm2)
				if hole.knobsc>0:
					# print(hole.id,hole.eli,prepackproduct, hole.knobsc,['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in hole.knob for knob in chain])

					sum2 = hole.knobsc+sum2

	# print('######',sum1,sum2,sum1+sum2)
	

	#################################SASA#####################

   #filter to process the protein with SASA method

	a_bb_atoms      =bbfull(A_atoms)
	b_bb_atoms      =bbfull(B_atoms)
	x_bb_atoms	    =bbfull(X_atoms)		
	abc_bb_atoms    =a_bb_atoms+b_bb_atoms
	abc_atoms       =A_atoms+B_atoms

	abcbb_x_atoms	= abc_bb_atoms + X_atoms
	xbb_abc_atoms	= abc_atoms + x_bb_atoms


	sasa=sasa_src.ShrakeRupley_nsc(probe_radius=1.2)   #ABCchain as hole chain
	sasa.compute(abc_bb_atoms)
	for chain in aholes:
		for hole in chain:	
			if hole.eli ==1:		
				sasa_hole      = sasahole(hole.hole)
				hole.npsasals  = sasa_hole
				hole.npsasa    = sum(sasa_hole).round(3)


	sasa.compute(abcbb_x_atoms)
	for chain in aholes:
		for hole in chain:
			if hole.eli ==1:
				sasa_hole	  = sasahole(hole.hole)
				hole.psasals  = sasa_hole
				hole.psasa    = sum(sasa_hole).round(3)
				hole.dsasa    = (hole.npsasa-hole.psasa).round(3)
				hole.dsasals  = (np.array(hole.npsasals)-np.array(hole.psasals)).round(3)
				# dsasa34		  =  hole.dsasals[1:3]
				# hole.patchdsasa= dsasa34
				if (hole.dsasals[0:2] == 0).any():
					# print('dsasa',hole.id,hole.dsasals)
					hole.knob=[]
					hole.knobsc=0

  
	sasa.compute(x_bb_atoms) # Xchain a hole chain
	for chain in xholes:
		for hole in chain:	
			if hole.eli == 1:		
				sasa_hole      = sasahole(hole.hole)
				hole.npsasals  = sasa_hole
				hole.npsasa    = sum(sasa_hole).round(3)


	sasa.compute(xbb_abc_atoms)
	for chain in xholes:
		for hole in chain:
			if hole.eli ==1:
				sasa_hole	  = sasahole(hole.hole)
				hole.psasals  = sasa_hole
				hole.psasa    = sum(sasa_hole).round(3)
				hole.dsasa    = (hole.npsasa-hole.psasa).round(3)
				hole.dsasals  = (np.array(hole.npsasals)-np.array(hole.psasals)).round(3)
				# dsasa34		  = np.array([dsasa for dsasa in hole.dsasals[1:3]])
				# hole.patchdsasa= dsasa34
				if ( hole.dsasals[0:2]== 0).any():
					# print('dsasa34',hole.id)
					hole.knob=[]
					hole.knobsc=0
	# print('##################################################')
	sum3=0
	for chain in aholes:
		for hole in chain:
			if hole.knobsc > 0:
				# print(hole.id, hole.eli, hole.knobsc, ['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in hole.knob for knob in chain])
				sum3+=hole.knobsc

	sum4=0
	for chain in xholes:
		for hole in chain:
			if hole.knobsc > 0:
				# print(hole.id, hole.eli, hole.knobsc, ['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in hole.knob for knob in chain])
				sum4+=hole.knobsc

	# print(sum3,sum4,sum3+sum4)

	# if (sum3+sum4)>=13:
	# 	print('Congrats, this design pass the first 2 filters')


	# ##################################################
	# 3rd filter, delta SASA befero/after the each knob packing into the hole chain
	
	sum5=0
	sum6=0
	finalscore=0
	if sum3+sum4 > y:


		if len(chainABC) >1 :
			bb_atoms 		=[a_bb_atoms,b_bb_atoms]
			abb_bfull       =a_bb_atoms + B_atoms
			bbb_afull       =A_atoms    + b_bb_atoms
			prepackatoms    =[abb_bfull , bbb_afull]
			xbb_afull		=x_bb_atoms + A_atoms
			xbb_bfull		=x_bb_atoms + B_atoms
			xprepackatoms	=[xbb_bfull , xbb_afull]



		for i,chain in enumerate(aholes):
			sasa.compute(prepackatoms[i])
			for hole in chain:
				if hole.knobsc >0:
					sasa_hole           = sasahole(hole.hole)
					hole.prepacksasals  = np.array(sasa_hole)


		for i, chain in enumerate(aholes):
			for hole in chain:
				if hole.knobsc >0:
					for chain in hole.knob[:]:
						for knob in chain[:]:
							knobpackatoms  = prepackatoms[i]+list(knob.get_atoms())
							sasa.compute(knobpackatoms)
							sasa_hole      = np.array(sasahole(hole.hole))
							ppknob         = (hole.prepacksasals-sasa_hole).round(3)
							ppknob34	   = ppknob[0:2]
							# print(hole.id,ppknob34,ppknob)
							if (ppknob34==0).any():
								chain.remove(knob)
								hole.knobsc -=1
							elif (ppknob34<1).any():
								hole.score +=0.5
							else:
								hole.score +=1
						if len(chain) == 0:
							hole.knob.remove(chain)
		


		for i,chain in enumerate(xholes):
			sasa.compute(xprepackatoms[i])
			for hole in chain:
				if hole.knobsc >0:
					sasa_hole           = sasahole(hole.hole)
					hole.prepacksasals  = np.array(sasa_hole)


		for i, chain in enumerate(xholes):
			for hole in chain:
				if hole.knobsc >0:
					for chain in hole.knob[:]:
						for knob in chain[:]:
							knobpackatoms  = xprepackatoms[i]+list(knob.get_atoms())
							sasa.compute(knobpackatoms)
							sasa_hole      = np.array(sasahole(hole.hole))
							ppknob         = (hole.prepacksasals-sasa_hole).round(3)
							ppknob34	   = ppknob[0:2]
							# print(hole.id,ppknob34,ppknob,hole.prepacksasals,sasa_hole)
							if (ppknob34==0).any():
								chain.remove(knob)
								hole.knobsc -=1
							elif (ppknob34<1).any():
								hole.score +=0.5
							else:
								hole.score +=1
						if len(chain) == 0:
							hole.knob.remove(chain)


		
		for chain in aholes:
			for hole in chain:
				if hole.knobsc>0:
					finalscore+=hole.score
					sum5+=hole.knobsc
					# print('For',hole.id,',',hole.knobsc,'knob(s) packing into the hole,', getknobname(hole.knob))
					#

		
		for chain in xholes:
			for hole in chain:
				if hole.knobsc>0:
					finalscore+=hole.score
					sum6+=hole.knobsc
					# print('For',hole.id,',',hole.knobsc,'knob(s) packing into the hole,', getknobname(hole.knob))

					# print(hole.id, hole.score,hole.knobsc,['knob'+str(knob.get_full_id()[2])+str(knob.get_id()[1]) for chain in hole.knob for knob in chain])

	result=[]
	if finalscore >x:
		holes=[aholes,xholes]
		for nholes in holes:
			for chain in nholes:
				for hole in chain:
					if hole.knobsc>0:
						# print('x')
						result.append([hole.id,getknobname(hole.knob),hole.score])

		# if finalscore >x:
		# 	print('*** good design,finalscore is above',x,'   ###', protein.get_parent().id)
		# 	print('*** The finalscore is ', finalscore, 'there are',sum6+sum5, 'knobs packing into holes')
		# 	holes=[aholes,xholes]
		# 	for nholes in holes:
		# 		for chain in nholes:
		# 			for hole in chain:
		# 				if hole.knobsc>0:
		# 					print(hole.id,getknobname(hole.knob))



	return sum5,sum6,finalscore,result


