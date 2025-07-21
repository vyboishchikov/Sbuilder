# (c) Sergei Vyboishchikov March-April 2024
import numpy as np
import sys, os, math, random
import networkx as nx


try:
	from rdkit import Chem
	from rdkit.Chem import AllChem
	RDKit_available = True
	print("RDKit is available")
except ImportError:
	RDKit_available = False
	print("RDKit is not installed. Some functionality may be unavailable.")

def MainProgram():
	global AtomicNamesDict, CovalentRadius, AtomicNames, N, NGen, AdjMatrixRequested, BackupRequested, RestoreRequested, Randomi, RandomizedGen, CheckStericsRequested, RunRDKIT, OptimizerRDKIT, ForceField
	NumberOfGenerations = 1    # default number of generations to be created
	InputLabels = []
	RestoreRequested = False
	CheckStericsRequested = False
	RunRDKIT = ""
	ResultFile = "Result.xyz"  # default name of the resulting xyz file
	CurrentFileName = os.path.basename(sys.argv[0])
	print(CurrentFileName+' '+' '.join(sys.argv[1:]))
	if len(sys.argv)<2:
		print("Usage:\n"+CurrentFileName+" Starting.xyz -gen number_of_generations -out Output.xyz  [-random gen1:share1[,gen2:share2]] [-opt uff|mmff] [-adjmatrix] [-backup] [-restore gen1,gen2,...]")
		print("or     \n"+CurrentFileName+" -label label -out Output.xyz")
		quit()
	for i in range(len(sys.argv)-1):
		if sys.argv[i]=="-gen":
			NumberOfGenerations = int(sys.argv[i+1])
			Randomi = {str(r):False for r in range(NumberOfGenerations+1)}
		elif sys.argv[i]=="-out": ResultFile = sys.argv[i+1]
		elif sys.argv[i]=="-label" or sys.argv[i]=="-labels": InputLabels = sys.argv[i+1].split(",")
		elif sys.argv[i]=="-restore":
			RestoreRequested = True
			RestoreGenerations = sys.argv[i+1].split(",")
		elif sys.argv[i]=="-opt":
			RunRDKIT = sys.argv[i+1]
			if RunRDKIT.lower() == 'uff':
				OptimizerRDKIT = AllChem.UFFOptimizeMolecule
				ForceField = "UFF-"
			elif RunRDKIT.lower() == 'mmff':
				OptimizerRDKIT = AllChem.MMFFOptimizeMolecule
				ForceField = "MMFF-"
			else:
				OptimizerRDKIT = None
				print("Warning: force field",RunRDKIT," is not understood")
				ForceField = "not "
		elif sys.argv[i]=="-random":
			RandomizedGen = {}
			for r in sys.argv[i+1].split(","):
				RandomizedGen[r.split(":")[0]] = float(r.split(":")[1]) 
				Randomi[r.split(":")[0]] = True
	AdjMatrixRequested = False
	BackupRequested = False
	for i in range(len(sys.argv)):
		AdjMatrixRequested = AdjMatrixRequested or sys.argv[i]=="-adjmatrix" or sys.argv[i]=="-adj_matrix"
		BackupRequested = BackupRequested or sys.argv[i]=="-backup"
		CheckStericsRequested = CheckStericsRequested or sys.argv[i]=="-checksterics" or sys.argv[i]=="-check_sterics"
		RunRDKIT = RunRDKIT or sys.argv[i]=="-rdkit"
	CovalentRadius=[0.0,0.38,0.32,1.34,0.9,0.82,0.77,0.75,0.73,0.71,0.69,1.54,1.3,1.18,1.11,1.06,1.02,0.99,0.97]
	CovalentRadius+=[2.27,1.97,1.62,1.47,1.34,1.28,1.27,1.26,1.26,1.24,1.38,1.34,1.35,1.22,1.19,1.16,1.14,1.1]
	CovalentRadius+=[2.48,2.15,1.8,1.6,1.46,1.39,1.36,1.34,1.35,1.37,1.53,1.51,1.67,1.41,1.38,1.35,1.33,1.3]
	CovalentRadius+=[2.65,2.22,1.87,1.818,1.824,1.814,1.834,1.804,1.804,1.804,1.773,1.781,1.762,1.761,1.759]
	CovalentRadius+=[1.76,1.738,1.59,1.46,1.46,1.59,1.35,1.37,1.385,1.44,1.51,1.7,1.47,1.46,1.29,1.38,1.45]
	CovalentRadius+=[1.53,1.59,1.4,1.79,1.63,1.56,1.55,1.59,1.73,1.74] + [1.75 for i in range(25)]

	AtomicNames = ['X','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',	'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',	'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',	'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',	'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',	'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds']
	AtomicNamesDict = {r:i for i,r in enumerate(AtomicNames)}
	for i,r in enumerate(AtomicNames):
		AtomicNamesDict[str(i)] = i
	if len(InputLabels)>0: # single-point calculation
		Titol = InputLabels[0].split(":")[0]
		XYZinput = Titol+".xyz"
		print("Input file: "+XYZinput)
		InitialGeometry = ReadXYZFile(XYZinput)
		StartingGeometry = MakeGrafFromGeometry(InitialGeometry,Titol)
		print("Starting geometry:")
		Fitxer = open(ResultFile,"w")
		PrintGeom(StartingGeometry)
		for InputLabel in InputLabels:
			print("Single-point mode: structure generation from label",InputLabel)
			CarbonPairs = [[InsertColon(r.split("-")[0]),InsertColon(r.split("-")[1])] for r in InputLabel.split(":")[1:]]
			Geometry = StartingGeometry.copy()
			print("CarbonPairs=",CarbonPairs)
			for i,carpa in enumerate(CarbonPairs):
				NGen = str(i+1)
				BenzeneGraf = MakeBenzeneGraf(NGen)
				r = FindHydrogens(Geometry,carpa)
				NumC1,NumC2,NumH1,NumH2 = r
				NewGraf = nx.union(Geometry,BenzeneGraf)
				BenzeneAtoms = ['C:1g'+NGen,'C:2g'+NGen,'C:3g'+NGen,'C:4g'+NGen,'H:1g'+NGen,'H:2g'+NGen,'H:3g'+NGen,'H:4g'+NGen]
				BenzeneXYZ = MakeBenzeneXYZ(Geometry,r)
				NewGraf.remove_nodes_from([NumH1,NumH2])
				NewGraf.add_edges_from([(NumC1,'C:4g'+NGen), (NumC2,'C:1g'+NGen)])
				nx.set_node_attributes(NewGraf, {node:xyz for node,xyz in zip(BenzeneAtoms,BenzeneXYZ)}, 'Cartesians')
				Geometry = NewGraf
				Geometry.graph.update({'Etiqueta': InputLabel})
				if RDKit_available and RunRDKIT:
					NewGraf = OptimizeGeometry(NewGraf)
			PrintGeom(Geometry)
			PrintGeomToFile(Geometry,Fitxer)
		quit()
	XYZinput = sys.argv[1]
	print("Input file: "+XYZinput)
	if '.' in XYZinput:
		Titol = XYZinput.split('.')[-2]
	else:
		Titol = XYZinput.strip()
	InitialGeometry = ReadXYZFile(XYZinput)
	print("InitialGeometry=",InitialGeometry)
	StartingGeometry = MakeGrafFromGeometry(InitialGeometry,Titol)
	print("Graph Etiqueta=", StartingGeometry.graph['Etiqueta'])
	print("StartingGeometry=",StartingGeometry)
# NumberOfGenerations generations are done below:
	print("Starting geometry:")
	Fitxer = open(ResultFile,"w")
	PrintGeom(StartingGeometry)
	PrintGeomToFile(StartingGeometry,Fitxer)
	print(NumberOfGenerations,"generations will be done")
	ResultingGeometries = [StartingGeometry]
	for NGen in [str(r+1) for r in range(NumberOfGenerations)]: # NGen is the generation number as a string
		if RestoreRequested and NGen in RestoreGenerations:
			# Read graphs from the binary file
			GrafsReadDict = nx.read_gpickle(Titol+"-g-"+NGen+".pickle")
			ResultingGeometriesNextGen = [GrafsReadDict[r] for r in GrafsReadDict.keys()]
			print("Generation "+NGen+":",len(ResultingGeometriesNextGen),"geometries were read in from file",Titol+"-g-"+NGen+".pickle")
		else:
			print("\nDoing generation "+NGen)
			BenzeneGraf = MakeBenzeneGraf(NGen)
			ResultingGeometriesNextGen = []
			for NewGeo in ResultingGeometries:
				PartialGeoms = DoJob(NewGeo,BenzeneGraf,ResultingGeometriesNextGen,Fitxer)
				ResultingGeometriesNextGen.extend(PartialGeoms)
				print(len(PartialGeoms),NGen+"-generation structures have been generated from",NewGeo.graph['Etiqueta'])
			print(NGen+"-th generation: A total of",len(ResultingGeometriesNextGen),"unique structures were generated")
			# Write the graphs to disk in binary format
			if BackupRequested:
				for Geo in ResultingGeometriesNextGen:
					nx.write_gpickle({r.graph['Etiqueta']:r for r in ResultingGeometriesNextGen}, Titol+"-g-"+NGen+".pickle")
		ResultingGeometries = ResultingGeometriesNextGen
	Fitxer.close()

def DoJob(StartingGeometry,BenzeneGraf,GeometriesAlreadyDone,Fitxer):
	SuitableSides = FindSuitableSides(StartingGeometry)
	Le = len(SuitableSides)
	print(Le,"SuitableSides=",SuitableSides)
	if Randomi[NGen]:
		Nrand = max(round(RandomizedGen[NGen]*len(SuitableSides)),1)
		print("Random",Nrand,"suitable sides out of",len(SuitableSides),"will be used")
		SuitableSides = [SuitableSides[r] for r in sorted(random.sample(range(Le),Nrand))]
		print(Nrand,"SuitableSides to be used:",SuitableSides)
	ResultingGrafs = []
	NewGrafUnlinked = nx.union(StartingGeometry,BenzeneGraf)
	BenzeneAtoms = ['C:1g'+NGen,'C:2g'+NGen,'C:3g'+NGen,'C:4g'+NGen,'H:1g'+NGen,'H:2g'+NGen,'H:3g'+NGen,'H:4g'+NGen]
	Etiqueta = StartingGeometry.graph['Etiqueta']
	for r in SuitableSides:
		NumC1,NumC2,NumH1,NumH2 = r
		NewGraf = NewGrafUnlinked.copy()
		BenzeneXYZ = MakeBenzeneXYZ(StartingGeometry,r)
		NewGraf.remove_nodes_from([NumH1,NumH2])
		NewGraf.add_edges_from([(NumC1,'C:4g'+NGen), (NumC2,'C:1g'+NGen)])
		nx.set_node_attributes(NewGraf, {node:xyz for node,xyz in zip(BenzeneAtoms,BenzeneXYZ)}, 'Cartesians')
		EtiquetaNova = Etiqueta+':'+NumC2.replace(':','',1)+'-'+NumC1.replace(':','',1)
		NewGraf.graph.update({'Etiqueta': EtiquetaNova})
		isom = CheckStericsRequested and Congested(NewGraf) # checking whether the current graph is sterically hindered
		Inchi = GraphToInchiKey(NewGraf)
		NewGraf.graph.update({'Inchi': Inchi})
		isom = isom or (Inchi in set(g.graph['Inchi'] for g in ResultingGrafs + GeometriesAlreadyDone))
		if isom: continue
		if RDKit_available and RunRDKIT:
			NewGraf = OptimizeGeometry(NewGraf)
		if not isom:
			print("Inchi = ",Inchi)
			print("The new graph "+EtiquetaNova+" is NOT found isomorphic to a previous graph, hence accepted")
			ResultingGrafs.append(NewGraf)
			PrintGeom(NewGraf)
			PrintGeomToFile(NewGraf,Fitxer)
	return ResultingGrafs

def GraphToInchiKey(Graf):
	from rdkit import Chem
	from rdkit.Chem import inchi
	Atoms = list(Graf.nodes)
	XYZ = nx.get_node_attributes(Graf, 'Cartesians')
	XYZ = [XYZ[r] for r in Atoms]
	ChemElement = [r.split(":")[0] for r in Atoms]
	AtomDictLocal = {r:i for i,r in enumerate(Atoms)}

	EdgeListRDKIT = []
	BondTypeRDKIT = []

	for p,q in Graf.edges():
		i = AtomDictLocal[p]
		j = AtomDictLocal[q]
		EdgeListRDKIT.append((i,j))
		if ChemElement[i] == 'H' or ChemElement[j] == 'H':
			BondTypeRDKIT.append('single')
		else:
			BondTypeRDKIT.append('aromatic')

	mol = construct_molecule(XYZ, ChemElement, EdgeListRDKIT, BondTypeRDKIT)
	try:
		Chem.SanitizeMol(mol)
		return inchi.MolToInchiKey(mol)
	except:
		return None

def MakeBenzeneGraf(NGen):
	BenzeneAtomsC = ['C:1g'+NGen,'C:2g'+NGen,'C:3g'+NGen,'C:4g'+NGen]
	BenzeneAtomsH = ['H:1g'+NGen,'H:2g'+NGen,'H:3g'+NGen,'H:4g'+NGen]
	Graf = nx.Graph(
    (('C:1g'+NGen,'C:2g'+NGen),('C:2g'+NGen,'C:3g'+NGen),('C:3g'+NGen,'C:4g'+NGen),
    ('C:1g'+NGen,'H:1g'+NGen),('C:2g'+NGen,'H:2g'+NGen),('C:3g'+NGen,'H:3g'+NGen),
    ('C:4g'+NGen,'H:4g'+NGen))
    )
	nx.set_node_attributes(Graf, {r:"C" for r in BenzeneAtomsC}, 'ChemElement')
	nx.set_node_attributes(Graf, {r:"H" for r in BenzeneAtomsH}, 'ChemElement')
	return Graf

def MakeBenzeneXYZ(Graf,r):
	RotMatrix, MidPoint = CalculateRotMatrix(Graf,r)
	BenzeneXYZ = np.array(
	[[-1.3963, 1.2092, 0.0],
	 [-0.6981, 2.4185, 0.0],
	 [ 0.6981, 2.4185, 0.0],
	 [ 1.3963, 1.2092, 0.0],
	 [-2.4826, 1.2092, 0.0],
	 [-1.2413, 3.3593, 0.0],
	 [ 1.2413, 3.3593, 0.0],
	 [ 2.4826, 1.2092, 0.0]]
	)
	BenzeneXYZ = np.matmul(BenzeneXYZ,RotMatrix)
	for i in range(len(BenzeneXYZ)):
		BenzeneXYZ[i] += MidPoint
	return BenzeneXYZ

def CalculateRotMatrix(Graf,r):
	NumC1,NumC2,NumH1,NumH2 = r
	XYZ = nx.get_node_attributes(Graf, 'Cartesians')
	VecA = XYZ[NumC1]-XYZ[NumC2]
	VecA /= np.linalg.norm(VecA)
	VecH = (XYZ[NumH1]-XYZ[NumC1]) + (XYZ[NumH2]-XYZ[NumC2])
	VecH /= np.linalg.norm(VecH)
	VecB = VecH - np.dot(VecA,VecH)*VecA
	VecB /= np.linalg.norm(VecB)
	if np.dot(VecB,VecH)<0:
		VecB = -VecB
	VecC = np.cross(VecA,VecB)
	MidPoint = (XYZ[NumC1]+XYZ[NumC2])/2
	return np.array([VecA,VecB,VecC]), MidPoint

def FindHydrogens(Geometry,carpa):
	p,q = carpa
	print("edges:",Geometry.edges())
	if Geometry.nodes[p]['ChemElement'] == 'H':
		print("Error in label: node "+p+" is hydrogen")
		quit()
	if Geometry.nodes[q]['ChemElement'] == 'H':
		print("Error in label: node "+p+" is hydrogen")
		quit()
	if not Geometry.has_edge(p,q):
		print("Error in label: node "+p+" is not linked to node "+q)
		quit()
	ptrue = False
	for r in Geometry.neighbors(p):
		ptrue = Geometry.nodes[r]['ChemElement'] == 'H'
		if ptrue: break
	qtrue = False
	for s in Geometry.neighbors(q):
		qtrue = Geometry.nodes[s]['ChemElement'] == 'H'
		if qtrue: break
	if not ptrue:
		print("Error in label: node "+p+" is not linked to a hydrogen")
		quit()
	elif not qtrue:
		print("Error in label: node "+q+" is not linked to a hydrogen")
		quit()
	else:
		print("(q,p,s,r )=",q,p,s,r)
		return q,p,s,r # C,C,H,H

def FindSuitableSides(Geometry):
	#for p,q in Geometry.__dict__.items(): print(p,q)
	HeavyListWithHydrogen = []
	SuitableSides = []
	for p,q in Geometry.edges():
		if Geometry.nodes[p]['ChemElement'] != 'H' and Geometry.nodes[q]['ChemElement'] == 'H':
			HeavyListWithHydrogen.append((p,q)) # C-H bond
		elif Geometry.nodes[p]['ChemElement'] == 'H' and Geometry.nodes[q]['ChemElement'] != 'H':
			HeavyListWithHydrogen.append((q,p)) # C-H bond

	for i,(p,q) in enumerate(HeavyListWithHydrogen):
		for k,(r,s) in enumerate(HeavyListWithHydrogen[:i]):
			if Geometry.has_edge(p,r):
				SuitableSides.append((p,r,q,s)) # C,C,H,H
	return SuitableSides

def MakeGrafFromGeometry(Geometry,Titol):
	global RunRDKIT
	AtomsHeavy,AtomsHyd,XYZHeavy,XYZHyd = Geometry
	N = len(AtomsHeavy) + len(AtomsHyd)
	XYZ = XYZHeavy + XYZHyd
	Atoms = AtomsHeavy + AtomsHyd
	print("Atoms =",Atoms)
	EdgeList = []
	#ChemElement = [AtomicNamesDict[r.split(":")[0]] for r in Atoms]
	for i in range(N):
		for j in range(i):
			d = np.linalg.norm(XYZ[i]-XYZ[j])
			#d = np.sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2+(Z[i]-Z[j])**2)
			if(d>0.85 and d<CovalentRadius[AtomicNamesDict[Atoms[i].split(":")[0]]]+CovalentRadius[AtomicNamesDict[Atoms[j].split(":")[0]]]+0.3):
				EdgeList.append((Atoms[j],Atoms[i]))
	print("Initial edge list:",EdgeList)
	Graf = nx.Graph()
	Graf.graph.update({'Etiqueta': Titol})
	Graf.add_nodes_from(Atoms)
	Graf.add_edges_from(EdgeList)
	nx.set_node_attributes(Graf, {r:r.split(":")[0] for r in Atoms}, 'ChemElement')
	nx.set_node_attributes(Graf, {r:xyz for r,xyz in zip(Atoms,XYZ)}, 'Cartesians')

	if RDKit_available and RunRDKIT:
		Graf = OptimizeGeometry(Graf)
	return Graf

def OptimizeGeometry(Graf): # Doing an MM optimization through RDkit
	global OptimizerRDKIT, ForceField
	Atoms = nx.nodes(Graf)
	print("OptimizeGeometry: Atoms=",Atoms)
	Etiqueta = Graf.graph['Etiqueta']
	AtomDictLocal = {r:i for i,r in enumerate(Atoms)}
	print("AtomDictLocal =",AtomDictLocal)
	XYZ = nx.get_node_attributes(Graf, 'Cartesians')
	XYZ = [XYZ[r] for r in Atoms]
	#ChemElement = list(nx.get_node_attributes(Graf, 'ChemElement').values())
	ChemElement = [r.split(":")[0] for r in Atoms]
	EdgeListRDKIT = []
	BondTypeRDKIT = []

	for p,q in Graf.edges():
		i = AtomDictLocal[p]
		j = AtomDictLocal[q]
		EdgeListRDKIT.append((i,j))
		if ChemElement[i] == 'H' or ChemElement[j] == 'H':
			BondTypeRDKIT.append('single')
		else:
			BondTypeRDKIT.append('aromatic')
	mol = construct_molecule(XYZ, ChemElement, EdgeListRDKIT, BondTypeRDKIT) # RDkit molecule construction
	try:
		Chem.SanitizeMol(mol)
		for bnd in mol.GetBonds():
			print("Bond type:",bnd.GetBeginAtomIdx()+1,"-",bnd.GetEndAtomIdx()+1,bnd.GetBondType())
		#AllChem.UFFOptimizeMolecule(mol)
		#AllChem.MMFFOptimizeMolecule(mol)
		OptimizerRDKIT(mol)
		print("Molecule "+Etiqueta+" was "+ForceField+"optimized")
		XYZ = np.array([mol.GetConformer().GetAtomPosition(i) for i in range(len(XYZ))])
		nx.set_node_attributes(Graf, {r:xyz for r,xyz in zip(Atoms,XYZ)}, 'Cartesians')
	except:
		print("Geometry optimization for "+Etiqueta+" failed. Non-optimized geometry retained.")
	return Graf


def construct_molecule(XYZ, ElementSymbols, Bindungen, BondTypes): # This is for RDkit
	mol = Chem.RWMol()  # Create a writable molecule object
	conf = Chem.Conformer(len(ElementSymbols))
	for r,xyz in zip(ElementSymbols,XYZ):
		atom = Chem.Atom(r)
		idx = mol.AddAtom(atom)
		conf.SetAtomPosition(idx, xyz)
	for bnd, bndtype in zip(Bindungen, BondTypes):
		if bndtype == 'single':
			mol.AddBond(*bnd, Chem.BondType.SINGLE)
		elif bndtype == 'double':
			mol.AddBond(*bnd, Chem.BondType.DOUBLE)
		elif bndtype == 'aromatic':
			mol.AddBond(*bnd, Chem.BondType.AROMATIC)
	mol.AddConformer(conf)
	return mol.GetMol()

def ReadXYZFile(XYZFileName):
	XYZFile = open(XYZFileName,"r")
	a = [r for r in XYZFile.readlines() if len(r.split())>=4]
	XYZFile.close()
	AtomsHeavy = []
	AtomsHyd = []
	XYZHeavy = []
	XYZHyd = []
	NHeavy = 0
	NHyd = 0
	for r in a:
		r_split = r.split()
		Ordnungszahl = AtomicNamesDict[r_split[0]]
		if Ordnungszahl>1:
			NHeavy += 1
			AtomsHeavy.append(AtomicNames[Ordnungszahl]+":"+str(NHeavy))
			XYZHeavy.append(np.array([float(r_split[1]),float(r_split[2]),float(r_split[3])]))
		else:
			NHyd += 1
			AtomsHyd.append("H:"+str(NHyd))
			XYZHyd.append(np.array([float(r_split[1]),float(r_split[2]),float(r_split[3])]))
	del a
	return AtomsHeavy, AtomsHyd, XYZHeavy, XYZHyd

def PrintGeom(Graf):
	#for p,q in Graf.__dict__.items(): print(p,q)
	global AdjMatrixRequested
	N = nx.number_of_nodes(Graf)
	ChemElements = nx.get_node_attributes(Graf, 'ChemElement')
	XYZ = nx.get_node_attributes(Graf, 'Cartesians')
	print
	print("\nName: ",Graf.graph['Etiqueta'])
	print("-----------------------------------------")
	print("%-3iAtoms             Coordinates            "%N)
	print("-----------------------------------------")
	for p in nx.nodes(Graf):
		print("%s%11.4f%14.4f%14.4f"%(ChemElements[p],*XYZ[p]))
	print("-----------------------------------------")
	if AdjMatrixRequested:
		print("Adjacency matrix:")
		for r in nx.to_numpy_array(Graf).astype(int):
			print(" ".join(["%s"%p for p in r]))
	if AdjMatrixRequested:
		print("Bonds:")
		adj_matrix = nx.to_numpy_array(Graf).astype(int)
		for i in range(N):
			for j in range(i):
				if adj_matrix[i,j]==1: print(j+1,i+1)

def PrintGeomToFile(Graf,Fitxer):
	N = nx.number_of_nodes(Graf)
	ChemElements = nx.get_node_attributes(Graf, 'ChemElement')
	XYZ = nx.get_node_attributes(Graf, 'Cartesians')
	Fitxer.write(str(N)+"\n")
	for p in Graf.graph['Etiqueta']:
		ToPrint = p[0]
		if len(p)>1:
			ToPrint = p[0]+":"+":".join([r for r in p[1:]])+" / "
		Fitxer.write(ToPrint)
	Fitxer.write("\n")
	for p in nx.nodes(Graf):
		Fitxer.write("%s%11.4f%14.4f%14.4f\n"%(ChemElements[p],*XYZ[p]))

def InsertColon(s):
	for i,char in enumerate(s):
		if char.isdigit(): break
	return s[:i]+":"+s[i:]

def Congested(Graf):
	Nodes = nx.nodes(Graf)
	XYZ = nx.get_node_attributes(Graf,'Cartesians')
	BadPairs =[]
	for i,p in enumerate(Nodes):
		for j,q in enumerate(Nodes):
			if j<i and np.linalg.norm(XYZ[p]-XYZ[q]) < 0.6:
				BadPairs.append((q,p))
		if len(BadPairs)>1: break
	if len(BadPairs)>1:
		print("More than one congested atom pairs in",Graf.graph['Etiqueta'])
		return True
	elif len(BadPairs)==1:
		print("One congested atom pair in",Graf.graph['Etiqueta'])
		return True
	else:
		return False

MainProgram()
