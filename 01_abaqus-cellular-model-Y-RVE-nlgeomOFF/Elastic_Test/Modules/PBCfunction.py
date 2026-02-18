# PBCfunction.py
# Author Alessia Ferrara, March 2021
# Last modified November 2023
#
# This script contains the following functions:
# - PeriodicBoundCond -> Apply periodic boundary conditions to periodic edges/faces of an RVE. The RVE can be in 2D or extruded in 3D.
# - planePBC -> Apply PBC in 2D model.
# - fullPBC -> -> Apply PBC in 3D model.
# - checkPairToConstrain -> Check if the number of couple of edges is equal or not to the number of lattice vectors
# - checkTolerance -> Check if the passed tolerance is less or not than the minimum distance in x-and y-direction between the nodes of the same set. Used in PeriodicBound before the pairing routine.

from part import *
from material import *
from section import *
from assembly import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from sketch import *
from visualization import *
from connectorBehavior import *
from math import *
import numpy as np
from copy import *

def PeriodicBoundCond(mdb, NameModel, periodicEdgesName, LatticeVector, strains, tol, dim):
# This function apply periodic boundary conditions on corresponding edges (2D) or faces (3D) of an RVE.
# First, it creates Reference Points to which apply the macro-strains or stresses, then it pairs the
# nodes of periodic edges/faces, and finally applies constraint equations between them.
# The function also requires the values of the macro-strain that will be applied since the equations
# change according to whether the strain is free or not.
# Args:
# - mdb = database of the model
# - NameModel = name of the model
# - periodicEdgesName = array of names of nodes sets to pair
# - LatticeVector = array of lattice vectors 
# - strains = array of strains
# - tol = tolerance to find paired nodes
# - dim = model dimension (2 for 2D or 3 for 3D)
# Note: The order of the sets names should be consistent w.r.t. the order of the lattice vectors.
# Return:
# - RefPointName = array of RP names
# - mdb = database of the model
	
	from part import DEFORMABLE_BODY
	myModel = mdb.models[NameModel]
	myAssembly = myModel.rootAssembly

	# Check modeling space
	if dim==2:
		from part import TWO_D_PLANAR
		dimensionality = TWO_D_PLANAR
	elif dim==3:
		from part import THREE_D
		dimensionality = THREE_D
	else:
		raise ValueError("Parameter dim must be 2 for 2D geometry or 3 for 3D geometry. You passed {}.".format(dim))
	# Check number of edges to pair Vs number of lattice vectors
	numCouples = checkPairToConstrain(periodicEdgesName, LatticeVector)

	# Create Reference Points (RP): one RP for each strain component and for each couple
	if len(strains)==3:
		StrainComp = ['X', 'Y', 'XY-X', 'XY-Y']
	elif len(strains)==6:
		StrainComp = ['X', 'Y', 'Z', 'YZ-Y', 'YZ-Z', 'XZ-X', 'XZ-Z', 'XY-X', 'XY-Y']
	else:
		raise ValueError("The vector of strains must contain 3 values for plane PBCs or 6 for full PBCs. You passed a vector of {} elements.".format(len(strains)))
	# Initialize vectors of RP names
	numStrainComp = len(StrainComp)
	RefPointName = [['']*numCouples for i in range(numStrainComp)]
	# Create RPs
	for i in range(numStrainComp): # component
		for k in range(numCouples): # couple
			RefPointName[i][k] = ('RefPoint-'+StrainComp[i]) # name of RP
		# Create part
		myModel.Part(dimensionality=dimensionality, name=RefPointName[i][k], type=DEFORMABLE_BODY)
		myModel.parts[RefPointName[i][k]].ReferencePoint(point=(0.0, 0.0, 0.0))
		# Create istance
		myAssembly.Instance(dependent=ON, name=RefPointName[i][k], part=myModel.parts[RefPointName[i][k]])
		# Create set
		myAssembly.Set(name=RefPointName[i][k], referencePoints=(myAssembly.instances[RefPointName[i][k]].referencePoints[1],))

	# Find paired nodes and apply constraint equations
	for couple in range(numCouples): # repeat for each couple of edges
		numEq = 0 # counter for number of equations
		nodes0 = myAssembly.sets[periodicEdgesName[2*couple]].nodes # nodes of edge 0
		nodes1 = myAssembly.sets[periodicEdgesName[2*couple+1]].nodes # nodes of edge 1
		# Check tolerance to find paired nodes Vs min distance between nodes
		#tol = checkTolerance(NameModel, couple, nodes1, tol) # only nodes1 is checked since the pairing routine is repeated for each node of nodes0 and the corresponding node is searched in nodes1
		checkNodes1=range(0,len(nodes1)) # sequence of number from 0 to len(nodes1)-1
		for nod0 in range(len(nodes0)):	# repeat for each node of nodes0
			
			stop=False # stop become true when paired nodes are found and the loop stops
			for nod1 in checkNodes1: # repeat for each node of nodes1
				if len(strains)==3:
					numEq, checkNodes1, stop = planePBC(mdb, NameModel, LatticeVector, couple, RefPointName, strains, numEq, checkNodes1, nodes0, nodes1, nod0, nod1, tol)
				else:
					numEq, checkNodes1, stop = fullPBC(mdb, NameModel, LatticeVector, couple, RefPointName, strains, numEq, checkNodes1, nodes0, nodes1, nod0, nod1, tol)
				
				if stop:
					break # exit from loop if paired nodes have been found

		# Check if constraint equations have been set, and if so check if among the constrained nodes there are out-of-edge nodes, i.e. internal nodes
		if numEq==0:
			print('WARNING: No constraint equation has been set for couple {}.'.format(couple+1))
			print('Check the corresponding sets of nodes {} and {}, and the lattice vector {}.'.format(periodicEdgesName[2*couple],periodicEdgesName[2*couple+1],couple+1))
			print('Check also the order in which the sets are passed to the function: the lattice vector goes from the first to the second nodes set passed to the function through the vector of sets name.')

	return (RefPointName, mdb)
# End: def PeriodicBoundCond

def planePBC(mdb, NameModel, LatticeVector, couple, RefPointName, strains, numEq, checkNodes1, nodes0, nodes1, nod0, nod1, tol):
# This function apply periodic boundary conditions for a 2D model.
# Args:
# - mdb = database of the model
# - NameModel = name of the model
# - LatticeVector = array of lattice vectors 
# - couple = index of the lattice vector
# - RefPointName = array of RP names
# - strains = array of strains
# - numEq = number of equations already created
# - checkNodes1 = array of index of nodes on edge 1
# - nodes0 = set of nodes of edge 0
# - nodes1 = set of nodes of edge 1
# - nod0 = index of node on edge 0
# - nod1 = index of node on edge 1
# - tol = tolerance to find paired nodes
# Note: The lattice vector goes from edge 0 to edge 1.
# Return:
# - numEq = new number of equations created
# - checkNodes1 = array of index of remaining nodes on edge 1
# - stop = boolean variable, True if a match has been found, False if not

	myModel = mdb.models[NameModel]
	myAssembly = myModel.rootAssembly
	LatticeVec = deepcopy(LatticeVector)
	stop=False # stop become true when paired nodes are found and the loop stops
	Coord0 = nodes0[nod0].coordinates # coord. of node on edge 0
	Coord1 = nodes1[nod1].coordinates # coord. of node on edge 1
	dx=Coord1[0]-Coord0[0]	# x-distance between node 1 and node 0
	dy=Coord1[1]-Coord0[1]	# y-distance between node 1 and node 0

	# Paired nodes found when the distance between node 1 and node 0 is equal to the lattice vector within the tolerance
	if abs(LatticeVec[couple][0]-dx)<=tol and abs(LatticeVec[couple][1]-dy)<=tol:
		# Create a set for each node to pair
		# Node name = Node - 'num. of the Couple (1, 2, 3,..)' - 'num. of the Node (num. of Eq.)' - 'num. of the Edge (0 or 1)' 
		myAssembly.Set(name='Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), nodes=nodes0[nod0:nod0+1])
		myAssembly.Set(name='Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1), nodes=nodes1[nod1:nod1+1])
		# Create constraint equations
		# Equation name = Eq - 'num. of the Couple (1, 2, 3,..)' - 'num. of the Node (=num. of Eq.)' - 'DOF'
		c1 = 1.0
		c2 = -1.0
		# DOF 1 (x-constraint)
		if strains[0]==UNSET and strains[2]==UNSET: # CASE UNSET, UNSET (free x-displacement)
			if strains[0]==UNSET and strains[1]==UNSET: # pure shear
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
			else: # pure axial
				LatticeVec[couple][1] = 0
		else:
			if strains[0]==UNSET:
				LatticeVec[couple][0] = 0
			if strains[2]==UNSET:
				LatticeVec[couple][1] = 0
		# u_xi - u_xj + dx * eps_x + dy * eps_xy + dz * eps_xz = 0
		myModel.Equation(name='Eq-'+str(couple+1)+'-'+str(numEq)+'-X',
			terms=((c1,'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), 1),(c2, 'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1),1),
			(LatticeVec[couple][0], RefPointName[0][couple], 1), (LatticeVec[couple][1], RefPointName[2][couple], 1)))
		
		# DOF 2 (y-constraint)
		LatticeVec = deepcopy(LatticeVector)
		if strains[1]==UNSET and strains[2]==UNSET: # Case UNSET, UNSET (free y-displacement)
			if strains[0]==UNSET and strains[1]==UNSET: # pure shear
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
			else: # pure axial
				LatticeVec[couple][0] = 0
		else:
			if strains[2]==UNSET:
				LatticeVec[couple][0] = 0
			if strains[1]==UNSET:
				LatticeVec[couple][1] = 0
		# u_yi - u_yj + dx * eps_xy + dy * eps_y + dz * eps_yz = 0
		myModel.Equation(name='Eq-'+str(couple+1)+'-'+str(numEq)+'-Y',
			terms=((c1,'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), 2),(c2, 'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1),2),
			(LatticeVec[couple][0], RefPointName[3][couple], 1), (LatticeVec[couple][1], RefPointName[1][couple], 1)))
		numEq = numEq+1 # increase number of equations
		checkNodes1.remove(nod1) # remove paired node from the sequence 
		stop=True # since paired nodes have been found, don't look further and go to following node

	return numEq, checkNodes1, stop
# End: planePBC

def fullPBC(mdb, NameModel, LatticeVector, couple, RefPointName, strains, numEq, checkNodes1, nodes0, nodes1, nod0, nod1, tol):
# This function apply periodic boundary conditions for a 3D model.
# Args:
# - mdb = database of the model
# - NameModel = name of the model
# - LatticeVector = array of lattice vectors 
# - couple = index of the lattice vector
# - RefPointName = array of RP names
# - strains = array of strains
# - numEq = number of equations already created
# - checkNodes1 = array of index of nodes on face 1
# - nodes0 = set of nodes of face 0
# - nodes1 = set of nodes of face 1
# - nod0 = index of node on face 0
# - nod1 = index of node on face 1
# - tol = tolerance to find paired nodes
# Note: The lattice vector goes from face 0 to face 1.
# Return:
# - numEq = new number of equations created
# - checkNodes1 = array of index of remaining nodes on face 1
# - stop = boolean variable, True if a match has been found, False if not

	myModel = mdb.models[NameModel]
	myAssembly = myModel.rootAssembly
	LatticeVec = deepcopy(LatticeVector)
	stop=False # stop become true when paired nodes are found and the loop stops
	Coord0 = nodes0[nod0].coordinates # coord. of node on face 0
	Coord1 = nodes1[nod1].coordinates # coord. of node on face 1
	dx=Coord1[0]-Coord0[0]	# x-distance between node 1 and node 0
	dy=Coord1[1]-Coord0[1]	# y-distance between node 1 and node 0
	dz=Coord1[2]-Coord0[2]	# z-distance between node 1 and node 0

	# Paired nodes found when the distance between node 1 and node 0 is equal to the lattice vector within the tolerance
	if abs(LatticeVec[couple][0]-dx)<=tol and abs(LatticeVec[couple][1]-dy)<=tol and abs(LatticeVec[couple][2]-dz)<=tol:
		# Create a set for each node to pair
		# Node name = Node - 'num. of the Couple (1, 2, 3,..)' - 'num. of the Node (num. of Eq.)' - 'num. of the Edge (0 or 1)' 
		myAssembly.Set(name='Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), nodes=nodes0[nod0:nod0+1])
		myAssembly.Set(name='Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1), nodes=nodes1[nod1:nod1+1])
		# Create constraint equations
		# Equation name = Eq - 'num. of the Couple (1, 2, 3,..)' - 'num. of the Node (=num. of Eq.)' - 'DOF'
		c1 = 1.0
		c2 = -1.0
		# DOF 1 (x-constraint)
		if strains[0]==UNSET and strains[5]==UNSET and strains[4]==UNSET: # CASE UNSET, UNSET, UNSET (free x-displacement)
			if strains[0]==UNSET and strains[1]==UNSET and strains[2]==UNSET: # pure shear
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
				LatticeVec[couple][2] = 0
			else: # pure axial
				LatticeVec[couple][1] = 0
				LatticeVec[couple][2] = 0
		else:
			if strains[0]==UNSET:
				LatticeVec[couple][0] = 0
			if strains[5]==UNSET:
				LatticeVec[couple][1] = 0
			if strains[4]==UNSET:
				LatticeVec[couple][2] = 0
		# u_xi - u_xj + dx * eps_x + dy * eps_xy + dz * eps_xz = 0
		myModel.Equation(name='Eq-'+str(couple+1)+'-'+str(numEq)+'-X',
			terms=((c1,'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), 1),(c2, 'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1),1),
			(LatticeVec[couple][0], RefPointName[0][couple], 1), (LatticeVec[couple][1], RefPointName[7][couple], 1),
			(LatticeVec[couple][2], RefPointName[5][couple], 1)))
		
		# DOF 2 (y-constraint)
		LatticeVec = deepcopy(LatticeVector)
		if strains[1]==UNSET and strains[3]==UNSET and strains[5]==UNSET: # Case UNSET, UNSET, UNSET (free y-displacement)
			if strains[0]==UNSET and strains[1]==UNSET and strains[2]==UNSET: # pure shear
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
				LatticeVec[couple][2] = 0
			else: # pure axial
				LatticeVec[couple][0] = 0
				LatticeVec[couple][2] = 0
		else:
			if strains[5]==UNSET:
				LatticeVec[couple][0] = 0
			if strains[1]==UNSET:
				LatticeVec[couple][1] = 0
			if strains[3]==UNSET:
				LatticeVec[couple][2] = 0
		# u_yi - u_yj + dx * eps_xy + dy * eps_y + dz * eps_yz = 0
		myModel.Equation(name='Eq-'+str(couple+1)+'-'+str(numEq)+'-Y',
			terms=((c1,'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), 2),(c2, 'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1),2),
			(LatticeVec[couple][0], RefPointName[8][couple], 1), (LatticeVec[couple][1], RefPointName[1][couple], 1),
			(LatticeVec[couple][2], RefPointName[3][couple], 1)))
			
		# DOF 3 (z-constraint)
		LatticeVec = deepcopy(LatticeVector)
		if strains[4]==UNSET and strains[3]==UNSET and strains[2]==UNSET: # Case UNSET, UNSET, UNSET (free z-displacement)
			if strains[0]==UNSET and strains[1]==UNSET and strains[2]==UNSET: # pure shear
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
				LatticeVec[couple][2] = 0
			else: # pure axial
				LatticeVec[couple][0] = 0
				LatticeVec[couple][1] = 0
		else:
			if strains[4]==UNSET:
				LatticeVec[couple][0] = 0
			if strains[3]==UNSET:
				LatticeVec[couple][1] = 0
			if strains[2]==UNSET:
				LatticeVec[couple][2] = 0
		# u_zi - u_zj + dx * eps_xz + dy * eps_yz + dz * eps_z = 0
		myModel.Equation(name='Eq-'+str(couple+1)+'-'+str(numEq)+'-Z',
			terms=((c1,'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(0), 3),(c2, 'Node-'+str(couple+1)+'-'+str(numEq)+'-'+str(1),3),
			(LatticeVec[couple][0], RefPointName[6][couple], 1), (LatticeVec[couple][1], RefPointName[4][couple], 1),
			(LatticeVec[couple][2], RefPointName[2][couple], 1)))

		numEq = numEq+1 # increase number of equations
		checkNodes1.remove(nod1) # remove paired node from the sequence 
		stop=True # since paired nodes have been found, don't look further and go to following node

	return numEq, checkNodes1, stop
# End: planePBC

def checkPairToConstrain(periodicEdgesName, LatticeVec):
# This function checks if the number of couples of edges is equal to the number of lattice vectors
# Args: 
# - periodicEdgesName = name of the edges
# - LatticeVec = array of lattice vectors
# Return:
# - numCouples = number of couple of edges to constrain

	numCouples = len(periodicEdgesName)/2 # number of couples of edges
	numLVs = len(LatticeVec)

	if numCouples > numLVs: # if there are more couple than lattice vectors, constrain only the couples corresponding to the given lattice vectors
		print('WARNING: The number of the couples of edges to pair is greater than the number of the given lattice vectors. Only the couples corresponding to the given lattice vectors will be constrained.')
		print('The unused set of nodes are:')
		for i in range(numLVs, numCouples):
			print(periodicEdgesName[2*i])
			print(periodicEdgesName[2*i+1])
		numCouples = numLVs # use only the couples corresponding to the given lattice vectors
	elif numCouples < numLVs: # if there are more lattice vectors than couple, use only the lattice vectors corresponding to the given couples
		print('WARNING: The number of the given lattice vectors is greater than the number of the couples of edges to pair. Only the lattice vectors corresponding to the given couples of edges will be used to constrain them.')	
		print('The unused lattice vector are:')
		for i in range(numCouples, numLVs):
			print("LV{}".format(i+1))

	return numCouples
# End: def checkPairToConstrain

def checkTolerance(NameModel, couple, nodes, tol):
# This function checks if the tolerance for the pairing routine is less or not than the minimum distance between the nodes
# Args: 
# - NameModel = name of the model
# - couple = number of the couple
# - nodes = array of nodes
# - tol = tolerance to find paired nodes
# Return:
# - min_tol = new values of tolerance to find paired nodes

	myModel = mdb.models[NameModel]
	myAssembly = myModel.rootAssembly
	# Initialize array of coordinates	
	Coord = np.zeros((len(nodes),3))
	# Store coordinates
	for nod in range(len(nodes)):	
			Coord[nod] = nodes[nod].coordinates

	# Check tolerance Vs min distance between nodes
	min_dist = tol*2
	min_tol = tol*2
	red_tol = False
	WarnNodeVec = []
	for i in range(len(nodes)):
		for j in range(i+1,len(nodes)):
			x_dist = abs(Coord[i][0]-Coord[j][0])
			y_dist = abs(Coord[i][1]-Coord[j][1])
			z_dist = abs(Coord[i][2]-Coord[j][2])			
			if 0<x_dist<min_dist and 0<y_dist<min_dist and 0<z_dist<min_dist:
				WarnNodeVec = WarnNodeVec+[nodes[i]]
				WarnNodeVec = WarnNodeVec+[nodes[j]]				
				red_tol = True
				if 0<min(x_dist, y_dist, z_dist)<min_tol:
					min_tol = min(x_dist, y_dist, z_dist)

	min_tol = min_tol/2. # new min tolerance

	if red_tol==True:
		print('Couple {}, edge 1: '.format(couple+1))
		print('there are some nodes distances that are less than the set tolerance, the tolerance will be set to {} to avoid mispicking nodes during the pairing process.'.format(min_tol))
		print('The concerned nodes can be found in the assembly set WarnNodeTolerance-{}.'.format(couple+1))
		myAssembly.Set(nodes=MeshNodeArray(WarnNodeVec), name='WarnNodeTol-{}'.format(couple+1))

	return min_tol
# End: def checkTolerance

