# Imports
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import regionToolset as regionToolset
from math import *
from copy import copy, deepcopy
import numpy as np
import glob

# Import PBC library from current directory
import imp
import os
cwd = os.getcwd()
import sys
sys.path.insert(0, cwd)
import PBCfunction as PBC
imp.reload(PBC)
from PBCfunction import *

class GR_RVE_Class:

    def __init__(self):
    # Constructor of the Test_Micromaterial_Class class.
    # Args: none

        # Model database
        self.modelName = 'Model-1'
        self.partName = 'GR'
        self.layersName = ['EW', 'TW', 'LW']

        # Geometry variables
        self.t = []
        self.t_tot = 0.
        self.s = 1.
        self.faceName = []
        self.LatticeVector = []
        self.LatticeVector_edge = []
        self.LatticeVector_vert = []
        self.RefPointName = []
        self.numCouples = 0

        # Material variables
        self.depvar = 126
        self.moistureRef = 0.12
        self.StiffMat = np.zeros((6,6))
        
        # Mesh variables
        self.meshSize = 0.
        self.minSizeFactor = 0.
        self.minMeshSize = 0.
        self.nlgeom = OFF
        
        return
    # end: def _init_

    def setup_growth_ring_part(self, t, s):
    # Create a square section of laminate composite material.
    # Args:
    # - t = array of ply thicknesses
    # - s = square section size
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb()
        self.t_original = deepcopy(t)     # e.g. [EW, TW, LW]
        self.layersName_original = self.layersName
        self.s = s

        # Split EW
        half_EW = self.t_original[0] / 2.
        # new thickness list: [EW/2, EW/2, TW, LW]
        self.t = [half_EW] + self.t_original[1:] + [half_EW]
        # new names, just so each half has its own part name
        self.layersName = [self.layersName_original[0] + '_0'] + self.layersName_original[1:] + [self.layersName_original[0] + '_1']
        
        # Recompute layers thickness
        self.t_tot = sum(self.t)
        myModel = self.mdb.models[self.modelName]
        myAssembly = myModel.rootAssembly

        # Build 4 layers
        y = [sum(self.t[:i]) for i in range(len(self.t))]    # [0, half_EW, half_EW+half_EW, half_EW+half_EW+TW]
        count = 0
        for i, name in enumerate(self.layersName):
            # sketch a rectangle from x=y[i] to x=y[i]+t[i], full height
            sk = myModel.ConstrainedSketch(name='__profile__', sheetSize=20.0)
            sk.rectangle(point1=(y[i], 0.), point2=(y[i] + self.t[i], self.t_tot))
            # make the part
            myModel.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
            part = myModel.parts[name]
            part.BaseSolidExtrude(sketch=sk, depth=self.s)
            del myModel.sketches['__profile__']

            # Create set
            myPart = myModel.parts[name]
            part.Set(cells=part.cells, name='Set-' + name)
            # Create instance
            myAssembly.Instance(name=name, part=myPart, dependent=ON)
            
            # merge into a single body as you go
            if i > 0:
                count += 1
                myInstance_temp = myInstance
                myInstance = myAssembly.instances[name]
                merged_name = 'Layer-{}'.format(count)
                myAssembly.InstanceFromBooleanMerge(name=merged_name, instances=(myInstance, myInstance_temp),
                    keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
                # reset for next loop
                myInstance = myAssembly.instances[merged_name + '-1']
            else:
                myInstance = myAssembly.instances[name]

            session.viewports['Viewport: 1'].setValues(displayedObject=myPart)


        # Delete temporary parts
        for i in range(count-1):
            del myModel.parts['Layer-{}'.format(i+1)]
        del myAssembly.instances['Layer-{}-1'.format(count)]
        myModel.parts.changeKey(fromName='Layer-{}'.format(count), toName=self.partName)

        # Create assembly
        myPart = myModel.parts[self.partName]
        myAssembly.Instance(name=self.partName, part=myPart, dependent=ON)

        # Set lattice vectors
        temp = self.t_tot
        self.t_tot = s
        self.s = temp
        self.LatticeVector = [[self.s, 0., 0.], [0., self.s, 0.], [0., 0., self.t_tot]] # periodic faces lattice vector
        self.LatticeVector_edge = [[self.s, 0., self.t_tot], [-self.s, 0., self.t_tot], # periodic edges lattice vector
                                   [self.s, -self.s, 0.], [self.s, self.s, 0.], 
                                   [0., self.s, self.t_tot], [0., self.s, -self.t_tot]]
        self.LatticeVector_vert = [[self.s, self.s, self.t_tot], [self.s, self.s, -self.t_tot], # periodic vertices lattice vector
                                   [self.s, -self.s, self.t_tot], [-self.s, self.s, self.t_tot]]
   
        return self.modelName
    # End: def create_part

    def setup_tissue_material(self):
    # Create and assign materials to composite plies.
    # Args:
    # - materialStiffMat = array of ply stiffness matrices
    # - alpha = array of ply orientations

        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[self.partName]
        
        # Setup materials
        for i in range(len(self.layersName)):
            materialName = self.layersName[i][:2]
            # Create material
            myModel.Material(name=materialName)
            myMaterial = myModel.materials[materialName]
            myMaterial.Depvar(n=self.depvar)
            myMaterial.UserMaterial(mechanicalConstants=(0.0, ))
            # Assign material orientation
            setName = 'Set-'+self.layersName[i]
            # Create and assign section
            sectionName = materialName
            myModel.HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)
            myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                region=myPart.sets[setName], sectionName=sectionName, thicknessAssignment=FROM_SECTION)
        
        myAssembly = myModel.rootAssembly        
        myAssembly.regenerate()

        # Create part set
        setName = np.char.add('Set-', self.layersName)
        myAssembly.SetByBoolean(name='Set-' + self.partName, sets=tuple(myAssembly.allInstances[self.partName].sets[set] for set in setName))

        return
    # End: def create_material

    def setup_mesh(self, meshSize, minSizeFactor):
    # Set up the mesh of the part.
    # Args:
    # - meshSize = global size of the mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the global element size

        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[self.partName]
        myAssembly = myModel.rootAssembly

        self.minMeshSize = meshSize * minSizeFactor
        self.meshSize = meshSize

        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=minSizeFactor, size=meshSize)

        # Set mesh control and element type
        myPart.setMeshControls(elemShape=HEX_DOMINATED, regions=myPart.cells, technique=SWEEP) # set mesh control
        myPart.setElementType(elemTypes=(ElemType(elemCode=C3D20R, elemLibrary=STANDARD),
            ElemType(elemCode=C3D15, elemLibrary=STANDARD),
            ElemType(elemCode=C3D10, elemLibrary=STANDARD)), regions=(myPart.cells, )) # set element type

        # Generate mesh
        myPart.generateMesh()
        myAssembly.regenerate()

        return
    # end: def setup_mesh

    def create_step(self, modelName='', stepName='', nlgeom=ON, previous='Initial', maxNumInc=100, initialInc=1e-5, minInc=1e-5, maxInc=1):
    # Create a static step in the given model.
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - nlgeom = ON/OFF to set non-linear geometry
    # - previous = name of the previous step, most likely "Initial"
    # - maxNumInc = maximum number of increments
    # - initialInc = initial increment size
    # - minInc = minimum increment size
    # - maxInc = maximum increment size

        self.nlgeom = copy(nlgeom)
        if modelName=='':
            modelName = self.modelName
        myModel = self.mdb.models[modelName]

        myModel.StaticStep(name=stepName, nlgeom=self.nlgeom, previous=previous, initialInc=initialInc,
            maxNumInc=maxNumInc, minInc=minInc, maxInc=maxInc)

        return
    # end: def create_step

    def setup_periodic_BC(self, modelName='', strains=[]):
    # Create Periodic Buondary Conditions (PBC) on coupled edges of the RVE.
    # Args:
    # - modelName = name of the model
    # - strains = array of macro-strains
    # - forces = array of macro-forces
    # Note: if strains is not empty, then the simulation will be in strain-control. If strains is empty, and forces is 
    # not empty then the simulation will be in load-control. If they are both empty, an error will be raised.
    
        if modelName=='': modelName = self.modelName
        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myInstance = myAssembly.instances[self.partName]

        # Create faces sets
        self.faceName = ('T1-0-nodes', 'T1-1-nodes', 'T2-0-nodes', 'T2-1-nodes', 'T3-0-nodes', 'T3-1-nodes') # names of face sets
        delta = 1e-3
        faceCoord = [[0.-delta, 0.+delta, 0.+delta, 0.+delta, self.s-delta, self.t_tot-delta], [self.s-delta, 0.+delta, 0.+delta, self.s+delta, self.s-delta, self.t_tot-delta],
                    [0.+delta, 0.-delta, 0.+delta, self.s-delta, 0.+delta, self.t_tot-delta], [0.+delta, self.s-delta, 0.+delta, self.s-delta, self.s+delta, self.t_tot-delta],
                    [0.+delta, 0.+delta, 0.-delta, self.s-delta, self.s-delta, 0.+delta], [0.+delta, 0.+delta, self.t_tot-delta, self.s-delta, self.s-delta, self.t_tot+delta]]
        for i in range(len(self.faceName)):
            myAssembly.Set(name=self.faceName[i], nodes=myInstance.nodes.getByBoundingBox(faceCoord[i][0],
                faceCoord[i][1], faceCoord[i][2], faceCoord[i][3], faceCoord[i][4], faceCoord[i][5]))
        # Set PBC on periodic faces
        self.RefPointName, self.mdb = PBC.PeriodicBoundCond(self.mdb, modelName, self.faceName, self.LatticeVector, strains, self.minMeshSize*0.1, 3)
        
        # Create edges sets
        self.edgeName = ('T13-0-nodes', 'T13-1-nodes', 'T13-2-nodes', 'T13-3-nodes',
                         'T21-0-nodes', 'T21-1-nodes', 'T21-2-nodes', 'T21-3-nodes', 
                         'T32-0-nodes', 'T32-1-nodes', 'T32-2-nodes', 'T32-3-nodes') # names of edge sets
        delta = 1e-3
        edgeCoord = [[0.-delta, 0.+delta, 0.-delta, 0.+delta, self.s-delta, 0.+delta], [self.s-delta, 0.+delta, self.t_tot-delta, self.s+delta, self.s-delta, self.t_tot+delta],
                     [self.s-delta, 0.+delta, 0.-delta, self.s+delta, self.s-delta, 0.+delta], [0.-delta, 0.+delta, self.t_tot-delta, 0.+delta, self.s-delta, self.t_tot+delta],
                     [0.-delta, self.s-delta, 0.+delta, 0.+delta, self.s+delta, self.t_tot-delta], [self.s-delta, 0.-delta, 0.+delta, self.s+delta, 0.+delta, self.t_tot-delta],
                     [0.-delta, 0.-delta, 0.+delta, 0.+delta, 0.+delta, self.t_tot-delta], [self.s-delta, self.s-delta, 0.+delta, self.s+delta, self.s+delta, self.t_tot-delta],
                     [0.+delta, 0.-delta, 0.-delta, self.s-delta, 0.+delta, 0.+delta], [0.+delta, self.s-delta, self.t_tot-delta, self.s-delta, self.s+delta, self.t_tot+delta],
                     [0.+delta, 0.-delta, self.t_tot-delta, self.s-delta, 0.+delta, self.t_tot+delta], [0.+delta, self.s-delta, 0.-delta, self.s-delta, self.s+delta, 0.+delta]]
        for i in range(len(self.edgeName)):
            myAssembly.Set(name=self.edgeName[i], nodes=myInstance.nodes.getByBoundingBox(edgeCoord[i][0],
                edgeCoord[i][1], edgeCoord[i][2], edgeCoord[i][3], edgeCoord[i][4], edgeCoord[i][5]))
        # Set PBC on periodic edges
        self.mdb = PBC.PeriodicBoundCond(self.mdb, modelName, self.edgeName, self.LatticeVector_edge, strains, self.minMeshSize*0.1, 3)[1]
        
        # Create vertices sets
        self.vertName = ('T123-0-node', 'T123-1-node', 'T123-2-node', 'T123-3-node',
                         'T123-4-node', 'T123-5-node', 'T123-6-node', 'T123-7-node') # names of vertex sets
        delta = 1e-3
        vertCoord = [[0.-delta, 0.-delta, 0.-delta, 0.+delta, 0.+delta, 0.+delta], [self.s-delta, self.s-delta, self.t_tot-delta, self.s+delta, self.s+delta, self.t_tot+delta],
                     [0.-delta, 0.-delta, self.t_tot-delta, 0.+delta, 0.+delta, self.t_tot+delta], [self.s-delta, self.s-delta, 0.-delta, self.s+delta, self.s+delta, 0.+delta],
                     [0.-delta, self.s-delta, 0.-delta, 0.+delta, self.s+delta, 0.+delta], [self.s-delta, 0.-delta, self.t_tot-delta, self.s+delta, 0.+delta, self.t_tot+delta],
                     [self.s-delta, 0.-delta, 0.-delta, self.s+delta, 0.+delta, 0.+delta], [0.-delta, self.s-delta, self.t_tot-delta, 0.+delta, self.s+delta, self.t_tot+delta]]
        for i in range(len(self.vertName)):
            myAssembly.Set(name=self.vertName[i], nodes=myInstance.nodes.getByBoundingBox(vertCoord[i][0],
                vertCoord[i][1], vertCoord[i][2], vertCoord[i][3], vertCoord[i][4], vertCoord[i][5]))
        # Set PBC on periodic vertices
        self.mdb = PBC.PeriodicBoundCond(self.mdb, modelName, self.vertName, self.LatticeVector_vert, strains, self.minMeshSize*0.1, 3)[1]
        
        return
    # end: def setup_periodic_BC

    def setup_BC_strain_control(self, modelName, stepName, strains=[], moisture=0.12):
    # Set boundary conditions. After establishing whether the simulation will be in strain- or load- control, the GPS boundary
    # conditions are set (if required), and then displacements or loads are appliede to RPs. The moisture content variation is
    # set as temperature variation.
    # Args:
    # - modelName = name of the model
    # - stepName = name of the step
    # - strains = array of strains to apply as displacement to RPs. If it's empty, the simulation is in load-control.
    # - forces = array of loads to apply to RPs
    # - dw = moisture content variation

        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly

        # APPLY PERIODIC BUONDARY CONDITIONS (constraint equations)
        self.setup_periodic_BC(modelName, strains=strains)
        
        # APPLY BOUNDARY CONDITIONS TO REF. POINTS FOR PBC        
        ampName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        bcName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        hoName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        for i in range(len(self.RefPointName)):
            for j in range(len(self.RefPointName[0])):
                ampName[i][j] = 'Amp-'+self.RefPointName[i][j]
                bcName[i][j] = 'BC-'+self.RefPointName[i][j]
                hoName[i][j] = 'HO-'+self.RefPointName[i][j]

        # Set displacements
        strains = [strains[0], strains[1], strains[2], strains[3], strains[3], strains[4], strains[4], strains[5], strains[5]]
        for i in range(len(self.RefPointName)): # for every component
            for j in range(len(self.RefPointName[0])): # for every point
                if strains[i]==UNSET:
                    myModel.DisplacementBC(amplitude=UNSET, createStepName=stepName, distributionType=UNIFORM,
                        fieldName='', fixed=OFF, localCsys=None, name=bcName[i][j], region=myAssembly.sets[self.RefPointName[i][j]], u1=strains[i])
                else:
                    myModel.EquallySpacedAmplitude(name=ampName[i][j], timeSpan=STEP,
                        smooth=SOLVER_DEFAULT, fixedInterval=1, begin=0.0, data=(0.0, strains[i]))
                    myModel.DisplacementBC(amplitude=ampName[i][j], createStepName=stepName, distributionType=UNIFORM,
                        fieldName='', fixed=OFF, localCsys=None, name=bcName[i][j], region=myAssembly.sets[self.RefPointName[i][j]], u1=1.)
                # Set Houtput History (U1, RF1, CF1)
                myModel.HistoryOutputRequest(name=hoName[i][j], createStepName=stepName, variables=('U1', 'RF1', 'CF1'), 
                    region=myAssembly.sets[self.RefPointName[i][j]], sectionPoints=DEFAULT, rebar=EXCLUDE)

        """# FIX NODE TO AVOID FREE SHIFTING
        delta = self.minMeshSize*0.1
        bcName = 'Fix-Node'
        myModel.EncastreBC(create_stepName=stepName, localCsys=None, name=bcName,
            region=myAssembly.Set(name=bcName,nodes=myInstance.nodes.getByBoundingBox(self.pointsVec[0][0]-delta,
            self.pointsVec[0][1]-delta,0.,self.pointsVec[0][0]+delta,self.pointsVec[0][1]+delta,0.)))"""
        
        # Set moisture field
        fieldName = 'Moisture-field'
        myModel.Temperature(name=fieldName, createStepName='Initial', region=myAssembly.sets['Set-' + self.RVEName], 
            distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(moisture, ))#self.moistureRef, ))
        """if moisture != self.moistureRef:
            myModel.predefinedFields[fieldName].setValuesInStep(stepName=stepName, magnitudes=(moisture, ))"""

        ### SET OUTPUT FIELD
        #myModel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='Field-Outputs')
        #myModel.historyOutputRequests.changeKey(fromName='H-Output-1', toName='History-Outputs')
        if self.nlgeom == ON:
            strainRequest = 'LE'
        else:
            strainRequest = 'E'
            myModel.fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'PE', 'PEEQ', 'PEMAG', strainRequest, 'U', 'RF', 'CF', 'CSTRESS', 'CDISP',
            'IVOL', 'EVOL', 'COORD', 'TEMP'))
        
        return
    # end: def setup_BC

    def setup_BC_load_control(self, inputs, moisture, stepName, modelName, stepTime=1, creep=False):
    # Set the boundary consitions (BC).
    # Args:
    # - inputs = array of loads/strains to apply
    # - stepName = name of the step
    # - modelName = name of the model
        
        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myInstance = myAssembly.instances[self.partName]
        
        # Create steps
        nlgeom=OFF
        # Create loading step
        stepName1 = 'Load'
        loadTime = 1e-4 if creep else 1
        initialInc = loadTime if creep else 0.1
        myModel.StaticStep(name=stepName1, nlgeom=nlgeom, previous="Initial", timePeriod=loadTime,
            initialInc=initialInc, maxNumInc=100, minInc=1e-5, maxInc=loadTime)
        # Create creep step
        if creep:
            stepName2 = 'Creep'
            self.stepTime = stepTime
            initialInc = 1e-2
            myModel.StaticStep(name=stepName2, nlgeom=nlgeom, previous=stepName1, timePeriod=self.stepTime,
                initialInc=initialInc, maxNumInc=100, minInc=1e-5, maxInc=10 if self.stepTime>10 else self.stepTime)
        
        # Initialize arrays of names       
        ampName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        bcName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        hoName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
        for i in range(len(self.RefPointName)):
            for j in range(len(self.RefPointName[0])):
                ampName[i][j] = 'Amp-'+self.RefPointName[i][j]
                bcName[i][j] = 'BC-'+self.RefPointName[i][j]
                hoName[i][j] = 'HO-'+self.RefPointName[i][j]

        # Set displacements
        inputs = [inputs[0], inputs[1], inputs[2], inputs[3], inputs[3], inputs[4], inputs[4], inputs[5], inputs[5]]
        for i in range(len(self.RefPointName)): # for every component
            for j in range(len(self.RefPointName[0])): # for every point
                if inputs[i]==UNSET or inputs[i]==0.:
                    myModel.DisplacementBC(amplitude=UNSET, createStepName=stepName1, distributionType=UNIFORM,
                        fieldName='', fixed=OFF, localCsys=None, name=bcName[i][j], region=myAssembly.sets[self.RefPointName[i][j]], u1=inputs[i])
                else:
                    myModel.EquallySpacedAmplitude(name=ampName[i][j], timeSpan=STEP,
                            smooth=SOLVER_DEFAULT, fixedInterval=loadTime, begin=0.0, data=(0.0, inputs[i]))                                                                  
                    myModel.ConcentratedForce(amplitude=ampName[i][j], createStepName=stepName1, distributionType=UNIFORM,
                        field='', localCsys=None, name=bcName[i][j], region=myAssembly.sets[self.RefPointName[i][j]], cf1=1.)
                    
                # Set Houtput History (U1, RF1, CF1)
                for stepName in ([stepName1, stepName2] if creep else [stepName1]):
                    myModel.HistoryOutputRequest(name=stepName+hoName[i][j], createStepName=stepName, variables=('U1', 'RF1', 'CF1'), 
                        region=myAssembly.sets[self.RefPointName[i][j]], sectionPoints=DEFAULT, rebar=EXCLUDE)
                
        # FIX NODE TO AVOID FREE SHIFTING
        delta = self.minMeshSize
        bcName = 'Fix-Node'
        myRegion = myInstance.nodes.getByBoundingBox(self.s/2.-delta,self.s/2.-delta, self.t_tot/2.-delta, 
                                                 self.s/2.+delta, self.s/2.+delta, self.t_tot/2.+delta)
        myModel.EncastreBC(createStepName=stepName1, localCsys=None, name=bcName,
            region=myAssembly.Set(name=bcName,nodes=myRegion))
        
        # Set moisture field
        fieldName = 'Moisture-field'
        myModel.Temperature(name=fieldName, createStepName='Initial', region=myAssembly.sets['Set-' + self.partName],
            distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(moisture, ))

        ### SET OUTPUT FIELD
        strainRequest = 'LE' if self.nlgeom == ON else 'E'
        myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 
            strainRequest, 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'IVOL', 'EVOL', 'COORD', 'TEMP'))
        for stepName in ([stepName1, stepName2] if creep else [stepName1]):
            myModel.FieldOutputRequest(name='STATEV', createStepName=stepName, variables=('SDV',))
        
        
        myAssembly.regenerate()

        return
    # End: def setup_BC_load_control

    def submit_job(self, jobName, modelName, umatPath):
    # Create and submit job.
    # Args:
    # - jobName = name of the job
    # - modelName = name of the model
        
        # Remove existing .lck file
        fileExt = '.lck'
        filePath = os.getcwd() + '/' + jobName + fileExt
        if os.path.exists(filePath):
            os.remove(filePath)
        
        # Create job
        self.mdb.Job(name=jobName, model = modelName, type=ANALYSIS, atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, 
            modelPrint=OFF, multiprocessingMode=DEFAULT, nodalOutputPrecision=SINGLE, numCpus=1, queue=None, scratch='', 
            waitHours=0, waitMinutes=0, userSubroutine=umatPath, resultsFormat=ODB, numGPUs=0)

        # Run job
        myJob = self.mdb.jobs[jobName]
        myJob.submit(consistencyChecking=OFF)
        myJob.waitForCompletion()

        # Remove simulation files except odb
        path_all = glob.glob(os.getcwd() + '/' + jobName + '.*')
        extension = ('.odb','.simdir') #, '.res', '.prt', '.mdl', '.stt') # extension to keep
        path_remove = [filename for filename in path_all if not filename.endswith(extension)]
        [os.remove(filePath) for filePath in path_remove]  #"""
        
        return
    # End: def submit_job

    def copy_model(self, nameOfNewModel, nameOfModelToCopy):
    #Create a copy of a model.
    # Important: The model to copy must NOT have a step nor job already defined!
    # Args:
    # - nameOfNewModel = name of the new model
    # - nameOfModelToCopy = name of the model to copy

        self.mdb.Model(name=nameOfNewModel, objectToCopy=self.mdb.models[nameOfModelToCopy])

        return
    # end: def copy_model

    def calculate_stress_strain(self, jobName, inc):
    # Calculate the macro-stresses and -strains or the RVE at a given increment.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = increment number
    # Return:
    # - stresses = array of macr-stresses
    # - strains = array of macro-strains

        odb = openOdb(jobName + ".odb")
        myStep = odb.steps[odb.steps.keys()[-1]]
        LV = deepcopy(self.LatticeVector)

        LVcomp = [0,1,2,2,1,2,0,1,0] # ['X', 'Y', 'Z', 'YZ-Y', 'YZ-Z', 'XZ-X', 'XZ-Z', 'XY-X', 'XY-Y']

        # Determine how to handle `inc`
        if isinstance(inc, int):  
            inc = [inc]  # convert single number to list
            multi_inc = False
        elif isinstance(inc, list) and len(inc) == 0:  
            # Auto-detect available increments
            ref_node = 'Node ' + self.RefPointName[0][0].upper() + '.1'
            inc = [i[0] for i in myStep.historyRegions[ref_node].historyOutputs['CF1'].data]
            multi_inc = True

        # Initialize arrays     
        strains_all = []
        stresses_all = []
        for k in range(len(inc)):
            if len(inc) == 1: k = inc[k] 
            # Initialize arrays     
            strains = [0.]*6
            stresses = [0.]*6
            RP_F = [[0.]*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
            RP_U = [[0.]*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
            nodeName = [['']*len(self.RefPointName[0]) for i in range(len(self.RefPointName))]
            # Read data
            for i in range(len(self.RefPointName)): # for each RP
                for j in range(len(self.RefPointName[0])): # for each lattice vec.
                    nodeName[i][j] = 'Node ' + self.RefPointName[i][j].upper() + '.1'
                    myHistory = myStep.historyRegions[nodeName[i][j]]
                    # Reaction Forces (RF1)
                    myOutput = myHistory.historyOutputs['CF1']
                    if LV[j][LVcomp[i]]!=0:
                        RP_F[i][j] = -myOutput.data[k][1]/LV[j][LVcomp[i]]
                    else:
                        RP_F[i][j] = 0 #-myOutput.data[k][1]
                    # Displacement (U1)
                    myOutput = myHistory.historyOutputs['U1']
                    RP_U[i][j] =  myOutput.data[k][1]

        
            # Calculate stresses
            stresses[0] = -RP_F[0][0]/(abs(LV[1][1]*LV[2][2])) # x
            stresses[1] = -RP_F[1][1]/(abs(LV[0][0]*LV[2][2])) # y
            stresses[2] = -RP_F[2][2]/(abs(LV[1][1]*LV[0][0])) # z
            stresses[3] = -RP_F[3][2]/abs(LV[0][0]*LV[1][1]) # yz
            stresses[4] = -RP_F[5][2]/abs(LV[0][0]*LV[1][1]) # xz
            stresses[5] = -RP_F[7][1]/abs(LV[0][0]*LV[2][2]) # xy
            # Calculate strains
            strains = [RP_U[0][0], RP_U[1][0], RP_U[2][0], (RP_U[3][2]+RP_U[4][1]),
                    (RP_U[5][2]+RP_U[6][0]), (RP_U[7][1]+RP_U[8][0])] # shear strains = engineering strains
            
            # Store results
            strains_all.append(strains)
            stresses_all.append(stresses)
        odb.close()

        # Return a single array if there was only one increment, otherwise return lists
        if not multi_inc:
            return stresses_all[0], strains_all[0], inc[0]
        else:
            return stresses_all, strains_all, inc
    # end: def calculate_stress_strain

    def calculate_compliance_matrix(self, jobName, inc=[]):
    # Calculate 3D compliance matrix from 9 elementary cases results.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - StiffMat = stiffness matrix of the laminate composite

        numElCase = 9 # Number of elementary cases
        # Initialize arrays
        stress_mat = []
        strain_mat = []
        inc_mat = []
        self.CompMat = [[None for _ in range(6)] for _ in range(6)]
    
        # Get stresses, strains, and inc_val from odb files for each elementary case.
        for i in range(numElCase):
            stress, strain, inc_val = self.calculate_stress_strain(
                jobName=jobName + "-{}".format(i+1), inc=inc)
            stress_mat.append(stress)
            strain_mat.append(strain)
            inc_mat.append(inc_val)
        
        # Calculate compliances
        if inc: # inc is provided as a number (nonempty)
            # Diagonal elements
            for i in range(6):
                stress = stress_mat[i]
                strain = strain_mat[i]
                self.CompMat[i][i] = strain[i]/stress[i]
            # Out-of-diagonal elements
            k = [6, 7, 8]
            i = [0, 0, 1]
            j = [1, 2, 2]
            for n in range(3):
                stress = stress_mat[k[n]]
                strain = strain_mat[k[n]]   
                Sij = (strain[i[n]]-self.CompMat[i[n]][i[n]]*stress[i[n]])/stress[j[n]]
                Sji = (strain[j[n]]-self.CompMat[j[n]][j[n]]*stress[j[n]])/stress[i[n]]           
                self.CompMat[i[n]][j[n]] = (Sij+Sji)/2.
                self.CompMat[j[n]][i[n]] = self.CompMat[i[n]][j[n]]
            result = {
                'time': inc,
                'D11': self.CompMat[0][0],
                'D22': self.CompMat[1][1],
                'D33': self.CompMat[2][2],
                'D44': self.CompMat[3][3],
                'D55': self.CompMat[4][4],
                'D66': self.CompMat[5][5],
                'D12': self.CompMat[0][1],
                'D13': self.CompMat[0][2],
                'D23': self.CompMat[1][2]
            }

        else: # inc==[] so that each elementary case returns arrays (one per inc)
            
            # Diagonal elements
            #print(stress_mat)
            for i in range(6):
                inc_val = inc_mat[i]
                n_inc = len(inc_val)
                self.CompMat[i][i] = [strain_mat[i][k][i] / stress_mat[i][k][i] for k in range(1,n_inc)]
            # Out-of-diagonal elements
            k = [6, 7, 8]
            i = [0, 0, 1]
            j = [1, 2, 2]
            for n in range(3):
                comp_list = []
                inc_val = inc_mat[k[n]]
                n_inc = len(inc_val)
                for inc_i in range(1,n_inc):
                    stress = stress_mat[k[n]]
                    strain = strain_mat[k[n]]
                    Sij = (strain[inc_i][i[n]] -
                        self.CompMat[i[n]][i[n]][inc_i-1] * stress[inc_i][i[n]]) / stress[inc_i][j[n]]
                    Sji = (strain[inc_i][j[n]] -
                        self.CompMat[j[n]][j[n]][inc_i-1] * stress[inc_i][j[n]]) / stress[inc_i][i[n]]
                    comp_list.append((Sij + Sji) / 2.)
                self.CompMat[i[n]][j[n]] = comp_list
                self.CompMat[j[n]][i[n]] = comp_list
                #inc_val = inc_val[1:]

            # Create the dictionary to return.
            result = {
                'time11': inc_mat[0][1:],
                'D11': self.CompMat[0][0],
                'time22': inc_mat[1][1:],
                'D22': self.CompMat[1][1],
                'time33': inc_mat[2][1:],
                'D33': self.CompMat[2][2],
                'time44': inc_mat[3][1:],
                'D44': self.CompMat[3][3],
                'time55': inc_mat[4][1:],
                'D55': self.CompMat[4][4],
                'time66': inc_mat[5][1:],
                'D66': self.CompMat[5][5],
                'time12': inc_mat[6][1:],
                'D12': self.CompMat[0][1],
                'time13': inc_mat[7][1:],
                'D13': self.CompMat[0][2],
                'time23': inc_mat[8][1:],
                'D23': self.CompMat[1][2]
            }
        return result
    # End: def calculate_stiffness_matrix

    def compute_3D_compliance_matrix(self, val, moisture, run, umatPath, stepTime=1, creepActivate=False):
    # Compute 3D stiffness matrix by running 9 elementary cases.
    # Args:
    # - val = strain to apply for each elementary case
    # - run =
    # Returns:
    # - StiffMat = stiffness matrix of the laminate composite

        # Set elementary cases
        val_x =  [val[0],   UNSET, UNSET, 0., UNSET, UNSET, val[0],   val[0],   UNSET]
        val_y =  [UNSET, val[1],   UNSET, UNSET, 0., UNSET, val[1],   UNSET, val[1]]
        val_z =  [UNSET, UNSET, val[2],   UNSET, UNSET, 0., UNSET, val[2],   val[2]]
        val_yz = [UNSET, UNSET, UNSET, val[3],   UNSET, UNSET, UNSET, UNSET, UNSET]
        val_xz = [UNSET, UNSET, UNSET, UNSET, val[4],   UNSET, UNSET, UNSET, UNSET]
        val_xy = [UNSET, UNSET, UNSET, UNSET, UNSET, val[5],   UNSET, UNSET, UNSET]

        # Activate/deactivate creep in umat
        parent_path = os.path.dirname(umatPath)
        filepath = os.path.join(parent_path, "activate_ve_creep.inc")
        self.update_umat_inc_file(creepActivate, filepath)

        if run:
            # Copy model
            numElCases = 9
            for i in range(numElCases-1):
                nameOfNewModel = 'Model-{}'.format(i+2)
                self.copy_model(nameOfNewModel, self.modelName)
        
            # Run elementary cases
            for k in range(numElCases):
                modelName = 'Model-{}'.format(k+1)
                stepName = 'Loading-{}'.format(k+1)
                jobName = 'Job-{}'.format(k+1)
                inputs = [val_x[k], val_y[k], val_z[k], val_yz[k], val_xz[k], val_xy[k]]
                # Set boundary conditions
                self.setup_periodic_BC(modelName, inputs)
                self.setup_BC_load_control(inputs, moisture, stepName, modelName, stepTime, creepActivate)
                # Run job
                self.submit_job(jobName, modelName, umatPath)
        else:
            # Build model to store variables and names
            k = 0
            modelName = 'Model-{}'.format(k+1)
            stepName = 'Loading-{}'.format(k+1)
            jobName = 'Job-{}'.format(k+1)
            inputs = [val_x[k], val_y[k], val_z[k], val_yz[k], val_xz[k], val_xy[k]]
            # Set boundary conditions
            self.setup_periodic_BC(modelName, inputs)
            self.setup_BC_load_control(inputs, moisture, stepName, modelName, stepTime, creepActivate)

        # Calculate stiffness matrix
        jobName = 'Job'
        stepName = 'Loading'
        inc = [] if creepActivate else -1
        comp_coeff = self.calculate_compliance_matrix(jobName, inc)
        """self.StiffMat = np.linalg.inv(self.CompMat)
        self.calculate_engineering_constants()"""

        return comp_coeff
    # End: def compute_3D_compliance_matrix

    def calculate_eng_const(self):
    # Calculate the engineering constants of the RVE from simulations results.
    # Args: none
    # Return:
    # - RVEEngConst = array of engineering constants calculated from simulations results
    #                 [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12]

        self.complianceMatrix = self.CompMat#np.linalg.inv(self.stiffnessMatrix)
        

        E1 = 1./self.complianceMatrix[0][0]
        E2 = 1./self.complianceMatrix[1][1]
        E3 = 1./self.complianceMatrix[2][2]
        G23 = 1./self.complianceMatrix[3][3]
        G13 = 1./self.complianceMatrix[4][4]
        G12 = 1./self.complianceMatrix[5][5]
        v23 = -self.complianceMatrix[2][1]*E2
        v13 = -self.complianceMatrix[2][0]*E1
        v12 = -self.complianceMatrix[1][0]*E1

        v32 = -self.complianceMatrix[1][2]*E3
        v31 = -self.complianceMatrix[0][2]*E3
        v21 = -self.complianceMatrix[0][1]*E2

        self.RVEEngConst = [E1,E2,E3,G23,G13,G12,v12,v13,v23,v21,v31,v32]

        return self.RVEEngConst#, self.complianceMatrix
    # end: def calculate_eng_const

    def update_umat_inc_file(self, creepActivate, filepath):
        if not os.path.isfile(filepath):
            raise FileNotFoundError("The file '{}' does not exist.".format(filepath))

        value = 1 if creepActivate else 0

        with open(filepath, "r") as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if "KW_ViscoIdent" in line:
                lines[i] = "      KW_ViscoIdent = {}\n".format(value)
                break

        with open(filepath, "w") as file:
            file.writelines(lines)
    # end: def update_umat_inc_file

