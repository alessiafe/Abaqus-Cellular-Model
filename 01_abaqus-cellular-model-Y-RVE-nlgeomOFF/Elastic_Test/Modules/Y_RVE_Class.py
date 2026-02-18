# Y_RVE_Class.py
# This script implements a class for 3D Y-shaped RVE in Abaqus.
# Author Alessia Ferrara, September 2021
# Last modified June 2025

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

class Y_RVE_Class:

    def __init__(self):
    # Constructor of the Y_RVE_Class.
    # Args: none

        # Model Data Base construction
        
        self.modelName = 'Model-1'
        self.RVEName = 'RVE'
        self.layersName = ["ML", "P", "S1+", "S1-", "S2", "S3"]

        # Geometry Parameters initialization
        self.dim = 3.
        self.h, self.l, self.theta, self.t, self.b = 0.,0.,0.,0.,0.
        self.x1, self.y1, self.x2, self.y2, self.x3, self.y3 = 0.,0.,0.,0.,0.,0.
        self.x4, self.y4, self.x10, self.y10, self.x6, self.y6 = 0.,0.,0.,0.,0.,0.
        self.x7, self.y7, self.x5, self.y5, self.x8, self.y8, self.x9, self.y9 = 0.,0.,0.,0.,0.,0.,0.,0.
        self.numLayer = 0.
        self.lt = []
        self.LatticeVec = []
        self.numPerBound = 0.
        self.pointsVec = np.zeros((13,2))
        self.pointsVec_1 = np.zeros((21,2))
        self.pointsVec_2 = np.zeros((21,2))
        self.pointsVec_3 = np.zeros((21,2))

        # Material variables
        self.alpha = [0., 0., 60., -60., 15., 75.]
        self.depvar = 126
        self.moistureRef = 0.12
        self.phi = []
        self.stiff_mat = []
        self.hygro = []
        self.stiffVec = np.zeros((1,10))

        # Reference points initialization
        self.RefPointName = []
        
        # Element features initialization
        self.nlgeom=OFF

        return
    # end: def _init_

    def setup_single_geometry(self, l, h, theta, t, b):
    # The geometry of the RVE is created according to the parametrisation scheme (refer to the documentation for the meaning of the
    # parameters). The function creates the  parts in 'Model-1'. A part is created for each layer of the cell walls, then they are
    # merged in one single assembly.
    # Args:
    # - l, h, theta, lt, b = geometry parameters of the RVE according to the parametrisation scheme
    # Return:
    # - modelName = name of the model

        self.mdb = Mdb()
        myModel = self.mdb.models[self.modelName]
        myAssembly = myModel.rootAssembly
        self.h = copy(h) # h
        self.l = copy(l) # l
        self.theta = copy(theta) # theta
        self.t = t # double-wall thickness
        self.b = copy(b) # longitudinal thickness
        
        # Vertices coordinates
        self.x0 = 0
        self.y0 = 0.
        self.x1 = 0.
        self.y1 = -h
        self.x2 = t/2.
        self.y2 = -h + t/2. * (1./cos(theta) - tan(theta))
        self.x3 = t/2.
        self.y3 = t/2. * (tan(theta) - 1./cos(theta))
        self.x4 = l * cos(theta)
        self.y4 = l * sin(theta) - t/(2.*cos(theta))
        self.x5 = l * cos(theta)
        self.y5 = l * sin(theta)
        self.x6 = l * cos(theta) - t/2.
        self.y6 = l * sin(theta) + t/2. * (1/cos(theta) - tan(theta))
        self.x7 = 0.
        self.y7 = t/(2. * cos(theta))
        self.x8 = -l * cos(theta) + t/2.
        self.y8 = l * sin(theta) + t/2. * (1/cos(theta) - tan(theta))
        self.x9 = -l * cos(theta)
        self.y9 = l * sin(theta)
        self.x10 = -l * cos(theta)
        self.y10 = l * sin(theta) - t/(2.*cos(theta))
        self.x11 = - t/2.
        self.y11 = t/2. * (tan(theta) - 1./cos(theta))
        self.x12 = -t/2.
        self.y12 = -h + t/2. * (1./cos(theta) - tan(theta))

        # Array of vertices coordinates taking into account geometry imperfection (horizontal translation of central triangle)
        self.pointsVec[0] = (self.x0, self.y0)
        self.pointsVec[1] = (self.x1, self.y1)
        self.pointsVec[2] = (self.x2, self.y2)
        self.pointsVec[3] = (self.x3+self.x0, self.y3)
        self.pointsVec[4] = (self.x4, self.y4)
        self.pointsVec[5] = (self.x5, self.y5)
        self.pointsVec[6] = (self.x6, self.y6)
        self.pointsVec[7] = (self.x7+self.x0, self.y7)
        self.pointsVec[8] = (self.x8, self.y8)
        self.pointsVec[9] = (self.x9, self.y9)
        self.pointsVec[10] = (self.x10, self.y10)
        self.pointsVec[11] = (self.x11+self.x0, self.y11)
        self.pointsVec[12] = (self.x12, self.y12)

        # Lattice Vector
        self.LatticeVec = [[self.x5 - self.x9, self.y5 - self.y9, 0], # T1
                           [self.x9 - self.x1, self.y9 - self.y1, 0], # T2
                           [self.x5 - self.x1, self.y5 - self.y1, 0]]  # T3
        self.numPerBound = len(self.LatticeVec) # number of periodic couples
    
        # Create part
        partName = self.RVEName
        myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        for k in range(1,len(self.pointsVec)-1):
            myModel.sketches['__profile__'].Line(point1=self.pointsVec[k], point2=self.pointsVec[k+1])
        myModel.sketches['__profile__'].Line(point1=self.pointsVec[-1], point2=self.pointsVec[1])
        myModel.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        myModel.parts[partName].BaseSolidExtrude(sketch=myModel.sketches['__profile__'], depth=self.b)
        del myModel.sketches['__profile__']

        # Create set
        setName = 'Set-' + partName
        self.create_set(partName, setName)

        # Create Assembly
        myAssembly.Instance(name=partName, part=myModel.parts[partName], dependent=ON)

        # Create reference point
        self.LatticeVec = self.LatticeVec + [[0, 0, self.b]] # TL
        self.numPerBound = len(self.LatticeVec) # number of periodic couples
        
        return self.modelName
    # end: def setup_single_geometry

    def setup_multilayers_geometry(self, l, h, theta, lt, b):
    # The geometry of the RVE is created according to the parametrisation scheme (refer to the documentation for the meaning of the
    # parameters). The function creates the  parts in 'Model-1'. A part is created for each layer of the cell wall, then they are
    # merged in one single assembly.
    # Args:
    # - l, h, theta, lt, b = geometry parameters of the RVE according to the parametrisation scheme
    # Return:
    # - modelName = name of the model

        self.mdb = Mdb()
        myModel = self.mdb.models[self.modelName]
        self.elementType = ''
        self.h = copy(h) # h
        self.l = copy(l) # l
        self.theta = copy(theta) # theta
        self.lt = copy(lt) # array of layers thickness
        self.numLayer = len(self.lt) # number of layers
        self.t = t = 2*sum(self.lt) # double-wall thickness
        self.b = copy(b) # longitudinal thickness

        # Vertices coordinates
        self.x0 = 0
        self.y0 = 0.
        self.x1 = 0.
        self.y1 = -h
        self.x2 = t/2.
        self.y2 = -h + t/2. * (1./cos(theta) - tan(theta))
        self.x3 = t/2.
        self.y3 = t/2 * (tan(theta) - 1./cos(theta))
        self.x4 = l * cos(theta)
        self.y4 = l * sin(theta) - t/(2. * cos(theta))
        self.x5 = l * cos(theta)
        self.y5 = l * sin(theta)
        self.x6 = l * cos(theta) - t/2.
        self.y6 = l * sin(theta) + t/2. * (1/cos(theta) - tan(theta))
        self.x7 = 0.
        self.y7 = t / (2. * cos(theta))
        self.x8 = -l * cos(theta) + t/2.
        self.y8 = l * sin(theta) + t/2. * (1/cos(theta) - tan(theta))
        self.x9 = -l * cos(theta)
        self.y9 = l * sin(theta)
        self.x10 = -l * cos(theta)
        self.y10 = l * sin(theta) - t / (2. * cos(theta))
        self.x11 = - t / 2
        self.y11 = t / 2 * (tan(theta) - 1 / cos(theta))
        self.x12 = -t / 2.
        self.y12 = -h + t / 2. * (1 / cos(theta) - tan(theta))

        # Array of vertices coordinates taking into account geometry imperfection (horizontal translation of central triangle)
        self.pointsVec[0] = (self.x0, self.y0)
        self.pointsVec[1] = (self.x1, self.y1)
        self.pointsVec[2] = (self.x2, self.y2)
        self.pointsVec[3] = (self.x3+self.x0, self.y3)
        self.pointsVec[4] = (self.x4, self.y4)
        self.pointsVec[5] = (self.x5, self.y5)
        self.pointsVec[6] = (self.x6, self.y6)
        self.pointsVec[7] = (self.x7+self.x0, self.y7)
        self.pointsVec[8] = (self.x8, self.y8)
        self.pointsVec[9] = (self.x9, self.y9)
        self.pointsVec[10] = (self.x10, self.y10)
        self.pointsVec[11] = (self.x11+self.x0, self.y11)
        self.pointsVec[12] = (self.x12, self.y12)

        # Lattice Vector
        self.LatticeVec = [[self.x5 - self.x9, self.y5 - self.y9, 0], # T1
                           [self.x9 - self.x1, self.y9 - self.y1, 0], # T2
                           [self.x5 - self.x1, self.y5 - self.y1, 0]]  # T3
        self.LatticeVec = self.LatticeVec + [[0, 0, self.b]]
        self.numPerBound = len(self.LatticeVec) # number of periodic couples
    
        ### Cell-1
        self.pointsVec_1[0] = (self.x0, self.y0)
        self.pointsVec_1[1] = (self.x1, self.y1)
        dir_x, dir_y = self.x2-self.x1, self.y2-self.y1
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(2,7):
            dist = sum(self.lt[0:i-1])/dir_xn
            self.pointsVec_1[i] = (self.x1+dist*dir_xn, self.y1+dist*dir_yn)
        self.pointsVec_1[7] = (self.x2, self.y2)
        self.pointsVec_1[8] = (self.x3 + self.x0, self.y3)
        self.pointsVec_1[9] = (self.x4, self.y4)
        dir_n = self.y5-self.y4
        dir_yn = self.t/2 / dir_n
        for i in range(10,15):
            dist = sum(self.lt[0:i-9])/dir_yn
            self.pointsVec_1[24-i] = (self.x4, self.y5-dist)
            self.pointsVec_1[15] = (self.x5, self.y5)
        dir_x, dir_y = self.x3, self.y3
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(16,21):
            dist = sum(self.lt[0:i-15])/dir_xn
            self.pointsVec_1[i] = (dist*dir_xn + self.x0, dist*dir_yn)
        partName = 'Cell-1'
        self.sketch_cell_wall(partName, *self.pointsVec_1)
        # Create external surface
        myFace1 = myModel.parts[partName+'-6'].faces.findAt((((self.x3+self.x4)/2.,(self.y3+self.y4)/2.,self.b/2.),),)
        myFace2 = myModel.parts[partName+'-6'].faces.findAt(((self.x2,(self.y2+self.y3)/2.,self.b/2.),),)
        myModel.parts[partName+'-6'].Surface(side1Faces=[myFace1, myFace2], name='Surf-1')

        ### Cell-2
        self.pointsVec_2[0] = (self.x0, self.y0)
        self.pointsVec_2[1] = (self.x5, self.y5)
        dir_x, dir_y = self.x6-self.x5, self.y6-self.y5
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(2,7):
            dist = sum(self.lt[0:i-1])/dir_xn
            self.pointsVec_2[i] = (self.x5-dist*dir_xn, self.y5-dist*dir_yn)
        self.pointsVec_2[7] = (self.x6, self.y6)
        self.pointsVec_2[8] = (self.x7 + self.x0, self.y7)
        self.pointsVec_2[9] = (self.x8, self.y8)
        dir_x, dir_y = self.x9-self.x8, self.y9-self.y8
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(10,15):
            dist = sum(self.lt[0:i-9])/dir_xn
            self.pointsVec_2[24-i] = (self.x9+dist*dir_xn, self.y9+dist*dir_yn)
            self.pointsVec_2[15] = (self.x9, self.y9)
        dir_n = self.y7
        dir_yn = self.t/2 / dir_n
        for i in range(16,21):
            dist = sum(self.lt[0:i-15])/dir_yn
            self.pointsVec_2[i] = (self.x0, dist)
        partName = 'Cell-2'
        self.sketch_cell_wall(partName, *self.pointsVec_2)
        # Create external surface
        myFace1 = myModel.parts[partName+'-6'].faces.findAt((((self.x7+self.x6)/2.,(self.y7+self.y6)/2.,self.b/2.),),)
        myFace2 = myModel.parts[partName+'-6'].faces.findAt((((self.x7+self.x8)/2.,(self.y7+self.y8)/2.,self.b/2.),),)
        myModel.parts[partName+'-6'].Surface(side1Faces=[myFace1, myFace2], name='Surf-2')

        ### Cell-3
        self.pointsVec_3[0] = (self.x0, self.y0)
        self.pointsVec_3[1] = (self.x1, self.y1)
        dir_x, dir_y = self.x12-self.x1, self.y12-self.y1
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(2,7):
            dist = sum(self.lt[0:i-1])/dir_xn
            self.pointsVec_3[i] = (self.x1-dist*dir_xn, self.y1-dist*dir_yn)
        self.pointsVec_3[7] = (self.x12, self.y12)
        self.pointsVec_3[8] = (self.x11 + self.x0, self.y11)
        self.pointsVec_3[9] = (self.x10, self.y10)
        dir_n = self.y9-self.y10
        dir_yn = self.t/2 / dir_n
        for i in range(10,15):
            dist = sum(self.lt[0:i-9])/dir_yn
            self.pointsVec_3[24-i] = (self.x10, self.y9-dist)
            self.pointsVec_3[15] = (self.x9, self.y9)
        dir_x, dir_y = self.x11, self.y11
        dir_n = (dir_x**2 + dir_y**2)**0.5
        dir_xn, dir_yn = dir_x / dir_n, dir_y / dir_n
        for i in range(16,21):
            dist = sum(self.lt[0:i-15])/dir_xn
            self.pointsVec_3[i] = (-dist*dir_xn + self.x0, -dist*dir_yn)
        partName = 'Cell-3'
        self.sketch_cell_wall(partName, *self.pointsVec_3)
        # Create external surface
        myFace1 = myModel.parts[partName+'-6'].faces.findAt((((self.x11+self.x10)/2.,(self.y11+self.y10)/2.,self.b/2.),),)
        myFace2 = myModel.parts[partName+'-6'].faces.findAt(((self.x12,(self.y12+self.y11)/2.,self.b/2.),),)
        myModel.parts[partName+'-6'].Surface(side1Faces=[myFace1, myFace2], name='Surf-3')

        # Create set for each Cell
        setName = []
        for k in range(3): # for each cell
            for i in range(len(self.layersName)): # for each layer
                partName = 'Cell-'+str(k+1)+'-'+str(i+1)
                setName.append('Set-'+ self.layersName[i] +str(k+1))
                self.create_set(partName, setName[-1])

        # Create assembly
        myAssembly = myModel.rootAssembly
        instMerge = [''] * (self.numLayer*3)
        for k in range(3): # for each cell
            for i in range(self.numLayer): # for each layer
                instName = 'Cell-'+str(k+1)+'-'+str(i+1)+'-'+str(1)
                instMerge[i+6*k] = myAssembly.instances[instName]
        myAssembly.InstanceFromBooleanMerge(name=self.RVEName, instances=(instMerge),
            keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
        myAssembly.features.changeKey(fromName=self.RVEName+'-1', toName=self.RVEName)

        # Create RVE set
        myAssembly.SetByBoolean(name='Set-' + self.RVEName, sets=tuple(myAssembly.allInstances[self.RVEName].sets[set] for set in setName))

        return self.modelName
    # end: def setup_multilayers_geometry

    def sketch_cell_wall(self, partName, *pointsVec):
    # Sketch layers of the cell wall.
    # Args:
    # - partName = name of the part
    # - *pointsVec = pointer to a vector of RVE coordinates

        self.sketch_layer(partName+'-1', *[pointsVec[0], pointsVec[1], pointsVec[2], pointsVec[16], pointsVec[14], pointsVec[15]]) # layer M
        self.sketch_layer(partName+'-2', *[pointsVec[16], pointsVec[2], pointsVec[3], pointsVec[17], pointsVec[13], pointsVec[14]]) # layer P
        self.sketch_layer(partName+'-3', *[pointsVec[17], pointsVec[3], pointsVec[4], pointsVec[18], pointsVec[12], pointsVec[13]]) # layer S1+
        self.sketch_layer(partName+'-4', *[pointsVec[18], pointsVec[4], pointsVec[5], pointsVec[19], pointsVec[11], pointsVec[12]]) # layer S1-
        self.sketch_layer(partName+'-5', *[pointsVec[19], pointsVec[5], pointsVec[6], pointsVec[20], pointsVec[10], pointsVec[11]]) # layer S2
        self.sketch_layer(partName+'-6', *[pointsVec[20], pointsVec[6], pointsVec[7], pointsVec[8], pointsVec[9], pointsVec[10]]) # layer S3

        return
    # end: def sketch_cell_wall

    def sketch_layer(self, partName, *sketchPoints):
    # Create sketch by connecting 6 given points.
    # Args:
    # - partName = name of the part
    # - *sketchPoints = pointer to a vector of 6 points coordinates
    
        myModel = self.mdb.models[self.modelName]
        myAssembly = myModel.rootAssembly

        myModel.ConstrainedSketch(name='__profile__', sheetSize=max(abs(self.x5-self.x9), abs(self.y6-self.y1)))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[0]), point2=(sketchPoints[1]))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[1]), point2=(sketchPoints[2]))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[2]), point2=(sketchPoints[3]))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[3]), point2=(sketchPoints[4]))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[4]), point2=(sketchPoints[5]))
        myModel.sketches['__profile__'].Line(point1=(sketchPoints[5]), point2=(sketchPoints[0]))

        myModel.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        myModel.parts[partName].BaseSolidExtrude(sketch=myModel.sketches['__profile__'], depth=self.b)
        del myModel.sketches['__profile__']

        myAssembly.Instance(name=partName+'-1', part=myModel.parts[partName], dependent=ON)

        return
    # end: def sketch_layer

    def create_set(self, partName, setName):
    # Create a set of the given part.
    # Args:
    # - partName: name of the part
    # - setName: name of the set

        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[partName]
        myModel.parts[partName].Set(cells=myPart.cells, name=setName)

        return
    # end: def create_set

    def create_layers_material(self, alpha=[]):
    # Create and assign materials to composite plies.
    # Args:
    # - materialStiffMat = array of ply stiffness matrices
    # - alpha = array of ply orientations

        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[self.RVEName]
        if alpha: self.alpha = alpha
        
        # Setup materials
        for i in range(len(self.layersName)):
            # Create material
            myModel.Material(name=self.layersName[i][:2])
            myMaterial = myModel.materials[self.layersName[i][:2]]
            myMaterial.Depvar(n=self.depvar)
            myMaterial.UserMaterial(mechanicalConstants=(0.0, ))
            # Create section
            sectionName = self.layersName[i]
            myModel.HomogeneousSolidSection(name=sectionName, material=self.layersName[i][:2], thickness=self.b)
            for k in range(3): # for each cell
                # Assign material orientation
                setName = 'Set-'+ self.layersName[i] +str(k+1)
                myRegion = myPart.sets[setName]
                normalAxisRegion = myPart.surfaces['Surf-{}'.format(k+1)]
                myPart.MaterialOrientation(region=myRegion, orientationType=DISCRETE, axis=AXIS_3,
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, flipNormalDirection=False, normalAxisDirection=AXIS_3,
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 0.0, 1.0), primaryAxisDirection=AXIS_1, flipPrimaryDirection=False, angle=self.alpha[i])
                # Assign section
                myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                    region=myPart.sets[setName], sectionName=sectionName, thicknessAssignment=FROM_SECTION)
        
        myAssembly = myModel.rootAssembly        
        myAssembly.regenerate()

        return
    # End: def create_layers_material

    def setup_isotropic_material(self, materialName='TEST', engconst=[1., 0.3]):
    # Create materials and assign material orientations to each layer of the RVE.
    # Args:
    # - stiff = array of stiffness matrices of layers
    # - phi = array of MFAs of layers
    # - hygro = array of hygroexpansion coefficients of layers
    # - thickness = section thickness

        myModel = self.mdb.models[self.modelName]
        myAssembly = myModel.rootAssembly
        myPart = myModel.parts[self.RVEName]
        self.isotropicEngConst = deepcopy(engconst)
        
        # Create material
        myModel.Material(name=materialName)
        myMaterial = myModel.materials[materialName]
        #myMaterial.Elastic(type=ISOTROPIC, table=((self.isotropicEngConst[0], self.isotropicEngConst[1]), ))
        myMaterial.Depvar(n=self.depvar)
        myMaterial.UserMaterial(mechanicalConstants=(0.0, ))

        # Create section
        sectionName = materialName
        myModel.HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)
        # Assign section
        setName = 'Set-' + self.RVEName
        myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                region=myPart.sets[setName], sectionName=sectionName, thicknessAssignment=FROM_SECTION)

        myAssembly.regenerate()

        return
    # end: def setup_isotropic_material

    def setup_mesh(self, meshSize, minSizeFactor):
    # Set up the mesh of the part.
    # Args:
    # - meshSize = global size of the mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the global element size

        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[self.RVEName]
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
        myInstance = myAssembly.instances[self.RVEName]
        
        # CREATE NODES SETS FOR EACH EDGE TO COUPLE
        delta = self.minMeshSize*0.1 # tolerance to find nodes (0.1 * min distance between nodes)
        periodNodeName = ('T1-0-nodes', 'T1-1-nodes', 'T2-0-nodes', 'T2-1-nodes', 'T3-0-nodes', 'T3-1-nodes') # sets names
        coordIdx = [10, 10, 9, 9, 4, 4, 5, 5, 1, 1, 2, 2, 9, 9, 8, 8, 12, 1, 1, 12, 6, 5, 5, 6]
        for i in range(len(coordIdx)/4):
            xmin, ymin, zmin = self.pointsVec[coordIdx[4*i]][0]-delta,self.pointsVec[coordIdx[4*i+1]][1]-delta,-delta
            xmax, ymax, zmax = self.pointsVec[coordIdx[4*i+2]][0]+delta,self.pointsVec[coordIdx[4*i+3]][1]+delta,self.b+delta
            face = myInstance.faces.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
            myAssembly.Set(name=periodNodeName[i], nodes=myAssembly.Set(name=periodNodeName[i], faces=face).nodes)
        
        
        periodNodeName_temp = ('TL-0-nodes', 'TL-1-nodes')
        coordZ = [0, self.b]
        for i in range(len(coordZ)):
            myAssembly.Set(name=periodNodeName_temp[i], nodes=myInstance.nodes.getByBoundingBox(
                self.pointsVec[10][0]-delta,self.pointsVec[1][1]-delta,coordZ[i]-delta,
                self.pointsVec[4][0]+delta,self.pointsVec[6][1]+delta,coordZ[i]+delta))
            myAssembly.SetByBoolean(name=periodNodeName_temp[i], operation=DIFFERENCE,
                sets=([myAssembly.sets[periodNodeName_temp[i]]]+[myAssembly.sets[periodNodeName[k]] for k in range(len(periodNodeName))]))

        periodNodeName = periodNodeName + periodNodeName_temp

        centerNodeName = ('A-nodes', 'B-nodes', 'C-nodes')
        coordIdx = [9, 5, 1]
        for i in range(len(centerNodeName)):
            center1 = (self.pointsVec[coordIdx[i]][0], self.pointsVec[coordIdx[i]][1], -delta)
            center2 = (self.pointsVec[coordIdx[i]][0], self.pointsVec[coordIdx[i]][1], self.b+delta)
            nodes = myInstance.nodes.getByBoundingCylinder(center1, center2, delta)
            myAssembly.Set(name=centerNodeName[i], nodes=nodes)
        for i in range(len(periodNodeName)-2):
            myAssembly.SetByBoolean(name=periodNodeName[i], operation=DIFFERENCE,
                sets=([myAssembly.sets[periodNodeName[i]]]+[myAssembly.sets[centerNodeName[k]] for k in range(len(centerNodeName))]))

        # SET PERIODIC BOUNDARY CONDITIONS
        self.RefPointName, self.mdb = PBC.PeriodicBoundCond(self.mdb, modelName, periodNodeName, self.LatticeVec, strains, self.minMeshSize*0.1, self.dim)
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
        myInstance = myAssembly.instances[self.RVEName]
        
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
                hoName[i][j] = '-HO-'+self.RefPointName[i][j]

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
        myModel.EncastreBC(createStepName=stepName1, localCsys=None, name=bcName,
            region=myAssembly.Set(name=bcName,nodes=myInstance.nodes.getByBoundingBox(self.pointsVec[0][0]-delta,
            self.pointsVec[0][1]-delta,self.b/2.-delta,self.pointsVec[0][0]+delta,self.pointsVec[0][1]+delta,self.b/2.+delta)))
        
        # Set moisture field
        fieldName = 'Moisture-field'
        myModel.Temperature(name=fieldName, createStepName='Initial', region=myAssembly.sets['Set-' + self.RVEName],
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
        [os.remove(filePath) for filePath in path_remove]
        
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
        LV = deepcopy(self.LatticeVec)
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
                for j in range(len(self.RefPointName[0])): # for each lattice vector
                    nodeName[i][j] = 'Node ' + self.RefPointName[i][j].upper() + '.1'
                    myHistory = myStep.historyRegions[nodeName[i][j]]
                    # Reaction Forces (RF1)
                    myOutput = myHistory.historyOutputs['CF1']
                    if LV[j][LVcomp[i]]!=0:
                        RP_F[i][j] = -myOutput.data[k][1]/LV[j][LVcomp[i]]
                    else:
                        RP_F[i][j] = -myOutput.data[k][1]
                    # Displacement (U1)
                    myOutput = myHistory.historyOutputs['U1']
                    RP_U[i][j] =  myOutput.data[k][1]

            # Calculate stresses
            stresses[0] = (RP_F[0][0]+RP_F[0][1])/(abs(LV[1][1]*LV[3][2])) # x
            stresses[1] = -RP_F[1][2]/(abs(LV[0][0])*LV[3][2]) # y
            stresses[2] = -RP_F[2][3]/(abs(LV[1][1]*LV[0][0])) # z
            stresses[3] = -RP_F[3][3]/(abs(LV[0][0]*LV[1][1])) # yz       
            stresses[4] = -RP_F[5][3]/(abs(LV[0][0]*LV[1][1])) # xz
            stresses[5] = -(RP_F[7][1]+RP_F[7][2])/(2*abs(LV[0][0])*LV[3][2]) # xy
            # Get strains
            strains = [RP_U[0][0], RP_U[1][1], RP_U[2][3], (RP_U[4][1]+RP_U[3][2]), (RP_U[6][1]+RP_U[5][2]), (RP_U[7][2]+RP_U[8][1])]

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

        self.RVEEngConst = [E1,E2,E3,G23,G13,G12,v23,v13,v12,v32,v31,v21]

        return self.RVEEngConst#, self.complianceMatrix
    # end: def calculate_eng_const

    def calculate_Ashby_eng_const(self):
    # Calculate the engineering constants of the RVE from from Gibson-Ashby 
    # equations for hexagonal honeycomb.
    # Args: none
    # Return:
    # - AshbyEngConst = array of engineering constants calculated by applying Gibson-Ashby
    #                   equations [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12]

        # In-plane-deformations
        E1 = self.isotropicEngConst[0] * (self.t/self.l)**3 * cos(self.theta)/((self.h/self.l+sin(self.theta))*sin(self.theta)**2)
        E2 = self.isotropicEngConst[0] * (self.t/self.l)**3 * (self.h/self.l+sin(self.theta))/cos(self.theta)**3
        v12 = cos(self.theta)**2/((self.h/self.l+sin(self.theta))*sin(self.theta))
        v21 = ((self.h/self.l+sin(self.theta))*sin(self.theta))/cos(self.theta)**2
        G12 = self.isotropicEngConst[0] * (self.t/self.l)**3 * (self.h/self.l+sin(self.theta))/((self.h/self.l)**2*(1+2*self.h/self.l)*cos(self.theta))

        # Out-of-plane deformations
        E3 = self.isotropicEngConst[0] * (self.h/self.l+2)/(2*(self.h/self.l+sin(self.theta))*cos(self.theta)) * self.t/self.l
        v31 = self.isotropicEngConst[1]
        v13 = E1/E3 * v31
        v32 = self.isotropicEngConst[1]
        v23 = E2/E3 * v32
        G13 = (self.isotropicEngConst[0]/(2*(1+self.isotropicEngConst[1]))) * 0.577 * self.t/self.l
        G23 = (self.isotropicEngConst[0]/(2*(1+self.isotropicEngConst[1]))) * 0.577 * self.t/self.l

        self.AshbyEngConst = [E1,E2,E3,G23,G13,G12,v23,v13,v12]

        return self.AshbyEngConst
    # end: def calculate_Ashby_eng_const

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