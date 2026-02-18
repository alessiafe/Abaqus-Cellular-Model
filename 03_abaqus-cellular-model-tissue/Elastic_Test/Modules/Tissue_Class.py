from abaqus import *
from abaqusConstants import *
from caeModules import *
from part import *
from mesh import *
from visualization import *
from copy import copy
from scipy.io import loadmat
import numpy as np
import time
import imp
from re import *
import glob
from mpmath import *
mp.dps = 150 # increase precision to stabilize matrix operations (default = 15)

class Tissue_Class:

    def __init__(self):
    # Constructor of Compaction_Class
    # Args: none

        # Model name       
        self.modelName = 'Model-1'

        # Geometry parameters initialization
        self.lthick = 1. # longitudinal thickness
        self.matrixName = 'Matrix' # matrix part name
        self.fibersName = [] # names of each fiber part
        self.fibersTotName  = 'Fibers' # name of merged fibers part
        self.tissueName = 'Tissue' # name of tissue part (= matrix + fibers)
        self.fibersSetsName = [] # names of fibers set
        self.lumenSurfName = [] # names of fibers inner surfaces
        self.fibersExtSurfName = [] # names of fibers outer surfaces
        self.fibersTotExtSurfName = ['ExtSurf-Fibers', 'LeftSurf-Fibers', 'RightSurf-Fibers'] # names of merged fibers external surfaces
        self.fibersTotInnerSurfName = 'InnerSurf-Fibers' # name of fibers total inner surface
        self.matrixInnerSurfName = 'InnerSurf-Matrix' # name of matrix inner surface
        self.matrixExtSurfName = ['ExtSurf-Matrix', 'BottomSurf-Matrix',
                                  'TopSurf-Matrix', 'LeftSurf-Matrix', 'RightSurf-Matrix'] # names of matrix external surfaces
        self.matrixFaceSurfName = ['Face-0-Matrix', 'Face-1-Matrix']
        self.fibersFaceSurfName = ['Face-0-Fibers', 'Face-1-Fibers']
        self.generalContactSurfName = ['ExtSurf-Tissue', 'LateralSurf-Tissue']#'LeftSurf-Tissue', 'RightSurf-Tissue'] # names of lateral surfaces for general contact interaction
        
        # Material parameters initialization
        self.depvar = 126
        self.moistureRef = 0.12
        self.matrix_eng_const = [] # array of matrix engineering constants [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12]
        self.fiber_eng_const = [] # array of fiber engineering constants [E1, E2, E3, nu23, nu13, nu12, G23, G13, G12]
        self.matrix_density = [] # matrix density [ton/microm3]
        self.fiber_density = [] # fibers density [ton/microm3]

    	# Element features initialization
        self.meshSize = 1. # global mesh size
        self.minMeshSize = 0.1 # minimum size factor

        # Output request
        self.addRPsRequest = False # request displacement and force output (only for RPs compation control)

        return
    # end: _init_

    def totuple(self, a, scale=1):
    # Turn point list in tuple of tuples
    # Args:
    # - a = point list
    # - scale = scale factor (=1 by default)
    # Returns:
    # - a = tuple of tuples
        try:
            aL=[]
            for i in range(a[0].size):
                aL.append((a[0,i]*scale,a[1,i]*scale))
            return tuple(aL)
        except TypeError:
            return a
    # end: totuple

    def create_single_part(self, structpath, depth=1.):
    # Setup geometry of the tissue structure: create parts, sets and surfaces
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - scalefib = scale factor of fibers
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb() # start model
        myModel=self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        self.depth = depth

        ############### LOAD DATA ###############
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab structure
        self.mstruct = data['mstruct'] # get data

        ############### CREATE MATRIX ###############
        print("Creating matrix...\n")
        mask = self.totuple(self.mstruct.mask) # get mask data
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
        s.setPrimaryObject(option=STANDALONE) # set sketch objetc in viewport
        s.rectangle(point1=(min_x, min_y), point2=(max_x, max_y))
        myPart = myModel.Part(name=self.matrixName, dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
        myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude mask
        s.unsetPrimaryObject() # unset sketch objetc in viewport
        del myModel.sketches['__profile__'] # close sketch        
        
        # Create mask assembly
        myPart = myModel.parts[self.matrixName]
        myAssembly.Instance(name=self.matrixName, part=myPart, dependent=ON) # create instance
        
        # Create front and back surface of matrix
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        delta = 1             
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, -1e-5, max_x+delta, max_y+delta, 1e-5)
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, depth-1e-5, max_x+delta, max_y+delta, depth+1e-5)      
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[1]) # create set

        # Create surface of matrix
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixExtSurfName[0]) # create outer surface
           
        # Create external surface of matrix
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[0], operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]],
            myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]],))
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[0]].faces, name=self.matrixExtSurfName[0]) # create set
        # Create matrix set
        myModel.parts[self.matrixName].Set(cells=myModel.parts[self.matrixName].cells, name='Set-'+self.matrixName) # create mask set        
        myAssembly.regenerate() # regenerate assembly

        return self.modelName
    # end: create_single_part

    def setup_geometry(self, structpath, depth=1.):
    # Setup geometry of the tissue structure: create parts, sets and surfaces
    # Args:
    # - structpath = path of matlab file containing the tissue structure data
    # - scalefib = scale factor of fibers
    # Returns:
    # - modelName = name of the model

        self.mdb = Mdb() # start model
        myModel=self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly
        self.depth = depth

        ############### LOAD DATA ###############
        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab structure
        self.mstruct = data['mstruct'] # get data

        ############### CREATE FIBERS ###############
        start = time.time() # start timer
        print("Creating fibers...\n")
        fiber_out = [] # initialize outer points array
        fiber_in = [] # initialize inner points array
        for i in range(self.mstruct.fibers.size): # repeat for each fiber
            print("    {}/{}\n".format(i+1,self.mstruct.fibers.size))
            # Get fiber data
            f1 = self.mstruct.fibers[i] # get fiber data
            a_in = f1.innerpts # internal points
            a_out = f1.outerpts # external points
            fiber_out.append(a_out) # append external points
            fiber_in.append(a_in) # append inner points
            s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
            s.setPrimaryObject(option=STANDALONE) # set sketch object in viewport
            a_out = self.totuple(a_out)
            a_in = self.totuple(a_in)
            s.Spline(points=a_out) # draw external fiber
            s.Spline(points=a_in) # draf inner fiber            
            self.fibersName.append('fiber-'+str(i+1)) # append fiber name            
            myPart = myModel.Part(name=self.fibersName[-1], dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
            myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude fiber
            s.unsetPrimaryObject() # unset sketch object in viewport
            del myModel.sketches['__profile__'] # close sketch
            # Create set of fiber
            self.fibersSetsName.append('Set-'+self.fibersName[-1]) # append set name
            myPart = myModel.parts[self.fibersName[-1]]
            myPart.Set(cells=myPart.cells, name=self.fibersSetsName[-1]) # create set            
            # Create inner surface of fiber
            surface = myPart.faces.findAt(((a_in[0][0], a_in[0][1], depth/2.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_in]  # find inner face
            self.lumenSurfName.append('Lumen-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Faces=surface, name=self.lumenSurfName[-1]) # create surface
            myPart.Set(faces=surface, name=self.lumenSurfName[-1]) # create set           
            # Create outer surface of fiber
            surface = myPart.faces.findAt(((a_out[0][0], a_out[0][1], depth/2.), )) #[myPart.edges.findAt(((point[0], point[1], 0.), )) for point in a_out] # find inner face
            self.fibersExtSurfName.append('ExtFiber-'+self.fibersName[-1]) # append surface name
            myPart.Surface(side1Faces=surface, name=self.fibersExtSurfName[-1]) # create surface
            myPart.Set(faces=surface, name=self.fibersExtSurfName[-1]) # create set
            # Create fiber instances
            myAssembly.Instance(name=self.fibersName[-1], part=myPart, dependent=ON)
        
        # Merge fibers in one single part
        myAssembly.InstanceFromBooleanMerge(name=self.fibersTotName , instances=([myAssembly.instances[self.fibersName[i]]
            for i in range(len(self.fibersName))] ), keepIntersections=ON, domain=GEOMETRY, originalInstances=DELETE)
        myAssembly.features.changeKey(fromName=self.fibersTotName +'-1', toName=self.fibersTotName)
        # Create fibers assembly set 
        myModel.parts[self.fibersTotName ].Set(cells=myModel.parts[self.fibersTotName].cells, name='Set-'+self.fibersTotName) # create set
        # Merge outer surfaces of fibers
        myPart = myModel.parts[self.fibersTotName]
        myPart.SurfaceByBoolean(name=self.fibersTotExtSurfName[0], surfaces=([myPart.surfaces[self.fibersExtSurfName[i]]
            for i in range(len(self.fibersExtSurfName))]))
        # Merge inner surfaces of fibers
        myPart.SurfaceByBoolean(name=self.fibersTotInnerSurfName, surfaces=([myPart.surfaces[self.lumenSurfName[i]]
            for i in range(len(self.lumenSurfName))]))
        
        ############### CREATE MATRIX ###############
        print("Creating matrix...\n")
        mask = self.mstruct.mask # get mask data
        s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0) # open sketch
        s.setPrimaryObject(option=STANDALONE) # set sketch objetc in viewport
        mask = self.totuple(mask)
        s.Spline(points=mask) # draw external fiber
        for i in range(self.mstruct.fibers.size): # for each fiber 
            f2 = self.mstruct.fibers[i] # get fiber data
            a_in = self.totuple(f2.innerpts) # inner points
            s.Spline(points=a_in) # mask inner points
        myPart = myModel.Part(name=self.matrixName, dimensionality=THREE_D, type=DEFORMABLE_BODY) # create part
        myPart.BaseSolidExtrude(sketch=s, depth=depth) # extrude mask
        s.unsetPrimaryObject() # unset sketch objetc in viewport
        del myModel.sketches['__profile__'] # close sketch        
        
        # Create mask assembly
        myPart = myModel.parts[self.matrixName]
        myAssembly.Instance(name=self.matrixName, part=myPart, dependent=ON) # create instance
        
        # Create front and back surface of matrix
        min_x = min(point[0] for point in mask)
        max_x = max(point[0] for point in mask)
        min_y = min(point[1] for point in mask)
        max_y = max(point[1] for point in mask)
        delta = 5             
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, -1e-5, max_x+delta, max_y+delta, 1e-5)
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x-delta, min_y-delta, depth-1e-5, max_x+delta, max_y+delta, depth+1e-5)      
        myPart.Surface(side1Faces=surface, name=self.matrixFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.matrixFaceSurfName[1]) # create set

        # Create surface of matrix
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixExtSurfName[0]) # create outer surface
        
        # Create front and back surface of fibers
        myPart = myModel.parts[self.fibersTotName]
        delta = 1e-5
        surface = myPart.faces.getByBoundingBox(min_x, min_y, 0.-delta, max_x, max_y, 0.+delta)
        myPart.Surface(side1Faces=surface, name=self.fibersFaceSurfName[0]) # create back surface
        myPart.Set(faces=surface, name=self.fibersFaceSurfName[0]) # create set
        surface = myPart.faces.getByBoundingBox(min_x, min_y, depth-delta, max_x, max_y, depth+delta)
        myPart.Surface(side1Faces=surface, name=self.fibersFaceSurfName[1]) # create front surface
        myPart.Set(faces=surface, name=self.fibersFaceSurfName[1]) # create set

        # Cut fibers from mask
        tempName = 'mask-temp'
        myAssembly.InstanceFromBooleanCut(cuttingInstances=(myAssembly.instances[self.fibersTotName],), 
            instanceToBeCut=myAssembly.instances[self.matrixName], name=tempName, originalInstances=SUPPRESS)
        myAssembly.features[self.fibersTotName ].resume()
        del myModel.parts[self.matrixName]
        del myAssembly.features[self.matrixName]
        myModel.parts.changeKey(fromName=tempName, toName=self.matrixName)
        myAssembly.features.changeKey(fromName=tempName+'-1', toName=self.matrixName)     
        myAssembly.regenerate()
        
        # Create inner surface of matrix
        myPart = myModel.parts[self.matrixName]
        surface = myPart.faces
        myPart.Surface(side1Faces=surface, name=self.matrixInnerSurfName) # create surface
        myPart.SurfaceByBoolean(name=self.matrixInnerSurfName, operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixInnerSurfName], myPart.surfaces[self.matrixExtSurfName[0]], ))
        myPart.Set(faces=myPart.surfaces[self.matrixInnerSurfName].faces, name=self.matrixInnerSurfName) # create set
        
        # Create external surface of matrix
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[0], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[0]],
            myPart.surfaces[self.matrixFaceSurfName[1]],))
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[0]].faces, name=self.matrixExtSurfName[0]) # create set
        # Create matrix set
        myModel.parts[self.matrixName].Set(cells=myModel.parts[self.matrixName].cells, name='Set-'+self.matrixName) # create mask set        
        myAssembly.regenerate() # regenerate assembly

        end = time.time() # stop timer
        if round((end-start)/60.,0)<1.:
            print("Elapsed time {} seconds\n".format(int(round(end-start,0))))
        else:
            print("Elapsed time {} minutes\n".format(int(round((end-start)/60.))))

        return self.modelName
    # end: setup_geometry

    def setup_tissue_material(self, tissue):
    # Create and assign materials to composite plies.
    # Args:
    # - materialStiffMat = array of ply stiffness matrices
    # - alpha = array of ply orientations

        myModel = self.mdb.models[self.modelName]
        
        # Setup materials
        for part in [self.matrixName, self.fibersTotName]:
            myPart = myModel.parts[part]
            # Create material
            materialName = tissue if part == self.fibersTotName else 'MATRIX'
            myModel.Material(name=materialName)
            myMaterial = myModel.materials[materialName]
            myMaterial.Depvar(n=self.depvar)
            myMaterial.UserMaterial(mechanicalConstants=(0.0, ))
            # Create and assign section
            setName = 'Set-'+part
            myModel.HomogeneousSolidSection(name=materialName, material=materialName, thickness=None)
            myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
                region=myPart.sets[setName], sectionName=materialName, thicknessAssignment=FROM_SECTION)
            
        # Assign material orientation to fibers
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        for i in range(len(self.fibersName)):
            myRegion = myPart.sets[self.fibersSetsName[i]] # fiber region (set)
            primaryAxisRegion = myPart.surfaces[self.lumenSurfName[i]] # primary axis normal surface       
            myPart.MaterialOrientation(region=myRegion, 
                orientationType=DISCRETE, axis=AXIS_3, normalAxisDefinition=SURFACE, primaryAxisVector=(0.0, 0.0, 1.0), 
                normalAxisRegion=primaryAxisRegion, flipNormalDirection=False, normalAxisDirection=AXIS_3, 
                primaryAxisDefinition=VECTOR, primaryAxisDirection=AXIS_1, flipPrimaryDirection=False, 
                additionalRotationType=ROTATION_NONE, angle=0.0, additionalRotationField='', stackDirection=STACK_3)
        # Assign material orientation to matrix
        myPart = myModel.parts[self.matrixName] # matrix part
        myRegion = myPart.sets['Set-'+self.matrixName] # matrix region
        myPart.MaterialOrientation(region=myRegion, 
            orientationType=GLOBAL, axis=AXIS_1, additionalRotationType=ROTATION_NONE, 
            localCsys=None, fieldName='', stackDirection=STACK_3)
        
        
        myAssembly = myModel.rootAssembly        
        myAssembly.regenerate()

        return
    # end: def setup_tissue_material
 
    def create_single_material(self):
    # Create material
        myModel = self.mdb.models[self.modelName]
        myPart = myModel.parts[self.tissueName]
        
        # Setup material
        part = self.matrixName
        # Create material
        myModel.Material(name=part)
        myMaterial = myModel.materials[part]
        myMaterial.Depvar(n=self.depvar)
        myMaterial.UserMaterial(mechanicalConstants=(0.0, ))
        # Assign material orientation
        setName = 'Set-'+part
        # Create and assign section
        sectionName = part
        myModel.HomogeneousSolidSection(name=sectionName, material=part, thickness=None)
        myPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE,
            region=myPart.sets[setName], sectionName=sectionName, thicknessAssignment=FROM_SECTION)
        
        myAssembly = myModel.rootAssembly        
        myAssembly.regenerate()

        # Create part set
        myAssembly.SetByBoolean(name='Set-' + self.matrixName, sets=tuple(myAssembly.allInstances[self.matrixName].sets[set] for set in setName))

        return
    # end: create_single_material

    def create_mesh(self, meshSizeFiber=2., meshSizeMatrix=1., minSizeFiber=0.5, minSizeMatrix=0.9):
    # Create part mesh
    # Args:
    # - meshSize = global element size of fibers mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the specified global element size
    # - scale = meshSize scale factor to get global element size of matrix mesh (1/5 by default)

        self.meshSizeFiber = meshSizeFiber # global element size of the mesh
        self.meshSizeMatrix = meshSizeMatrix # global element size of the mesh
        self.minSizeFiber = minSizeFiber # fraction of minimum mesh size
        self.minSizeMatrix = minSizeMatrix # fraction of minimum mesh size
        #self.minMeshSize = self.meshSize * self.minSizeFactor # minimum mesh size
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        print("Create mesh...\n")
        # Element type
        elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
            secondOrderAccuracy=ON, distortionControl=DEFAULT)
        
        ## SET FIBERS MESH
        myPart = myModel.parts[self.fibersTotName] # merged fibers part
        # Set element type
        myRegion = myPart.sets['Set-'+self.fibersTotName ] # fibers region (set)
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        #myPart.setMeshControls(elemShape=HEX, technique=SWEEP, algorithm=ADVANCING_FRONT, regions=myPart.cells)
        myPart.setMeshControls(elemShape=HEX_DOMINATED, regions=myPart.cells, technique=SWEEP)      
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeFiber, size=self.meshSizeFiber)
        # Generate mesh
        myPart.generateMesh()

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=TET, technique=FREE, regions=myPart.cells)        
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSizeMatrix, size=self.meshSizeMatrix)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly

        return
    # end: create_mesh

    def create_single_mesh(self, meshSize=2., minSize=0.5):
    # Create part mesh
    # Args:
    # - meshSize = global element size of fibers mesh
    # - minSizeFactor = size of the smallest allowable element as a fraction of the specified global element size
    # - scale = meshSize scale factor to get global element size of matrix mesh (1/5 by default)

        self.meshSize = meshSize # global element size of the mesh
        self.minSize = minSize # fraction of minimum mesh size
        myModel = self.mdb.models[self.modelName] # model
        myAssembly = myModel.rootAssembly # assembly

        # Element type
        elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
        elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
            secondOrderAccuracy=ON, distortionControl=DEFAULT)

        ## SET MATRIX MESH
        myPart = myModel.parts[self.matrixName]   # matrix part
        # Set element type
        myRegion = myPart.sets['Set-'+self.matrixName]
        myPart.setElementType(elemTypes=(elemType1, elemType2, elemType3), regions=myRegion)
        # Set mesh control
        myPart.setMeshControls(elemShape=TET, technique=FREE, regions=myPart.cells)        
        # Seed part
        myPart.seedPart(deviationFactor=0.1, minSizeFactor=self.minSize, size=self.meshSize)
        # Generate mesh
        myPart.generateMesh()

        myAssembly.regenerate() # regenerate assembly
        
        return
    # end: create_single_mesh

    def find_points_within_tolerance(self, mask, tolerance):
    #    
        mask_array = np.array(mask) # convert mask to npArray
        mean_index = np.mean(mask_array[:, 0]) # calculate mean index for the x-coordinates
        close_points = [] # initialize a list to store points within the specified tolerance
        
        # Iterate through the points in the mask
        for point in mask_array:
            distance = abs(point[0] - mean_index) # calculate the distance from the mean index
            # Check if the distance is within the tolerance
            if distance <= tolerance:
                close_points.append(point) # append the point to the list if it's within tolerance
        
        return np.array(close_points)
    # end: find_points_within_tolerance

    def find_max_min_y(self, close_points):
        # Extract y-coordinates (position j) from close_points
        if close_points.size == 0:
            return None, None  # Handle the case with no close points        
        y_coordinates = close_points[:, 1]
        
        # Find maximum and minimum y values
        max_y = np.max(y_coordinates)
        min_y = np.min(y_coordinates)
        
        return max_y, min_y
    # end: find_max_min_y

    def create_surface_partition(self, meshSize, sides=False):
    #
        myModel = self.mdb.models[self.modelName] # model
        myPart = myModel.parts[self.matrixName] # matrix part
        mask = self.mstruct.mask # get mask data
        mask = self.totuple(mask)
        delta = 1.

        # Find points on surfaces
        tolerance = 25.  # Define your tolerance
        close_points = self.find_points_within_tolerance(mask, tolerance)
        ymax, ymin = self.find_max_min_y(close_points)
        xmin = min(row[0] for row in mask)
        xmax = max(row[0] for row in mask)
        
        # Partition top surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymax-meshSize/2.)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[13], faces=myPart.faces)
        myPart.Surface(name=self.matrixExtSurfName[2], side1Faces=myPart.faces.getByBoundingBox(xmin-delta, ymax-delta, 0., xmax+delta, ymax+delta))
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[2], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[2]], myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]], ))
        # Partition bottom surface
        myPart.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=ymin+meshSize/2.)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[17], faces=myPart.faces)
        myPart.Surface(name=self.matrixExtSurfName[1], side1Faces=myPart.faces.getByBoundingBox(xmin-delta, ymin-delta, 0., xmax+delta, ymin+delta))
        myPart.SurfaceByBoolean(name=self.matrixExtSurfName[1], operation=DIFFERENCE, surfaces=(
            myPart.surfaces[self.matrixExtSurfName[1]], myPart.surfaces[self.matrixFaceSurfName[0]], myPart.surfaces[self.matrixFaceSurfName[1]], ))
        
        # Partition sides surfaces
        if sides:
            surfName = 'SidesSurf-Matrix'
            myPart.SurfaceByBoolean(name=surfName, operation=DIFFERENCE, surfaces=(myPart.surfaces[self.matrixExtSurfName[0]], 
                myPart.surfaces[self.matrixExtSurfName[1]],  myPart.surfaces[self.matrixExtSurfName[2]], ))        
            faces = myPart.surfaces[surfName].faces
            xmax = max(mask, key=lambda x: x[0])[0]
            xmin = min(mask, key=lambda x: x[0])[0]
            leftSide = ()
            rightSide = ()
            for i in range(len(faces)):
                if faces[i].pointOn[0][0] < (xmin+xmax)/2.:
                    leftSide += faces[i],
                else:
                    rightSide += faces[i],
            myPart.Surface(name=self.matrixExtSurfName[4], side1Faces=FaceArray(rightSide))
            myPart.Surface(name=self.matrixExtSurfName[3], side1Faces=FaceArray(leftSide))
     
        return
    # end:create_surface_partition

    def create_model_copy(self, numCopy):
    # Create copies of the first model
    # Args:
    # - numCopy = number of the copies of the model to create
    # Returns:
    # - modelName = array of model names
        
        if numCopy==0:
            # if 0 has been passed as number of copies, warn the user and return a 1-element array with the name of the first and only model
            modelName = [self.modelName] # array of model name
            print('You passed 0 as number of model copies to create. No copy has been created.')
        else:
            newModelName = ['']*numCopy # initialize array of model names
            for k in range(numCopy): # repeat numCopy times
                newModelName[k] = 'Model-' + str(k+2) # new model name
                self.mdb.Model(name=newModelName[k], objectToCopy=self.mdb.models[self.modelName]) # copy model
            modelName = [self.modelName] + newModelName # append model name to array

        return modelName
    # end: create_model_copy

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

    def save_model(self, savepath):
    # Save .cae file
    # Args:
    # - savepath = file path to save the .cae file
    # Returns:
    # - path = full path to the file

        self.mdb.saveAs(pathName=savepath) # save file
    
        return
    # end: save_model

    def open_model(self, modelpath, structpath):
    # Open .cae file with structure and materials already created.
    # Args:
    # - savepath = folder path where the .cae file is saved
    # - filename = name of the .cae file to open
    # Returns:
    # - path = full path to the file
        
        self.mdb = openMdb(pathName=modelpath) # open file
        self.modelName = mdb.models.keys()[0]
        myModel=self.mdb.models[self.modelName]
        myPart = myModel.parts['fiber-1']
        self.b = myPart.featuresById[1].depth

        data = loadmat(structpath, squeeze_me=True, struct_as_record=False) # import matlab data structure
        self.mstruct = data['mstruct'] # get points data
        
        self.fibersName = []
        self.fibersSetsName = []
        self.lumenSurfName = []
        self.fibersExtSurfName = []
        # Initialize name arrays (setup_geometry)
        for i in range(len(myModel.parts)-2): # repeat for each fiber
            self.fibersName.append('fiber-'+str(i+1)) # append fiber name
            self.fibersSetsName.append('Set-'+self.fibersName[-1]) # append set name
            self.lumenSurfName.append('Lumen-'+self.fibersName[-1]) # append surface name
            self.fibersExtSurfName.append('ExtFiber-'+self.fibersName[-1]) # append surface name

        return self.modelName
    # end: open_model

    def set_general_constraints(self, inputs, moisture, modelName, stepTime=1, creep=False):
    # Set periodic boundary consitions (PBC)
    # Args:
    # - strains = array of strains to apply
    # - stepName = name of the step
        
        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        myInstance = myAssembly.instances[self.matrixName]
        myPart = myModel.parts[self.matrixName]

        # Create steps
        nlgeom=ON
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
                initialInc=initialInc, maxNumInc=100, minInc=1e-5, maxInc=10)        

        ## CREATE INTERACTION CONSTRAINTS
        # Tie matrix (inner) and fibers (outer)
        contIntName = 'FiberMatrix-Tie' # interaction name
        region1 = myAssembly.instances[self.fibersTotName ].surfaces[self.fibersTotExtSurfName[0]] # region 1 = fibers external surface
        region2 = myInstance.surfaces[self.matrixInnerSurfName] # region 2 = matrix inner surface
        myModel.Tie(name=contIntName, master=region1, slave=region2, constraintEnforcement=SURFACE_TO_SURFACE, 
            positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, thickness=ON) # create tie interaction
    
        myAssembly.regenerate() # regenerate assembly

        ## CREATE RPs
        # Find bounding box
        myPart.Set(nodes=myPart.nodes, name=self.matrixName+'-nodes') # create matrix nodes set
        nodeSet = myPart.sets[self.matrixName+'-nodes'] # matrix nodes set
        xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
        yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
        self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord)] # store matrix bounding box coordinates
        self.refPointName = []
        t = 10
        # Create bottom RP
        self.refPointName.append('RP-Bottom') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[2]-t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create top RP
        self.refPointName.append('RP-Top') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., self.tissueBox[3]+t, 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create left RP
        self.refPointName.append('RP-Left') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[0]-t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create right RP
        self.refPointName.append('RP-Right') # append RP name        
        myAssembly.ReferencePoint(point=(self.tissueBox[1]+t, (self.tissueBox[3]+self.tissueBox[2])/2., 0.))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create back RP
        self.refPointName.append('RP-Back') # append RP name        
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., -t))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
        # Create front RP
        self.refPointName.append('RP-Front') # append RP name
        depth = 1.     
        myAssembly.ReferencePoint(point=((self.tissueBox[1]+self.tissueBox[0])/2., (self.tissueBox[3]+self.tissueBox[2])/2., depth+t))
        myAssembly.features.changeKey(fromName=myAssembly.features.keys()[-1], toName=self.refPointName[-1])
        myAssembly.Set(name=self.refPointName[-1], referencePoints=(myAssembly.referencePoints[myAssembly.referencePoints.keys()[0]],))
               
        # Create sides node sets
        setName1 = ['Bottom-nodes', 'Top-nodes', 'Left-nodes', 'Right-nodes']
        for i in range(len(setName1)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName1[i], nodes=myPart.surfaces[self.matrixExtSurfName[1+i]].nodes)
        myPart.SetByBoolean(name=setName1[2], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[2]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        myPart.SetByBoolean(name=setName1[3], operation=DIFFERENCE, sets=(
            myPart.sets[setName1[3]], myPart.sets[setName1[0]], myPart.sets[setName1[1]],))
        # Create back/front node sets
        setName2 = ['Back-nodes', 'Front-nodes']
        for i in range(len(setName2)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName2[i], nodes=myPart.surfaces[self.matrixFaceSurfName[i]].nodes)       
        myPart.Set(name=self.matrixInnerSurfName+'-nodes', nodes=myPart.sets[self.matrixInnerSurfName].nodes)

        myPart.SetByBoolean(name=setName2[0], operation=DIFFERENCE, sets=(myPart.sets[setName2[0]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]], myPart.sets[self.matrixInnerSurfName+'-nodes'],))    
        myPart.SetByBoolean(name=setName2[1], operation=DIFFERENCE, sets=(myPart.sets[setName2[1]], myPart.sets[setName1[0]],
                            myPart.sets[setName1[1]], myPart.sets[setName1[2]], myPart.sets[setName1[3]], myPart.sets[self.matrixInnerSurfName+'-nodes'],))        
        # Merge back/front matrix nodes with fibers nodes
        myPart = myModel.parts[self.fibersTotName]
        setName2 = ['Back-nodes', 'Front-nodes']
        for i in range(len(setName2)):
            region1 = myAssembly.sets[self.refPointName[i]]
            myPart.Set(name=setName2[i], nodes=myPart.surfaces[self.fibersFaceSurfName[i]].nodes)
            myAssembly.SetByBoolean(name=setName2[i], sets=(
                myAssembly.allInstances[self.matrixName].sets[setName2[i]], 
                myAssembly.allInstances[self.fibersTotName].sets[setName2[i]], ))
            
        # Set moisture field
        fieldName = 'Moisture-field'
        for partName in [self.matrixName, self.fibersTotName]:
            setName = 'Set-'+partName
            myPart = myModel.parts[partName]
            myPart.Set(cells=myPart.cells, name=setName)
            myModel.Temperature(name=fieldName, createStepName='Initial', region=myAssembly.instances[partName].sets[setName],
                distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(moisture, ))
               
        # Apply periodic boundary conditions
        if creep:
            
            self.set_periodic_boundary_conditions(modelName, [stepName1, stepName2], inputs, loadTime)
        else:
            self.set_periodic_boundary_conditions(modelName, [stepName1], inputs, loadTime)

        return
    # end: set_general_constraints

    def set_periodic_boundary_conditions(self, modelName, steps, load, loadTime):

        myModel = self.mdb.models[modelName]
        myAssembly = myModel.rootAssembly
        stepName = steps[0]
        
        # Initialize load matrix and face mapping
        loads = [[0.0]*6 for _ in range(6)]
        face_map = {0:(2,3), 1:(0,1), 2:(4,5)}
        # Distribute normal loads
        for axis in (0,1,2):
            if load[axis] != 0.0:
                neg_face, pos_face = face_map[axis]
                loads[neg_face][axis] = -load[axis]
                loads[pos_face][axis] =  load[axis]
        # Distribute shear loads
        shear_map = {
            3: {'facesA': (0,1), 'axisA':2, 'facesB':(4,5), 'axisB':1},
            4: {'facesA': (2,3), 'axisA':2, 'facesB':(4,5), 'axisB':0},
            5: {'facesA': (0,1), 'axisA':0, 'facesB':(2,3), 'axisB':1},}
        for sid, props in shear_map.items():
            if load[sid] != 0.0:
                for idx, face in enumerate(props['facesA']):
                    sign = -1.0 if idx == 0 else 1.0
                    loads[face][props['axisA']] = sign * load[sid]
                for idx, face in enumerate(props['facesB']):
                    sign = -1.0 if idx == 0 else 1.0
                    loads[face][props['axisB']] = sign * load[sid]
        # Identify normal and shear axes
        normal_axes = [i for i in (0,1,2)   if load[i] != 0.0]
        shear_axes  = [i for i in (3,4,5)   if load[i] != 0.0]

        # Create face sets
        myPart = myModel.parts[self.matrixName]
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[2]].faces, name=self.matrixExtSurfName[2])
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[1]].faces, name=self.matrixExtSurfName[1])
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[3]].faces, name=self.matrixExtSurfName[3])
        myPart.Set(faces=myPart.surfaces[self.matrixExtSurfName[4]].faces, name=self.matrixExtSurfName[4])

        # Set names
        const_names = ['Bottom-Coupling','Top-Coupling','Left-Coupling',
                       'Right-Coupling','Back-Coupling','Front-Coupling']
        set_names = ['Bottom-nodes','Top-nodes','Left-nodes',
                    'Right-nodes','Back-nodes','Front-nodes']
        
        # Create amplitude for forces
        amp_name = 'Amp-Force'
        myModel.EquallySpacedAmplitude(name=amp_name, timeSpan=STEP,
            smooth=SOLVER_DEFAULT, fixedInterval=loadTime,
            begin=0.0, data=(0.0, 1.0))

        ########### Axial cases ###########
        if len(normal_axes) >= 1 and not shear_axes:
            # Set X-symmetry on left face
            left_set = myAssembly.instances[self.matrixName].sets[self.matrixExtSurfName[3]]
            myModel.XsymmBC(name='SYM-X-left', createStepName='Initial',
                            region=left_set, localCsys=None)
            # Set Y-symmetry on bottom face
            bottom_set = myAssembly.instances[self.matrixName].sets[self.matrixExtSurfName[1]]
            myModel.YsymmBC(name='SYM-Y-bottom', createStepName='Initial',
                            region=bottom_set, localCsys=None)
            # Merge back faces and apply Z-symmetry
            myAssembly.SetByBoolean(name='Face-0', sets=(
                myAssembly.instances[self.matrixName].sets[self.matrixFaceSurfName[0]],
                myAssembly.instances[self.fibersTotName].sets[self.fibersFaceSurfName[0]],))
            back_region = myAssembly.sets['Face-0']
            myModel.ZsymmBC(name='SYM-Z-back', createStepName='Initial',
                            region=back_region, localCsys=None)

            # For each normal axis, apply force on positive face and pin rotations
            for axis in normal_axes:
                # Map faces
                neg_face, pos_face = face_map[axis]
                ref_name = self.refPointName[pos_face]
                # Read loads
                fx, fy, fz = loads[pos_face][0:3]
                # Couple ref. point to nodes set (only loaded directions)               
                coupx = ON if fx != 0. else OFF
                coupy = ON if fy != 0. else OFF
                coupz = ON if fz != 0. else OFF
                couprx = OFF if fx != 0. else ON
                coupry = OFF if fy != 0. else ON
                couprz = OFF if fz != 0. else ON
                region_cp = myAssembly.sets[ref_name]
                if fz != 0.:
                    region_ns = myAssembly.sets[set_names[pos_face]]
                else:
                    region_ns = myAssembly.instances[self.matrixName].sets[set_names[pos_face]]
                myModel.Coupling(name=const_names[pos_face], controlPoint=region_cp, 
                    surface=region_ns, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,                   
                    localCsys=None, u1=coupx, u2=coupy, u3=coupz, ur1=couprx, ur2=coupry, ur3=couprz)
                # Apply concentrated force
                myModel.ConcentratedForce(name='CF-%s' % ref_name,
                    createStepName=stepName, amplitude=amp_name,
                    region=myAssembly.sets[ref_name], cf1=fx, cf2=fy, cf3=fz)
                # Pin rotations
                ux = UNSET if fx != 0. else 0.
                uy = UNSET if fy != 0. else 0.
                uz = UNSET if fz != 0. else 0.
                myModel.DisplacementBC(name='BC-%s' % ref_name,
                    createStepName=stepName, region=myAssembly.sets[ref_name],
                    u1=ux, u2=uy, u3=uz, ur1=0.0, ur2=0.0, ur3=0.0,
                    amplitude=UNSET, distributionType=UNIFORM,
                    fieldName='', fixed=OFF, localCsys=None)
                # Create history output
                for step_name in steps:
                    myModel.HistoryOutputRequest(name='HO-%s' % ref_name,
                        createStepName=step_name, variables=('U1','U2','U3','CF1','CF2','CF3'),
                        region=myAssembly.sets[ref_name], sectionPoints=DEFAULT, rebar=EXCLUDE)
        
        ########### Shear cases ###########
        elif shear_axes and len(normal_axes) == 0:
            # Select active face
            active_faces = [i for i,row in enumerate(loads) if any(v != 0.0 for v in row)]
            apply_xsymm = not any(f in active_faces for f in (2,3))
            apply_ysymm = not any(f in active_faces for f in (0,1))
            apply_zsymm = not any(f in active_faces for f in (4,5))
            # Set symmetry
            if apply_xsymm:
                # Set X-symmetry on left face
                left_set = myAssembly.instances[self.matrixName].sets[self.matrixExtSurfName[3]]
                myModel.XsymmBC(name='SYM-X-left', createStepName='Initial',
                                region=left_set, localCsys=None)
            if apply_ysymm:
                # Set Y-symmetry on bottom face
                bottom_set = myAssembly.instances[self.matrixName].sets[self.matrixExtSurfName[1]]
                myModel.YsymmBC(name='SYM-Y-bottom', createStepName='Initial',
                                region=bottom_set, localCsys=None)
            if apply_zsymm:
                # Merge back faces and apply Z-symmetry
                myAssembly.SetByBoolean(name='Face-0', sets=(
                    myAssembly.instances[self.matrixName].sets[self.matrixFaceSurfName[0]],
                    myAssembly.instances[self.fibersTotName].sets[self.fibersFaceSurfName[0]],))
                back_region = myAssembly.sets['Face-0']
                myModel.ZsymmBC(name='SYM-Z-back', createStepName='Initial',
                                region=back_region, localCsys=None)

            # Couple ref. point to nodes set (only loaded directions)
            for i in active_faces[2:]:
                # Read loads
                fx, fy, fz = loads[i][0:3]
                region1 = myAssembly.sets[self.refPointName[i]]
                region2 = (myAssembly.instances[self.matrixName].sets[set_names[i]]
                        if i < 4 else myAssembly.sets[set_names[i]])
                coupx = ON
                coupy = ON
                coupz = ON
                couprx = ON
                coupry = ON
                couprz = ON
                myModel.Coupling(name=const_names[i], controlPoint=region1, surface=region2, 
                    influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None,
                    u1=coupx, u2=coupy, u3=coupz, ur1=couprx, ur2=coupry, ur3=couprz)
                # Apply concentrated force
                ref_name = self.refPointName[i]
                myModel.ConcentratedForce(name='CF-%s' % ref_name,
                    createStepName=stepName, amplitude=amp_name,
                    region=myAssembly.sets[ref_name],
                    cf1=fx, cf2=fy, cf3=fz)
                # Pin displacements and rotations
                ux = UNSET if fx != 0. else 0.
                uy = UNSET if fy != 0. else 0.
                uz = UNSET if fz != 0. else 0.
                myModel.DisplacementBC(name='BC-%s' % ref_name,
                    createStepName=stepName, region=myAssembly.sets[ref_name],
                    u1=ux, u2=uy, u3=uz, ur1=0., ur2=0., ur3=0.,
                    amplitude=UNSET, distributionType=UNIFORM,
                    fieldName='', fixed=OFF, localCsys=None)
                # Create history output
                for step_name in steps:
                    myModel.HistoryOutputRequest(name='HO-%s' % ref_name,
                        createStepName=step_name, variables=('U1','U2','U3','CF1','CF2','CF3'),
                        region=myAssembly.sets[ref_name], sectionPoints=DEFAULT, rebar=EXCLUDE)
                    
            for i in active_faces[:2]:
                # Read loads
                fx, fy, fz = loads[i][0:3]
                region1 = myAssembly.sets[self.refPointName[i]]
                region2 = (myAssembly.instances[self.matrixName].sets[set_names[i]]
                        if i < 4 else myAssembly.sets[set_names[i]])
                coupx = ON if (not any(f in active_faces for f in (2,3)) or fx != 0.) else OFF
                coupy = ON if (not any(f in active_faces for f in (0,1)) or fy != 0.) else OFF
                coupz = ON if (not any(f in active_faces for f in (4,5)) or fz != 0.) else OFF
                couprx = ON if any(f in active_faces for f in (2,3)) else OFF
                coupry = ON if any(f in active_faces for f in (0,1)) else OFF
                couprz = ON if any(f in active_faces for f in (4,5)) else OFF
                myModel.Coupling(name=const_names[i], controlPoint=region1, surface=region2, 
                    influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None,
                    u1=coupx, u2=coupy, u3=coupz, ur1=couprx, ur2=coupry, ur3=couprz)
                # Pin displacements and rotations
                ux = 0.
                uy = 0.
                uz = 0.
                ref_name = self.refPointName[i]
                myModel.DisplacementBC(name='BC-%s' % ref_name,
                    createStepName=stepName, region=myAssembly.sets[ref_name],
                    u1=ux, u2=uy, u3=uz, ur1=0., ur2=0., ur3=0.,
                    amplitude=UNSET, distributionType=UNIFORM,
                    fieldName='', fixed=OFF, localCsys=None)
                # Create history output
                for step_name in steps:
                    myModel.HistoryOutputRequest(name='HO-%s' % ref_name,
                        createStepName=step_name, variables=('U1','U2','U3','CF1','CF2','CF3'),
                        region=myAssembly.sets[ref_name], sectionPoints=DEFAULT, rebar=EXCLUDE)
            
            #print(steps)
            # FIX NODE TO AVOID FREE SHIFTING
            delta = 5
            bcName = 'Fix-Node'
            x_coord = (self.tissueBox[1]+self.tissueBox[0])/2. #253
            y_coord = (self.tissueBox[3]+self.tissueBox[2])/2. #211
            # Get matrix nodes
            myInstance = myAssembly.instances[self.matrixName]
            nodes1 = myInstance.nodes.getByBoundingBox(x_coord-delta,
                y_coord-delta, 0.2*self.b, x_coord+delta, y_coord+delta, 0.8*self.b)
            # Get fibers nodes
            myInstance = myAssembly.instances[self.fibersTotName]
            nodes2 = myInstance.nodes.getByBoundingBox(x_coord-delta,
                y_coord-delta, 0.2*self.b, x_coord+delta, y_coord+delta, 0.8*self.b)
            # Apply encastre
            myModel.EncastreBC(createStepName=stepName, localCsys=None, name=bcName,
                region=myAssembly.Set(name=bcName,nodes=nodes1+nodes2))
            
        #for step_name in steps:
        #myModel.FieldOutputRequest(name='STATEV', createStepName=stepName, variables=('SDV',))
            
        myAssembly.regenerate()

        return
    # end: set_periodic_boundary_conditions

    def calculate_stress_strain(self, jobName, inc=-1):
    # Calculate stress and strain from the reaction forces and displacements of the reference points.
    # Args:
    # - jobName = name of the job
    # - stepName = name of the step
    # - inc = number of the job increment
    # Returns:
    # - stress = array of stress components
    # - strain = array of strain components

        # Open odb file
        odb = openOdb(jobName + ".odb")
        stepName = odb.steps.keys()[-1]
        myStep = odb.steps[stepName]
        myAssembly = odb.rootAssembly
        numComp = 6
     
        # Determine how to handle `inc`
        if isinstance(inc, int):  
            inc = [inc]  # convert single number to list
            multi_inc = False
        elif isinstance(inc, list) and len(inc) == 0:  
            # Auto-detect available increments
            ref_node = 'Assembly ASSEMBLY'
            inc = [i[0] for i in myStep.historyRegions[ref_node].historyOutputs['ALLAE'].data]
            multi_inc = True

        # Initialize arrays     
        strains_all = []
        stresses_all = []
        for k in range(len(inc)):
            if len(inc) == 1: k = inc[k] 
            # Initialize arrays
            RF = []
            U = []
            strain = [0.]*numComp
            stress = [0.]*numComp
            
            
            # Store displacements and forces
            for i in range(numComp):
                myRegion = 'Node ASSEMBLY.{}'.format(i+1)
                if myRegion not in myStep.historyRegions:
                    RF.append([0, 0, 0])
                    U.append([0, 0, 0])
                else:
                    myHistory = myStep.historyRegions[myRegion]
                    RF.append([myHistory.historyOutputs['CF1'].data[k][1], myHistory.historyOutputs['CF2'].data[k][1], myHistory.historyOutputs['CF3'].data[k][1]])
                    U.append([myHistory.historyOutputs['U1'].data[k][1], myHistory.historyOutputs['U2'].data[k][1], myHistory.historyOutputs['U3'].data[k][1]])

            
            # Find bounding box
            nodeSet = myAssembly.instances['MATRIX'].nodeSets['SET-MATRIX']
            xCoord = [curNode.coordinates[0] for curNode in nodeSet.nodes]
            yCoord = [curNode.coordinates[1] for curNode in nodeSet.nodes]
            zCoord = [curNode.coordinates[2] for curNode in nodeSet.nodes]
            self.tissueBox = [min(xCoord), max(xCoord), min(yCoord), max(yCoord), min(zCoord), max(zCoord)] # store matrix bounding box coordinates       
            self.LatticeVector = [[self.tissueBox[1]-self.tissueBox[0], 0., 0.], [0., self.tissueBox[3]-self.tissueBox[2], 0.], [0.0, 0.0, self.tissueBox[5]-self.tissueBox[4]]]
            
            # Calculate stresses from reaction forces
            stress[0] = (RF[2][0]+RF[3][0])/(self.LatticeVector[1][1]*self.LatticeVector[2][2])
            stress[1] = (RF[0][1]+RF[1][1])/(self.LatticeVector[0][0]*self.LatticeVector[2][2])
            stress[2] = (RF[4][2]+RF[5][2])/(self.LatticeVector[0][0]*self.LatticeVector[1][1])
            stress[3] = ((abs(RF[0][2])+abs(RF[1][2]))/(self.LatticeVector[0][0]*self.LatticeVector[2][2])+
                        (abs(RF[4][1])+abs(RF[5][1]))/(self.LatticeVector[1][1]*self.LatticeVector[0][0]))/2. # yz       
            stress[4] = ((abs(RF[2][2])+abs(RF[3][2]))/(self.LatticeVector[2][2]*self.LatticeVector[1][1])+
                        (abs(RF[4][0])+abs(RF[5][0]))/(self.LatticeVector[0][0]*self.LatticeVector[1][1]))/2. # xz       
            stress[5] = ((abs(RF[0][0])+abs(RF[1][0]))/(self.LatticeVector[0][0]*self.LatticeVector[2][2])+
                        (abs(RF[2][1])+abs(RF[3][1]))/(self.LatticeVector[1][1]*self.LatticeVector[2][2]))/2. # xy

            # Calculate strains from displacements
            strain[0] = (U[2][0]+U[3][0])/self.LatticeVector[0][0]
            strain[1] = (U[0][1]+U[1][1])/self.LatticeVector[1][1]
            strain[2] = (U[4][2]+U[5][2])/self.LatticeVector[2][2]

            strain[3] = ((abs(U[0][2])+abs(U[1][2]))/self.LatticeVector[1][1]+
                        (abs(U[4][1])+abs(U[5][1]))/self.LatticeVector[2][2]) # engineering yz            
            strain[4] = ((abs(U[2][2])+abs(U[3][2]))/self.LatticeVector[0][0]+
                        (abs(U[4][0])+abs(U[5][0]))/self.LatticeVector[2][2]) # engineering xz
            strain[5] = ((abs(U[0][0])+abs(U[1][0]))/self.LatticeVector[1][1]+
                        (abs(U[2][1])+abs(U[3][1]))/self.LatticeVector[0][0]) # engineering xy

        # Store results
            strains_all.append(strain)
            stresses_all.append(stress)
        odb.close()

        # Return a single array if there was only one increment, otherwise return lists
        if not multi_inc:
            return stresses_all[0], strains_all[0], inc[0]
        else:
            return stresses_all, strains_all, inc
    # end: def calculate_stress

    def compute_3D_compliance_matrix(self, value, moisture, run, umatPath, stepTime=1, creepActivate=False):
    # Compute 3D stiffness matrix by running 9 elementary cases.

        # Set elementary cases strains
        val_x = [value[0], 0., 0., 0., 0., 0., value[0], value[0], 0.]
        val_y = [0., value[1], 0., 0., 0., 0., value[1], 0., value[1]]
        val_z = [0., 0., value[2], 0., 0., 0., 0., value[2], value[2]]
        val_yz = [0., 0., 0., value[3], 0., 0., 0., 0., 0.]
        val_xz = [0., 0., 0., 0., value[4], 0., 0., 0., 0.]
        val_xy = [0., 0., 0., 0., 0., value[5], 0., 0., 0.]

        # Activate/deactivate creep in umat
        parent_path = os.path.dirname(umatPath)
        filepath = os.path.join(parent_path, "activate_ve_creep.inc")
        self.update_umat_inc_file(creepActivate, filepath)

        if run:
            # Copy model
            numElCases = 9
            modelName = self.create_model_copy(numElCases-1)
        
            # Run elementary cases
            for k in range(numElCases):
                modelName = 'Model-{}'.format(k+1)
                inputs = [val_x[k], val_y[k], val_z[k], val_yz[k], val_xz[k], val_xy[k]]
                # Set boundary conditions
                self.set_general_constraints(inputs, moisture, modelName, stepTime, creepActivate)
                # Run job
                jobName = "Job-{}".format(k+1)
                self.submit_job(jobName, modelName, umatPath)
        else:
            # Build model to store variables and names
            k = 0
            modelName = 'Model-{}'.format(k+1)
            inputs = [val_x[k], val_y[k], val_z[k], val_yz[k], val_xz[k], val_xy[k]]
            # Set boundary conditions
            self.set_general_constraints(inputs, moisture, modelName, stepTime, creepActivate)

        # Calculate stiffness matrix
        jobName = 'Job'
        inc = [] if creepActivate else -1
        comp_coeff = self.calculate_compliance_matrix(jobName, inc)

        return comp_coeff
    # end: def compute_3D_compliance_matrix

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
            stress, strain, inc_val = self.calculate_stress_strain(jobName=jobName + "-{}".format(i+1), inc=inc)
            stress_mat.append(stress)
            strain_mat.append(strain)
            inc_mat.append(inc_val)
            print('Case {}'.format(i+1))
            print('Stress: ', stress)
            print('Strain: ',strain)
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

        else: # inc==[] so that each elementary case returns arrays (one per inc)
            
            # Diagonal elements
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
            'time11': inc if isinstance(inc, int) else inc_mat[0][1:],
            'D11': self.CompMat[0][0],
            'time22': inc if isinstance(inc, int) else inc_mat[1][1:],
            'D22': self.CompMat[1][1],
            'time33': inc if isinstance(inc, int) else inc_mat[2][1:],
            'D33': self.CompMat[2][2],
            'time44': inc if isinstance(inc, int) else inc_mat[3][1:],
            'D44': self.CompMat[3][3],
            'time55': inc if isinstance(inc, int) else inc_mat[4][1:],
            'D55': self.CompMat[4][4],
            'time66': inc if isinstance(inc, int) else inc_mat[5][1:],
            'D66': self.CompMat[5][5],
            'time12': inc if isinstance(inc, int) else inc_mat[6][1:],
            'D12': self.CompMat[0][1],
            'time13': inc if isinstance(inc, int) else inc_mat[7][1:],
            'D13': self.CompMat[0][2],
            'time23': inc if isinstance(inc, int) else inc_mat[8][1:],
            'D23': self.CompMat[1][2]
        }
        return result
    # end: def calculate_stiffness_matrix

    def calculate_eng_const(self):
    # Calculate the engineering constants of the RVE from simulations results.
    # Args: none
    # Return:
    # - RVEEngConst = array of engineering constants calculated from simulations results
    #                 [E1,E2,E3,nu23,nu13,nu12,G23,G13,G12]

        self.complianceMatrix = self.CompMat
        
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

        return self.RVEEngConst
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