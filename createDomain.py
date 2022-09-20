# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
import csv

def getCurrentDisplayedObjectAndModel():
    vp = session.viewports[session.currentViewportName]
    p = vp.displayedObject # possibly a part
    if p is not None:
        model = mdb.models[p.modelName]
    else:
        # if there is no displayed object, then return the first model
        model = mdb.models.values()[0]
    return p, model

def createPreliminaryModel(lengths, DRM_depth, PML_depth, partName, materialName, 
    youngsModulus, poissonsRatio, density, timeIncrement, duration, jobName, 
    noDamping=False, applyRampForce=False, stepApplication=None):
    # Creating the part
    p, model = getCurrentDisplayedObjectAndModel()
    s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.rectangle(point1=(0.0, 0.0), point2=(lengths['x'], lengths['y']))
    p = model.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=lengths['z'])
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=p)
    del model.sketches['__profile__']

    # Creating datum planes for partitioning
    principalPlanes = {'x': YZPLANE, 'y': XZPLANE, 'z': XYPLANE}
    for direction, length in lengths.items():
        if direction != 'z':
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth)
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth+DRM_depth)
        p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=length-PML_depth)
        p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=length-PML_depth-DRM_depth)
        
    # Partitioning
    for plane in p.datums.values():
        p.PartitionCellByDatumPlane(datumPlane=plane, cells=p.cells)

    # Set up corresponding cell array
    cells = {}
    cells['internal'] = p.cells.getByBoundingBox(xMin=PML_depth+DRM_depth, yMin=PML_depth+DRM_depth, zMax=lengths['z']-PML_depth-DRM_depth,
        xMax=lengths['x']-PML_depth-DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
    cells['DRM'] = p.cells.getByBoundingBox(xMin=PML_depth, yMin=PML_depth, zMax=lengths['z']-PML_depth,
        xMax=lengths['x']-PML_depth, yMax=lengths['y']-PML_depth)
    cells['DRM'] = part.CellArray([cell for cell in cells['DRM'] if cell not in cells['internal']])
    cells['PML'] = part.CellArray([cell for cell in p.cells if cell not in cells['internal'] and cell not in cells['DRM']])
    # Create sets
    sets = {cellName: p.Set(cells=cells[cellName], name=cellName) for cellName in cells.keys()}
    
    # Create material properties
    material = model.Material(name='soil')
    material.Density(table=((density, ), ))
    material.Elastic(table=((youngsModulus, poissonsRatio), ))
    # Create section and section assignments
    section = model.HomogeneousSolidSection(name='soil', material='soil', thickness=None)
    [p.SectionAssignment(region=cellSet, sectionName='soil') for cellName, cellSet in sets.items() if cellName != 'PML'] 

    # Create the instance
    a = model.rootAssembly
    domainInstance = a.Instance(name=p.name+'-1', part=p, dependent=ON)

    # Create the BC
    externalFaces = p.faces.getByBoundingBox(zMin=lengths['z']) + \
        p.faces.getByBoundingBox(xMax=0) + p.faces.getByBoundingBox(xMin=lengths['x']) + \
        p.faces.getByBoundingBox(yMax=0) + p.faces.getByBoundingBox(yMin=lengths['y'])
    externalFacesSet = p.Set(faces=externalFaces, name='externalFaces')
    model.DisplacementBC(name='fixExteralFaces', createStepName='Initial', 
        region=domainInstance.sets['externalFaces'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
    
    # Mesh
    p.seedPart(size=DRM_depth, deviationFactor=0.1, minSizeFactor=0.1)
    # NOTE: not sure whether it's necessary to set up elemType2 and elemType3 for Wedge and Tet elements
    # NOTE 2: For wave propagation problems, reduced integration elements are generally not used. (There are no bending problem.)
    elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
    p.setElementType(regions=(p.cells, ), elemTypes=(elemType1, ))
    # NOTE: This is dummy element 
    dummyElemType = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
    p.setElementType(regions=p.sets['PML'], elemTypes=(dummyElemType, ))
    # NOTE: Generate mesh set by set can make the node numbering on each set continuous
    [p.generateMesh(regions=cellSet.cells) for cellSet in sets.values()]

    # Create element/node sets
    eleDRM = p.Set(elements=p.sets['DRM'].elements, name='eleDRM')
    nodes_on_DRM_external_sides = p.sets['DRM'].nodes.getByBoundingBox(zMin=lengths['z']-PML_depth) + \
        p.sets['DRM'].nodes.getByBoundingBox(xMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(xMin=lengths['x']-PML_depth) + \
        p.sets['DRM'].nodes.getByBoundingBox(yMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(yMin=lengths['y']-PML_depth)
    sideDRM = p.Set(nodes=nodes_on_DRM_external_sides, name='sideDRM')
    nodes_on_DRM_internal_sides = p.sets['DRM'].nodes.getByBoundingBox(zMax=lengths['z']-PML_depth-DRM_depth,
        xMin=PML_depth+DRM_depth, xMax=lengths['x']-PML_depth-DRM_depth,
        yMin=PML_depth+DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
    inDRM = p.Set(nodes=nodes_on_DRM_internal_sides, name='inDRM')
    # Prepare sets for history output
    node_at_center = p.nodes.getClosest((lengths['x']/2, lengths['y']/2, lengths['z']/2))
    # NOTE: Abaqus Set only takes MeshNodeArray instead of MeshNode.
    centerNode = p.Set(nodes=mesh.MeshNodeArray((node_at_center,)), name='centerNode')
    node_at_top_center = p.nodes.getClosest((lengths['x']/2, lengths['y']/2, 0))
    topCenterNode = p.Set(nodes=mesh.MeshNodeArray((node_at_top_center,)), name='topCenterNode')

    # Create step
    step1 = model.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=duration, 
        maxNumInc=10000, initialInc=timeIncrement, maxInc=timeIncrement, minInc=duration/100000.0)
    if stepApplication is not None:
        if stepApplication.lower() == 'transient fidelity':
            step1.setValues(application=TRANSIENT_FIDELITY, initialConditions=ON)
        elif stepApplication.lower() == 'moderate dissipation':
            step1.setValues(application=MODERATE_DISSIPATION, initialConditions=OFF)
        elif stepApplication.lower() == 'quasi static':
            step1.setValues(application=QUASI_STATIC, amplitude=RAMP)
        else:
            raise ValueError("The step application is not valid.")
    # else: # NOTE: This is not necessary. Just for reference.
    #     step1.setValues(application=DEFAULT, amplitude=STEP, initialConditions=DEFAULT)
    if noDamping:
        step1.setValues(alpha=0.0) # NOTE: the default value is calculated automatically
    if applyRampForce:
        step1.setValues(amplitude=RAMP) # NOTE: the default value is STEP
    # model.ExplicitDynamicsStep(name='Step-1', previous='Initial', 
    #     timePeriod=duration, improvedDtMethod=ON)
    # NOTE: Change the output request frequency. Be aware if the model is not using a default output request setting.
    for fieldOutputRequest in model.fieldOutputRequests.values():
        fieldOutputRequest.setValues(numIntervals=50, timeMarks=OFF)
    for historyOutputRequest in model.historyOutputRequests.values():
        historyOutputRequest.setValues(frequency=1)
    # Create history outputs
    HOutputRequests = {'H-Output-Center': 'centerNode', 'H-Output-TopCenter': 'topCenterNode'}
    for HOutputName, nodeSetName in HOutputRequests.items():
        model.HistoryOutputRequest(name=HOutputName, createStepName=step1.name, 
            variables=('U1', 'U2', 'U3', 'V1', 'V2', 'V3', 'A1', 'A2', 'A3'), frequency=1, 
            region=domainInstance.sets[nodeSetName])

    # Create job and write the input file
    job = mdb.Job(name=jobName+'_pre', model=model.name)
    job.writeInput()
    return


# def modifyInputAndCheck(lengths, PML_depth, partName, materialName, alpha, beta, jobName, subroutineFileName):
#     modifyInput(lengths, PML_depth, partName, materialName, alpha, beta, jobName)
#     # Submit a new job from the created input file with the subroutine
#     mdb.JobFromInputFile(name=jobName, inputFileName=jobName+'.inp', 
#         userSubroutine=subroutineFileName)
#     # mdb.jobs[jobName].submit(datacheckJob=True) # Check the input. Not sure whether it's necessary.
#     return

# def modifyInput(lengths, PML_depth, partName, materialName, alpha, beta, jobName, cLoadFileName='Cload.txt'):
#     density, youngsModulus, poissonsRatio = getMaterialPropertiesFromInputFile(jobName, materialName)
#     # Read the input file and make modifications for UEL use
#     with open(jobName+'_pre.inp', 'r') as f:
#         lines = f.readlines()
#         dummyElementLine = getLineIndex(lines, '*Element, type=C3D8R\n', isConversionNeeded=True)
#         lines[dummyElementLine] = '*Element, type=U3\n'
#         # [NOTE] Info for parameters under *USER ELEMENT:
#         #   TYPE: Must be 'Un' where n is a positive integer less than 10000, and it must be the same as the element type key used to identify this element on the *ELEMENT option.
#         #   NODE: Number of nodes associated with an element of this type.
#         #   COORDINATES: 3 for 3D, 2 for 2D.
#         #   PROPERTY: The number of property values needed as data in UEL to define such an element.
#         #   VARIABLES: 360 for 3D problems and 80 for 2D problems. This is determined by Wenyang's UEL subroutine.
#         #   [TODO (maybe?)] So how to calculate the number of or variables (VARIABLES) that are needed for *USER ELEMENT?
#         #   Numbers in the second line: The active DoF
#         lines.insert(dummyElementLine, '*User Element, Type=U3, Nodes=8, Coordinates=3, Properties=12, Variables=360, Unsymm\n')
#         lines.insert(dummyElementLine+1, '1, 2, 3, 21, 22, 23, 24, 25, 26\n')
#         # NOTE: Since lines are changed above, we cannot create lowerLines like in other functions and reuse it. As a result, we need to create the lowerLines list again (with isConversionNeeded=True).
#         endPartLine = getLineIndex(lines, '*End Part\n', isConversionNeeded=True)
#         # [NOTE] Info for parameters used for UEL
#         #   E: Young's Modulus
#         #   xnu: Poisson Ratio
#         #   rho: Density
#         #   EleType_pos: 0.0 (A trivial property)
#         #   PML_L: Depth of PML (distance from the boundary of DRM layer to the exterior of PML region)
#         #   afp: Described as `m` and is equal to 2.0 in Wenyang's paper. It is the polynomial degree for alpha and beta functions.
#         #   PML_Rcoef: Described as `R` and is equal to 10^(-10) in Wenyang's paper. It is a user-tunable reflection coefficient.
#         #   RD_half_width_x: Half width of the domain in x direction without PML layers (i.e., the distance in x direction from the center to the boundary of DRM layer)
#         #   RD_half_width_y: Half width of the domain in y direction without PML layers (i.e., the distance in y direction from the center to the boundary of DRM layer)
#         #   RD_depth: The depth (in z direction) of the domain without PML (i.e., the depth of the interested region + the thickness of DRM layer)
#         #   Damp_alpha, Damp_beta: alpha and beta used in Rayleigh Damping
#         lines[endPartLine:endPartLine] = ['*Parameter\n', 'E=%e\n'%youngsModulus, 'xnu=%f\n'%poissonsRatio, 'rho=%f\n'%density, 
#             'EleType_pos=0.0\n', 'PML_L=%f\n'%PML_depth, 'afp=2.0\n', 'PML_Rcoef=1e-10\n', 
#             'RD_half_width_x=%f\n'%(lengths['x']/2-PML_depth), 'RD_half_width_y=%f\n'%(lengths['y']/2-PML_depth),
#             'RD_depth=%f\n'%(lengths['z']-PML_depth), 'Damp_alpha=%f\n'%alpha, 'Damp_beta=%f\n'%beta, 
#             '*UEL Property, elset=PML\n', '<E>, <xnu>, <rho>, <EleType_pos>, <PML_L>, <afp>, <PML_Rcoef>, <RD_half_width_x>,\n', 
#             '<RD_half_width_y>, <RD_depth>, <Damp_alpha>, <Damp_beta>\n']
#         # NOTE: Again, the lines are changed above, so we need to create the lowerLines list again.
#         endStepLine = getLineIndex(lines, '*End Step\n', isConversionNeeded=True)
#         lines.insert(endStepLine, '*Include, input=%s\n'%cLoadFileName)
#     with open(jobName+'.inp', 'w') as f:
#         f.writelines(lines)
#     return

def getLineIndex(lines, target, startLine=0, isConversionNeeded=False):
    if isConversionNeeded: # NOTE: make the search case insensitive due to different input file writing preferences
        lines = [line.lower() for line in lines]
    return startLine + lines[startLine:].index(target.lower())

def getMaterialPropertiesFromInputFile(jobName, materialName):
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    materialLine = getLineIndex(lowerLines, '*Material, name=%s\n'%materialName)
    densityLine = getLineIndex(lowerLines, '*Density\n', startLine=materialLine)
    density = float(lines[densityLine+1].split(',')[0])
    elasticLine = getLineIndex(lowerLines, '*Elastic\n', startLine=materialLine)
    elasticProperties = lines[elasticLine+1].split(',')
    youngsModulus = float(elasticProperties[0])
    poissonsRatio = float(elasticProperties[1])
    return density, youngsModulus, poissonsRatio

def setViews():
    arguments = ['nearPlane', 'farPlane', 'width', 'height', 'projection', 
        'cameraPosition', 'cameraUpVector', 'cameraTarget', 
        'viewOffsetX', 'viewOffsetY']
    # NOTE: this function supports up to 4 view vectors
    viewVectors = [(1, 1, -1), (-1, -1, -1)]
    viewport = session.viewports[session.currentViewportName]
    for i, viewVector in enumerate(viewVectors):
        viewport.view.setViewpoint(viewVector=viewVector, cameraUpVector=(0, 0, -1))
        session.View('User-%d'%(i+1), *[getattr(viewport.view, arg) for arg in arguments], autoFit=ON)

def readDataForCreatingAbaqusModel(fileName='dataForAbaqusModel.csv'):
    with open(fileName, 'r') as f:
        reader = csv.DictReader(f)
        for line in reader:
            data = {}
            for key in line.keys():
                try:
                    data[key] = float(line[key])
                except ValueError:
                    data[key] = line[key]
            return data

if __name__ == '__main__':
    data = readDataForCreatingAbaqusModel()
    lengths = {key: data['length_'+key] for key in ['x', 'y', 'z']}
    dataKeysForPreliminaryModel = ['DRM_depth', 'PML_depth', 'partName', 'materialName', 'youngsModulus', 'poissonsRatio', 'density', 'timeIncrement', 'duration', 'jobName']
    createPreliminaryModel(lengths, *[data[key] for key in dataKeysForPreliminaryModel], noDamping=True, stepApplication='moderate dissipation')
    setViews()
    # dataKeysForModifyingInput = ['PML_depth', 'partName', 'materialName', 'alpha', 'beta', 'jobName', 'subroutineFileName']
    # modifyInputAndCheck(lengths, *[data[key] for key in dataKeysForModifyingInput])