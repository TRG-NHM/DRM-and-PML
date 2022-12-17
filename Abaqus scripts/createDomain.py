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
    gravityScale = 1.515, isStaticStepIncluded=True, 
    isCoordinateConverted=False, isOriginAtCenter=False,
    noDamping=False, applyRampForce=False, stepApplication=None):
    if isCoordinateConverted:
        # NOTE: should we swap x and y?
        lengths = {'x': lengths['y'], 'y': lengths['x'], 'z': lengths['z']}
        zSign = -1
    else:
        zSign = 1
    # ===== Creating the part =====
    p, model = getCurrentDisplayedObjectAndModel()
    s = model.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.rectangle(point1=(0.0, 0.0), point2=(lengths['x'], lengths['y']))
    p = model.Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=lengths['z'])
    vp = session.viewports[session.currentViewportName]
    vp.setValues(displayedObject=p)
    del model.sketches['__profile__']
    # ===== Creating datum planes for partitioning =====
    principalPlanes = {'x': YZPLANE, 'y': XZPLANE, 'z': XYPLANE}
    for direction, length in lengths.items():
        if direction != 'z': # x and y
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth)
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth+DRM_depth)
        if direction == 'z' and isCoordinateConverted:
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth)
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=PML_depth+DRM_depth)
        else: # x and y + z if coordinate is not converted
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=length-PML_depth)
            p.DatumPlaneByPrincipalPlane(principalPlane=principalPlanes[direction], offset=length-PML_depth-DRM_depth)
    # ===== Partitioning =====
    for plane in p.datums.values():
        p.PartitionCellByDatumPlane(datumPlane=plane, cells=p.cells)
    # ===== Set up corresponding cell array =====
    cells = {}
    if isCoordinateConverted:
        cells['internal'] = p.cells.getByBoundingBox(xMin=PML_depth+DRM_depth, yMin=PML_depth+DRM_depth, zMin=PML_depth+DRM_depth,
            xMax=lengths['x']-PML_depth-DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
        cells['DRM'] = p.cells.getByBoundingBox(xMin=PML_depth, yMin=PML_depth, zMin=PML_depth,
            xMax=lengths['x']-PML_depth, yMax=lengths['y']-PML_depth)
    else:
        cells['internal'] = p.cells.getByBoundingBox(xMin=PML_depth+DRM_depth, yMin=PML_depth+DRM_depth, zMax=lengths['z']-PML_depth-DRM_depth,
            xMax=lengths['x']-PML_depth-DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
        cells['DRM'] = p.cells.getByBoundingBox(xMin=PML_depth, yMin=PML_depth, zMax=lengths['z']-PML_depth,
            xMax=lengths['x']-PML_depth, yMax=lengths['y']-PML_depth)
    cells['DRM'] = part.CellArray([cell for cell in cells['DRM'] if cell not in cells['internal']])
    cells['PML'] = part.CellArray([cell for cell in p.cells if cell not in cells['internal'] and cell not in cells['DRM']])
    # ===== Create sets =====
    sets = {cellName: p.Set(cells=cells[cellName], name=cellName) for cellName in cells.keys()}
    
    # ===== Create material properties =====
    material = model.Material(name=materialName)
    material.Density(table=((density, ), ))
    material.Elastic(table=((youngsModulus, poissonsRatio), ))
    # ===== Create section and section assignments =====
    section = model.HomogeneousSolidSection(name=materialName, material=materialName, thickness=None)
    [p.SectionAssignment(region=cellSet, sectionName=materialName) for cellName, cellSet in sets.items() if cellName != 'PML'] 

    # ===== Create the instance =====
    a = model.rootAssembly
    domainInstance = a.Instance(name=p.name+'-1', part=p, dependent=ON)
    if isCoordinateConverted:
        a.translate(instanceList=(domainInstance.name, ), vector=(0.0, 0.0, -lengths['z']))
    if isOriginAtCenter:
        a.translate(instanceList=(domainInstance.name, ), vector=(-lengths['x']/2, -lengths['y']/2, 0.0))
    
    # ===== Mesh =====
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
    # ===== Create element/node sets =====
    eleDRM = p.Set(elements=p.sets['DRM'].elements, name='eleDRM')
    if isCoordinateConverted:
        nodes_on_DRM_external_sides = p.sets['DRM'].nodes.getByBoundingBox(zMax=PML_depth) + \
            p.sets['DRM'].nodes.getByBoundingBox(xMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(xMin=lengths['x']-PML_depth) + \
            p.sets['DRM'].nodes.getByBoundingBox(yMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(yMin=lengths['y']-PML_depth)
        nodes_on_DRM_internal_sides = p.sets['DRM'].nodes.getByBoundingBox(zMin=PML_depth+DRM_depth,
            xMin=PML_depth+DRM_depth, xMax=lengths['x']-PML_depth-DRM_depth,
            yMin=PML_depth+DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
    else:
        nodes_on_DRM_external_sides = p.sets['DRM'].nodes.getByBoundingBox(zMin=lengths['z']-PML_depth) + \
            p.sets['DRM'].nodes.getByBoundingBox(xMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(xMin=lengths['x']-PML_depth) + \
            p.sets['DRM'].nodes.getByBoundingBox(yMax=PML_depth) + p.sets['DRM'].nodes.getByBoundingBox(yMin=lengths['y']-PML_depth)
        nodes_on_DRM_internal_sides = p.sets['DRM'].nodes.getByBoundingBox(zMax=lengths['z']-PML_depth-DRM_depth,
            xMin=PML_depth+DRM_depth, xMax=lengths['x']-PML_depth-DRM_depth,
            yMin=PML_depth+DRM_depth, yMax=lengths['y']-PML_depth-DRM_depth)
    outDRM = p.Set(nodes=nodes_on_DRM_external_sides, name='outDRM')
    inDRM = p.Set(nodes=nodes_on_DRM_internal_sides, name='inDRM')
    # ===== Prepare sets for history output =====
    node_at_center = p.nodes.getClosest((lengths['x']/2, lengths['y']/2, lengths['z']/2))
    # NOTE: Abaqus Set only takes MeshNodeArray instead of MeshNode.
    centerNode = p.Set(nodes=mesh.MeshNodeArray((node_at_center,)), name='centerNode')
    if isCoordinateConverted:
        node_at_top_center = p.nodes.getClosest((lengths['x']/2, lengths['y']/2, lengths['z']))
    else:
        node_at_top_center = p.nodes.getClosest((lengths['x']/2, lengths['y']/2, 0))
    topCenterNode = p.Set(nodes=mesh.MeshNodeArray((node_at_top_center,)), name='topCenterNode')

    # ===== Create the BC =====
    externalFaces = p.faces.getByBoundingBox(xMax=0) + p.faces.getByBoundingBox(xMin=lengths['x']) + \
        p.faces.getByBoundingBox(yMax=0) + p.faces.getByBoundingBox(yMin=lengths['y'])
    if isCoordinateConverted:
        externalFaces += p.faces.getByBoundingBox(zMax=0)
    else:
        externalFaces += p.faces.getByBoundingBox(zMin=lengths['z'])
    externalFacesSet = p.Set(faces=externalFaces, name='externalFaces')
    model.DisplacementBC(name='fixExteralFaces', createStepName='Initial', 
        region=domainInstance.sets['externalFaces'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
    
    if isStaticStepIncluded:
        # ===== Create static step =====
        staticStep = model.StaticStep(name='Step-Static', previous='Initial')
        prevStepName = staticStep.name
        # ===== Create additional BC =====
        BC_fixDRMNodes = model.DisplacementBC(name='fixDRMNodes', createStepName='Initial',
            region=domainInstance.sets['outDRM'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
        # ===== Create the Load =====
        amp = model.TabularAmplitude(name='Amp-1', timeSpan=TOTAL, 
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, gravityScale)))
        model.Gravity(name='Gravity', createStepName=staticStep.name, 
            comp3=zSign*9.81, amplitude=amp.name, distributionType=UNIFORM, field='')
            # region=domainInstance.sets['DRM_and_internal'])
        # ===== Create job and write the input file =====
        job = mdb.Job(name=jobName+'_static_pre', model=model.name)
        job.writeInput()
        mdb.saveAs(pathName='soil_static')
        # ===== Suppress the BC for DRM nodes for the complete analysis =====
        BC_fixDRMNodes.suppress()
        # NOTE: Since we are going to apply reaction forces extracted from the static analysis to 
        # the complete model in the static step, BC_fixDRMNodes should be released at the beginning. 
        # If we only want to deactivate it in the dynamic step, use the following line instead:
        # BC_fixDRMNodes.deactivate(dynamicStep.name) # Deactivate the BC for DRM nodes in dynamic step.
        # ===== Name for the complete_pre job =====
        completeJobName = jobName+'_complete_pre'
    else:
        prevStepName = 'Initial'
        completeJobName = jobName+'_pre'

    # ===== Create/Add dynamic implicit step =====
    dynamicStep = model.ImplicitDynamicsStep(name='Step-Dynamic', previous=prevStepName, timePeriod=duration, 
        maxNumInc=int(duration/timeIncrement*10), initialInc=timeIncrement, maxInc=timeIncrement, minInc=duration/100000.0)
    if stepApplication is not None:
        if stepApplication.lower() == 'transient fidelity':
            dynamicStep.setValues(application=TRANSIENT_FIDELITY, initialConditions=ON)
        elif stepApplication.lower() == 'moderate dissipation':
            dynamicStep.setValues(application=MODERATE_DISSIPATION, initialConditions=OFF)
        elif stepApplication.lower() == 'quasi static':
            dynamicStep.setValues(application=QUASI_STATIC, amplitude=RAMP)
        else:
            raise ValueError("The step application is not valid.")
    # else: # NOTE: This is not necessary. Just for reference.
    #     dynamicStep.setValues(application=DEFAULT, amplitude=STEP, initialConditions=DEFAULT)
    if noDamping:
        dynamicStep.setValues(alpha=0.0) # NOTE: the default value is calculated automatically
    if applyRampForce:
        dynamicStep.setValues(amplitude=RAMP) # NOTE: the default value is STEP

    # ===== Output request settings =====
    # NOTE: Change the output request frequency. Be aware if the model is not using a default output request setting.
    for fieldOutputRequest in model.fieldOutputRequests.values():
        fieldOutputRequest.setValuesInStep(stepName=dynamicStep.name, numIntervals=50, timeMarks=OFF)
    for historyOutputRequest in model.historyOutputRequests.values():
        historyOutputRequest.setValues(frequency=1)
    # ===== Create history outputs =====
    HOutputRequests = {'H-Output-Center': 'centerNode', 'H-Output-TopCenter': 'topCenterNode'}
    for HOutputName, nodeSetName in HOutputRequests.items():
        model.HistoryOutputRequest(name=HOutputName, createStepName=dynamicStep.name, 
            variables=('U1', 'U2', 'U3', 'V1', 'V2', 'V3', 'A1', 'A2', 'A3'), frequency=1, 
            region=domainInstance.sets[nodeSetName])

    # ===== Create job and write the input file =====
    job = mdb.Job(name=completeJobName, model=model.name)
    job.writeInput()
    mdb.saveAs(pathName='soil')
    return

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
                    if line[key].lower() == 'true':
                        data[key] = True
                    elif line[key].lower() == 'false':
                        data[key] = False
                    else:
                        data[key] = line[key]
            return data

if __name__ == '__main__':
    data = readDataForCreatingAbaqusModel()
    lengths = {key: data['length_'+key] for key in ['x', 'y', 'z']}
    dataKeysForPreliminaryModel = ['DRM_depth', 'PML_depth', 'partName', 'materialName', 
        'youngsModulus', 'poissonsRatio', 'density', 'timeIncrement', 'duration', 'jobName', 
        'noDamping', 'stepApplication', 'isStaticStepIncluded', 'isCoordinateConverted', 'isOriginAtCenter']
    createPreliminaryModel(lengths, **{key: data[key] for key in dataKeysForPreliminaryModel if key in data.keys()})
    setViews()