import csv
from materialDatabase.getMaterialProperties import getMaterialProperties
from func.getDataFromInputFile import getDataFromHerculesInputFile, getMaterialPropertiesFromInputFile
from func.getLine import getLineIndex, getNextKeywordLine, getNextLineIndexStartsWith, getAllLinesIndicesStartsWith
from func.getPMLElementInputLine import getUserElementLines, getParametersLines
from func.getDistance import getDistanceBetweenTwoCoordinates
from func.getEquivalentForces import getEquivalentForces
from func.getHistoryOutputForDRM import getAndWriteDisplacementHistoryForDRM
from func.getNodes import writeNodeTable
from writeMaterialProperties import writeMaterialPropertiesForElements, writeMaterialPropertiesForPMLElements

def modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, 
    isHomogeneous=True, stepName=None, dummyElementType='C3D8R', userElementType='U3',
    preInputFileName=None, modelType=None, cLoadFileName='Cload.txt',
    materialFileName='ElementMaterial.txt', elementSetFileName='ElementSet.txt', 
    sectionFileName='ElementSection.txt', PMLMaterialFileName='ElementPML.txt'):
    ''' Modify Abaqus input file. '''
    if preInputFileName is None:
        if modelType is None:
            preInputFileName = jobName+'_pre.inp'
        else:
            preInputFileName = jobName+'_'+modelType+'_pre.inp'
    # Read the input file and make modifications for UEL use
    with open(preInputFileName, 'r') as f:
        lines = f.readlines()
    partLine = getLineIndex(lines, '*Part, name=%s\n'%partName)
    dummyElementLine = getLineIndex(lines, '*Element, type=%s\n'%dummyElementType, startLine=partLine+1)
    if isHomogeneous:
        density, youngsModulus, poissonsRatio = getMaterialPropertiesFromInputFile(materialName, fileName=jobName+'.inp')[:3]
        lines[dummyElementLine] = '*Element, type=%s\n'%userElementType
        lines[dummyElementLine:dummyElementLine] = getUserElementLines(userElementType)
        endPartLine = getLineIndex(lines, '*End Part\n', startLine=partLine+1)
        lines[endPartLine:endPartLine] = getParametersLines(youngsModulus, poissonsRatio, density, PML_depth, lengths, alpha, beta)
    else: # heterogeneous
        # Remove dummy element section since it would be replaced with PMLMaterialFile
        endDummyElementLine = getNextKeywordLine(lines, startLine=dummyElementLine+1)
        del lines[dummyElementLine:endDummyElementLine]
        # Remove material section assignments since it would be replaced with elementSetFile and sectionFile
        endPartLine = getLineIndex(lines, '*End Part\n', startLine=partLine+1)
        sectionLines = getAllLinesIndicesStartsWith(lines[:endPartLine], heading='*Solid Section,', startLine=partLine+1)
        for sectionLine in reversed(sectionLines):
            endSectionLine = getNextKeywordLine(lines, startLine=sectionLine+1)
            del lines[sectionLine:endSectionLine]
        # Insert part-related files
        partLine = getLineIndex(lines, '*Part, name=%s\n'%partName)
        endPartLine = getLineIndex(lines, '*End Part\n', startLine=partLine+1)
        lines.insert(endPartLine, '*Include, input=%s\n'%elementSetFileName)
        lines.insert(endPartLine+1, '*Include, input=%s\n'%sectionFileName)
        lines.insert(endPartLine+2, '*Include, input=%s\n'%PMLMaterialFileName)
        # Insert material property file
        endAssemblyLine = getLineIndex(lines, '*End Assembly\n')
        lines.insert(endAssemblyLine+1, '*Include, input=%s\n'%materialFileName)
    # Insert the step-related file (cLoad file)
    if modelType != 'static':
        if stepName is None:
            stepLine = getNextLineIndexStartsWith(lines, heading='*Step')
        else:
            stepLine = getNextLineIndexStartsWith(lines, heading='*Step, name=%s'%stepName)
        endStepLine = getLineIndex(lines, '*End Step\n', startLine=stepLine+1)
        lines.insert(endStepLine, '*Include, input=%s\n'%cLoadFileName)
    # /// Save inp file
    if modelType is None:
        inpName = jobName+'.inp'
    else:
        inpName = jobName+'_'+modelType+'.inp'
    with open(inpName, 'w') as f:
        f.writelines(lines)
    return

def getMaterialPropertiesAtCentroid(centroid, isAdjustmentNeeded=True, origin=None, 
    istanbulMaterialModelOrigin=(40.9565, 28.675), minVs=None):
    if isAdjustmentNeeded and origin is not None:
        distance, x_dis, y_dis = getDistanceBetweenTwoCoordinates(istanbulMaterialModelOrigin, origin)
        centroid = [centroid[0]+x_dis, centroid[1]+y_dis, centroid[2]]
    # NOTE: The input of `getMaterialProperties` should be a dict. The key is designed 
    # to be an element's tag (label). Here we just use a dummy 1 to make it works.
    materialProperties = getMaterialProperties({1: centroid})
    Vs, Vp, rho = materialProperties[1]
    # NOTE: Applying this may make the material properties weird (even with some negative Poisson's Ratios). Use with cautions.
    if minVs is not None and Vs < minVs:
        Vs = minVs
    poissonsRatio = (Vp**2 - 2*Vs**2)/(2*(Vp**2 - Vs**2)) # Ref: Eq. (5.34) in Geotechnical Earthquake Engineering (Kramer)
    G = rho*Vs**2 # Shear Modulus
    youngsModulus = 2*G*(1+poissonsRatio)
    return youngsModulus, poissonsRatio, rho

def writeDataForCreatingAbaqusModel(dataForAbaqusModel, fileName='dataForAbaqusModel.csv'):
    ''' dataForAbaqusModel should be a dict '''
    with open(fileName, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=dataForAbaqusModel.keys())
        writer.writeheader()
        writer.writerow(dataForAbaqusModel)
    return

def prospectusModel(stepNum: int) -> None:
    ''' step 1: writeDataForCreatingAbaqusModel
        NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare 
        the model (run `createDomain.py`. abaqus cae noGUI=createDomain.py can do the trick).
        After the _pre.inp file generated, move the _pre.inp file back.
        step 2: modifyInput 
        step 3: getAndWriteDisplacementHistoryForDRM (MPI can be used in this step)
        step 4: getEquivalentForces 
        NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file '''
    dbPath = '../HerculesDatabaseInquirySystem/database/planedisplacements.hdf5'
    # ===== Abaqus model information =====
    # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
    # total length of the part (including the interested domain, DRM layer, and PML layer).
    lengths = {'x': 100, 'y': 100, 'z': 100}
    DRM_depth = 5 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
    PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    partName = 'Part-Soil'
    materialName = 'soil'
    # ===== Rayleigh Damping =====
    alpha = 0
    beta = 0
    # =====
    jobName = 'DRM_PML_Analysis'
    subroutineFileName = 'PML3dUEL912.f'
    timeIncrement = 0.05 # unit: sec. This should be read from the site response file
    duration = 29.95 # unit: sec. This should be read from the site response file
    origin = (41.022024713821054, 28.88575077152881)
    if stepNum == 1:
        # ===== Hercules model information =====
        HerculesInputFilePath = '../Mw5p83/inputfiles/parameters_ZDistrictDRM.in'
        HerculesData = getDataFromHerculesInputFile(HerculesInputFilePath)
        minVs = float(HerculesData['simulation_shear_velocity_min'])
        # ===== Get material properties =====
        centroid = [x/2 for x in lengths.values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(centroid, origin=origin, minVs=minVs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_x': lengths['x'], 'length_y': lengths['y'], 'length_z': lengths['z'], 
            'DRM_depth': DRM_depth, 'PML_depth': PML_depth,
            'partName': partName, 'materialName': materialName, 
            'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density, 
            'alpha': alpha, 'beta': beta, 
            'jobName': jobName, 'subroutineFileName': subroutineFileName,
            'timeIncrement': timeIncrement, 'duration': duration}
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # For homogeneous model without static step:
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta)
    elif stepNum == 3:
        # ===== Getting the displacement histories for DRM nodes =====
        # NOTE: mpirun works here
        import timeit
        numExec = 1
        print('Averaged Elapsed Time: %.2f secs' % (timeit.timeit(lambda: getAndWriteDisplacementHistoryForDRM(dbPath, jobName, partName, origin), number=numExec)/numExec))
        # df = getAndWriteDisplacementHistoryForDRM(jobName, partName, origin)
    elif stepNum == 4:
        # ===== Compute the equivalent forces =====
        getEquivalentForces(jobName, partName, isHomogeneous=True, dispHistoryFileName='DispHistory.csv', materialName=materialName)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file
    else:
        raise ValueError('There are only 4 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., DispHistory.csv generated based on HDF5 database).
    return

def IstanbulModel(stepNum: int) -> None:
    ''' step 1: writeDataForCreatingAbaqusModel
        NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare 
        the model (run `createDomain.py`. abaqus cae noGUI=createDomain.py can do the trick).
        After the _pre.inp files generated, move the _pre.inp files back.
        step 2: modifyInput 
        NOTE: After this step, move the .inp files (for inhomogeneous model, also move ElementSet.txt, 
        ElementSection.txt, ElementMaterial.txt, ElementPML.txt, and PML3dUEL_Inhomogeneous.for manually 
        modified with maxVp in maxVpInPMLDomain.txt) to the Abaqus working directory and run the static 
        analysis and get RF.txt by running getReactionForce.py, then move RF.txt back.
        In addition, the generated nodeTable.csv should be used to specify stations in Hercules' input 
        file and run Hercules analysis to generate station files.
        step 3: getEquivalentForces 
        NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file '''
    # ===== Abaqus model information =====
    # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
    # total length of the part (including the interested domain, DRM layer, and PML layer).
    lengths = {'x': 100, 'y': 100, 'z': 30}
    DRM_depth = 2 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
    PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    partName = 'Part-Soil'
    materialName = 'soil'
    # ===== Rayleigh Damping =====
    alpha = 0
    beta = 0
    # =====
    jobName = 'Istanbul_model'
    subroutineFileName = 'PML3dUEL912.f'
    timeIncrement = 0.01 # unit: sec. This should be read from the site response file
    duration = 39.99 # unit: sec. This should be read from the site response file
    origin = (41.0318, 28.9417)
    # ===== Hercules model information =====
    HerculesInputFilePath = '../Abaqus Steel Building Model from Bulent/Istanbul_sim55/inputfiles/parameters_FullRegion_all_station_topo_final.in'
    HerculesData = getDataFromHerculesInputFile(HerculesInputFilePath)
    minVs = float(HerculesData['simulation_shear_velocity_min'])
    if stepNum == 1:
        # ===== Get material properties =====
        centroid = [x/2 for x in lengths.values()]
        youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(centroid, origin=origin, minVs=minVs)
        # ===== Write data for creating Abaqus model =====
        # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
        dataForAbaqusModel = {'length_x': lengths['x'], 'length_y': lengths['y'], 'length_z': lengths['z'], 
            'DRM_depth': DRM_depth, 'PML_depth': PML_depth,
            'partName': partName, 'materialName': materialName, 
            'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density, 
            'alpha': alpha, 'beta': beta, 
            'jobName': jobName, 'subroutineFileName': subroutineFileName,
            'timeIncrement': timeIncrement, 'duration': duration}
        writeDataForCreatingAbaqusModel(dataForAbaqusModel)
        # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
        # After the _pre.inp file generated, move the _pre.inp file back.
    elif stepNum == 2:
        # # For homogeneous model:
        # # ===== Modify the preliminary Abaqus input file =====
        # modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, modelType='static')
        # modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, modelType='complete')

        # For heterogeneous model
        # ===== Writing Material Properties =====
        writeMaterialPropertiesForElements(jobName+'_static_pre', partName, origin, isCoordinateConverted=True, minVs=minVs)
        writeMaterialPropertiesForPMLElements(jobName+'_static_pre', partName, origin, PML_depth, 
            lengths, isCoordinateConverted=True, minVs=minVs)
        # ===== Modify the preliminary Abaqus input file =====
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, isHomogeneous=False, modelType='static')
        modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, isHomogeneous=False, modelType='complete')
        # ===== Write nodeTable =====
        writeNodeTable(jobName+'_static', partName, coordinateCorrection=True, AbaqusModelOrigin=origin, isCoordinateConverted=True)
        # NOTE: This node table should be used to specify stations in Hercules' input file and run Hercules analysis to generate station files.
    elif stepNum == 3:
        # ===== Compute the equivalent forces =====
        stationFolder = '../Abaqus Steel Building Model from Bulent/Stations_flat/'
        getEquivalentForces(jobName+'_complete', partName, stationFolder=stationFolder, RFFileName='RF.txt', isCoordinateConverted=True)
        # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file.
    else:
        raise ValueError('There are only 3 steps.')
    # NOTE: After completing the Abaqus analysis, open the .odb file and run plotAndSave.py in Abaqus to get the desired results.
    # Additionally and optionally, verification.py can be used to plot the results above and compare it with the ground truths (e.g., Hercules' station files).
    return

if __name__ == '__main__':
    # NOTE: Change working diretory to the folder of this script
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # prospectusModel(1)
    IstanbulModel(3)