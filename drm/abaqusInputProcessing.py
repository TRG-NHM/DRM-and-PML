import csv
from material.getMaterialProperties import getMaterialProperties
from func.getDataFromInputFile import getMaterialPropertiesFromInputFile
from func.getLine import getLineIndex, getNextKeywordLine, getNextLineIndexStartsWith, getAllLinesIndicesStartsWith
from func.getPMLElementInputLine import getUserElementLines, getParametersLines
from func.getDistance import getDistanceBetweenTwoCoordinates

def modifyInput(jobName, partName, materialName, lengths, PML_depth, alpha, beta, 
    isHomogeneous=True, stepName=None, dummyElementType='C3D8R', userElementType='U3',
    preInputFileName=None, modelType=None, cLoadFileName='Cload.txt',
    materialFileName='ElementMaterial.txt', elementSetFileName='ElementSet.txt', 
    sectionFileName='ElementSection.txt', PMLMaterialFileName='ElementPML.txt',
    boundaryName=None, **kwargs):
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
        density, youngsModulus, poissonsRatio = getMaterialPropertiesFromInputFile(materialName, fileName=preInputFileName)[:3]
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
    if modelType != 'static':
        # Insert the step-related file (cLoad file)
        if stepName is None:
            stepLine = getNextLineIndexStartsWith(lines, heading='*Step')
        else:
            stepLine = getNextLineIndexStartsWith(lines, heading='*Step, name=%s'%stepName)
        endStepLine = getLineIndex(lines, '*End Step\n', startLine=stepLine+1)
        lines.insert(endStepLine, '*Include, input=%s\n'%cLoadFileName)
        # Deactivate the boundary conditions in the dynamic step
        if boundaryName is None:
            boundaryPreLine = getNextLineIndexStartsWith(lines, heading='** BOUNDARY CONDITIONS')
            boundaryLine = getNextLineIndexStartsWith(lines, heading='*Boundary', startLine=boundaryPreLine)
        else:
            boundaryNameLine = getNextLineIndexStartsWith(lines, heading='** Name: %s'%boundaryName)
            boundaryPreLine = getNextLineIndexStartsWith(lines, heading='** BOUNDARY CONDITIONS', startLine=boundaryNameLine, isReversed=True)
            boundaryLine = getNextLineIndexStartsWith(lines, heading='*Boundary', startLine=boundaryNameLine)
        boundaryEndLine = getNextKeywordLine(lines, startLine=boundaryLine+1)
        boundaryLines = lines[boundaryPreLine:boundaryEndLine].copy()
        boundaryLine = getNextLineIndexStartsWith(boundaryLines, heading='*Boundary')
        boundaryLines[boundaryLine] = '*Boundary, op=NEW\n'
        dynamicStepLine = getNextLineIndexStartsWith(lines, heading='*Dynamic', startLine=stepLine+1)
        endDynamicStepLine = getNextKeywordLine(lines, startLine=dynamicStepLine+1)
        lines[endDynamicStepLine:endDynamicStepLine] = boundaryLines
    # /// Save inp file
    if modelType is None:
        inpName = jobName+'.inp'
    else:
        inpName = jobName+'_'+modelType+'.inp'
    with open(inpName, 'w') as f:
        f.writelines(lines)
    return

def getMaterialPropertiesAtCentroid(centroid, isAdjustmentNeeded=True, origin=None, 
    istanbulMaterialModelOrigin=(40.9565, 28.675), minVs=None, **kwargs):
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
    # NOTE: Negative Poisson's ratio is usually caused by the replacement of a very small Vs with minVs.
    if poissonsRatio < 0:
        raise ValueError('Computed Poisson\'s Ratio for %s is negative: %.4f'%(str(origin), poissonsRatio))
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
