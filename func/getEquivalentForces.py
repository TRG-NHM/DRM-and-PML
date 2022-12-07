import numpy as np
# import pandas as pd
# import io
from time import time
from getNodes import getNodeCoordinates
from getLine import getNextLineIndexStartsWith
from getElements import getElements
from getLabelsInSet import getLabelsInSet
from getHistoryOutputForDRM import getHistoryOutputForDRMFromDispHistoryFile, getHistoryOutputForDRMFromStationFiles
from getInstanceName import getInstanceName
from FEM import getGlobalMatrices

def getQuantityMatrixAtPoints(histories: dict[int, dict[str, np.array]], pointLabelList: list[int], quantity: str) -> np.array:
    ''' getQuantityMatrixAtPoints returns the response (defined by the input `quantity`) histories for given points (defined by the input `pointLabelList`). 
    The input `quantity` can be 'u', 'v', or 'a'. '''
    matrix = [histories[pointLabel][quantity+direction] for pointLabel in pointLabelList for direction in ['x', 'y', 'z']]
    return np.array(matrix)

def getEquivalentForces(jobName: str, partName: str, cLoadFileName='Cload.txt', elementTypeOnPart='C3D8', 
        DRM_ElSetName='DRM', inDRM_NSetName='inDRM', outDRM_NSetName='outDRM', stationFolder='Stations',
        dispHistoryFileName=None, RFFileName=None, isCoordinateConverted=False, isHomogeneous=False,
        materialName=None, nodeTableFileName='nodeTable.csv', truncateTime=None) -> None:
    ''' getEquivalentForces write a file named `cLoadFileName` containing equivalent forces at DRM nodes. 
    If isHomogeneous==True, dispHistoryFileName and materialName must be defined. 
    If isHomogeneous==False, stationFolder and nodeTableFileName must be defined. '''
    tic = time()
    nodes = getNodeCoordinates(jobName, partName)
    elements = getElements(jobName, partName, elementType=elementTypeOnPart)
    DRM_elementLabels = getLabelsInSet(jobName, setName=DRM_ElSetName, setType='element')
    DRM_elements = {label: elements[label] for label in DRM_elementLabels}
    DRM_interiorNodeLabels = getLabelsInSet(jobName, setName=inDRM_NSetName, setType='node')
    DRM_exteriorNodeLabels = getLabelsInSet(jobName, setName=outDRM_NSetName, setType='node')
    DRM_sideNodeLabels = DRM_interiorNodeLabels + DRM_exteriorNodeLabels
    # TODO: Implement higher order quadrature rules
    # NOTE: The simplest Gaussian quadrature is used here. Should be able to change it if needed.
    # For more information: https://en.wikipedia.org/wiki/Gaussian_quadrature
    integrationPoints = 1/np.sqrt(3) * np.array([[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1]])
    weighting = [[1, 1, 1] for i in range(8)]
    print(time()-tic)
    drmK, drmM, drmC = getGlobalMatrices(integrationPoints, weighting, DRM_elements, DRM_sideNodeLabels, nodes, 
        isHomogeneous=isHomogeneous, jobName=jobName, materialName=materialName)
    # NOTE: when dispHistoryFileName is not None, it will turn to the old way to get history for DRM elements. 
    # Note that the old way does not consider whether the coordinate is converted.
    print(time()-tic)
    if dispHistoryFileName is None:
        timePoints, histories = getHistoryOutputForDRMFromStationFiles(stationFolder, isCoordinateConverted=isCoordinateConverted, nodeTableFileName=nodeTableFileName, truncateTime=truncateTime)
    else:
        timePoints, histories = getHistoryOutputForDRMFromDispHistoryFile(dispHistoryFileName)
    print(time()-tic)
    # TODO: Can MPI be implemented here?
    ub = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'u')
    ue = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'u')
    vb = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'v')
    ve = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'v')
    ab = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'a')
    ae = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'a')
    print(time()-tic)
    Kbe = drmK[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Keb = drmK[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    Mbe = drmM[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Meb = drmM[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    Cbe = drmC[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Ceb = drmC[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    print(time()-tic)
    P_eff_b = -np.matmul(Mbe, ae) - np.matmul(Cbe, ve) - np.matmul(Kbe, ue)
    P_eff_e = np.matmul(Meb, ab) + np.matmul(Ceb, vb) + np.matmul(Keb, ub)
    # P_eff_bx = P_eff_b[range(0, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_by = P_eff_b[range(1, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_bz = P_eff_b[range(2, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_ex = P_eff_e[range(0, 3*len(DRM_exteriorNodeLabels), 3)]
    # P_eff_ey = P_eff_e[range(1, 3*len(DRM_exteriorNodeLabels), 3)]
    # P_eff_ez = P_eff_e[range(2, 3*len(DRM_exteriorNodeLabels), 3)]
    print(time()-tic)
    startIndices = {'x': 0, 'y': 1, 'z': 2}
    DoF = {'x': 1, 'y': 2, 'z': 3}
    P_effs = np.vstack((P_eff_b, P_eff_e))
    instanceName = getInstanceName(jobName, partName)
    print(time()-tic)
    # /// Read reaction forces from the file
    RF = {}
    if RFFileName is not None:
        with open(RFFileName, 'r') as f:
            lines = f.readlines()
        lowerLines = [line.lower() for line in lines]
        CloadLineIndices = [i for i, line in enumerate(lowerLines) if line.startswith('*cload')]
        amplitudeNames = {}
        for i in CloadLineIndices:
            amplitudeName = lines[i].strip('\n').split(',')[1].split('=')[1]
            region, directionIndex, magnitude = lines[i+1].strip('\n').split(',')
            instanceNameInLine, nodeLabel = region.split('.')
            if instanceNameInLine.lower() != instanceName.lower():
                continue
            nodeLabel = int(nodeLabel)
            if nodeLabel not in amplitudeNames.keys():
                amplitudeNames[nodeLabel] = {}
            amplitudeNames[nodeLabel][int(directionIndex)] = amplitudeName
        for nodeLabel in amplitudeNames.keys():
            RF[nodeLabel] = {}
            for directionIndex, amplitudeName in amplitudeNames[nodeLabel].items():
                heading = '*Amplitude, name='+amplitudeName
                index = getNextLineIndexStartsWith(lowerLines, heading=heading, isConversionNeeded=False)
                staticStepDuration, RF[nodeLabel][directionIndex] = [float(x) for x in lines[index+1].strip('\n').split(',')[-2:]]
    else:
        staticStepDuration = 0
    print(time()-tic)
    # /// Build Cload.txt NOTE: this process takes a lot of time. Not sure what can be done to make it faster
    lines = []
    for direction, startIndex in startIndices.items():
        for nodeIndex, nodeLabel in enumerate(DRM_sideNodeLabels):
            lines += ['*Cload, amplitude=Amp-%d%s\n'%(nodeLabel, direction),
                '%s.%d, %s, 1\n'%(instanceName, nodeLabel, DoF[direction])]
            if staticStepDuration != 0:
                lines += ['*Amplitude, name=Amp-%d%s, time=TOTAL TIME\n'%(nodeLabel, direction),
                    '%.4f, %e\n'%(0, 0)]
            else:
                lines.append('*Amplitude, name=Amp-%d%s\n'%(nodeLabel, direction))
            # P_eff = pd.DataFrame({'timePoints': timePoints, 'P_eff': P_effs[range(startIndex, 3*len(DRM_sideNodeLabels), 3)][nodeIndex]})
            # if nodeLabel in RF.keys():
            #     P_eff['P_eff'] = P_eff['P_eff'] + RF[nodeLabel][DoF[direction]]
            # lines.append(P_eff.to_csv(header=False, index=False))
            # s = io.StringIO()
            # np.savetxt(s, P_eff.values, fmt=['%.4f', '%e'], delimiter=', ')
            # lines.append(s.getvalue())

            # P_eff = P_effs[range(startIndex, 3*len(DRM_sideNodeLabels), 3)][nodeIndex]
            # if nodeLabel in RF.keys():
            #     P_eff = P_eff + RF[nodeLabel][DoF[direction]]
            # lines += ['%.4f, %e\n'%(timePoints[timeIndex]+staticStepDuration, P) for timeIndex, P in enumerate(P_eff)]

            for timeIndex, P_eff in enumerate(P_effs[range(startIndex, 3*len(DRM_sideNodeLabels), 3)][nodeIndex]):
                if nodeLabel in RF.keys():
                    P_eff = P_eff + RF[nodeLabel][DoF[direction]]
                lines.append('%.4f, %e\n'%(timePoints[timeIndex]+staticStepDuration, P_eff))
    with open(cLoadFileName, 'w') as f:
        f.writelines(lines)
    print(time()-tic)
    return

if __name__ == '__main__':
    # NOTE: write material files first before run this script.
    # Meanwhile, run the static analysis and get RF.txt by running getReactionForce.py
    # ===== Abaqus model information =====
    jobName = 'Model_complete'
    partName = 'Part-Soil'
    nodeTableFileName = 'nodeTable.csv'
    # ===== =====
    getEquivalentForces(jobName, partName, elementTypeOnPart='C3D8', 
        DRM_ElSetName='Set-DRM-element', RFFileName='RF.txt',
        inDRM_NSetName='Set-DRM-inner-node', outDRM_NSetName='Set-DRM-outer-node',
        isCoordinateConverted=True, stationFolder='Stations_topo')
        # stationFolder='Stations_flat', truncateTime=[6.0, 11.99])
        # stationFolder='Istanbul_sim55/outputfiles/stations_fullDomain', truncateTime=[3.0, 8.99])

    # FOR DEBUGGING
    # import matplotlib.pyplot as plt
    # timePoints, histories = getHistoryOutputForDRMFromStationFiles('Stations_flat', isCoordinateConverted=True)
    
    # with open('ElementMaterial.txt', 'r') as f:
    #     lines = f.readlines()
    # indices = getAllLinesIndicesStartsWith(lines, heading='*sdjsidjo', startLine=5000, isConversionNeeded=True)
    # print(indices)