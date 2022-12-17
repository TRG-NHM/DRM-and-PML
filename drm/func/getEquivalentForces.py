import numpy as np
# from scipy.linalg.blas import dgemm
from time import time
from getNodes import getNodeCoordinates
# from getLine import getNextLineIndexStartsWith
from getElements import getElements
from getLabelsInSet import getLabelsInSet
from getHistoryOutputForDRM import getHistoryOutputForDRMFromDispHistoryFile, getHistoryOutputForDRMFromStationFiles
from getInstanceName import getInstanceName
from FEM import getGlobalMatrices

def getQuantityMatrixAtPoints(histories: dict[int, dict[str, np.array]], pointLabelList: list[int], quantity: str) -> np.array:
    ''' getQuantityMatrixAtPoints returns the response (defined by the input `quantity`) histories for given points (defined by the input `pointLabelList`). 
    The input `quantity` can be 'u', 'v', or 'a'. '''
    return np.array([histories[pointLabel][quantity+direction] for pointLabel in pointLabelList for direction in ['x', 'y', 'z']])

def getEquivalentForces(jobName: str, partName: str, cLoadFileName='Cload.txt', elementTypeOnPart='C3D8', 
        DRM_ElSetName='DRM', inDRM_NSetName='inDRM', outDRM_NSetName='outDRM', stationFolder='Stations',
        dispHistoryFileName=None, RFFileName=None, isCoordinateConverted=False, isHomogeneous=False,
        materialName=None, nodeTableFileName='nodeTable.csv', truncateTime=None, isResponseRecalculationNeeded=False) -> None:
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
    print('Preliminaries: %f sec'%(time()-tic))
    drmK, drmM, drmC = getGlobalMatrices(integrationPoints, weighting, DRM_elements, DRM_sideNodeLabels, nodes, 
        isHomogeneous=isHomogeneous, jobName=jobName, materialName=materialName)
    # NOTE: when dispHistoryFileName is not None, it will turn to the old way to get history for DRM elements. 
    # Note that the old way does not consider whether the coordinate is converted.
    print('Get global matrices: %f sec'%(time()-tic))
    if dispHistoryFileName is None:
        timePoints, histories = getHistoryOutputForDRMFromStationFiles(stationFolder, isCoordinateConverted=isCoordinateConverted, 
            nodeTableFileName=nodeTableFileName, truncateTime=truncateTime, isResponseRecalculationNeeded=isResponseRecalculationNeeded)
    else:
        timePoints, histories = getHistoryOutputForDRMFromDispHistoryFile(dispHistoryFileName)
    print('Get history output for DRM: %f sec'%(time()-tic))
    ub = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'u')
    ue = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'u')
    vb = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'v')
    ve = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'v')
    ab = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'a')
    ae = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'a')
    # NOTE: After this point, we are not using histories variable anymore.
    del histories # release some memory
    Kbe = drmK[:3*len(DRM_interiorNodeLabels), 3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels)]
    Keb = drmK[3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels), :3*len(DRM_interiorNodeLabels)]
    Mbe = drmM[:3*len(DRM_interiorNodeLabels), 3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels)]
    Meb = drmM[3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels), :3*len(DRM_interiorNodeLabels)]
    Cbe = drmC[:3*len(DRM_interiorNodeLabels), 3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels)]
    Ceb = drmC[3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels), :3*len(DRM_interiorNodeLabels)]
    print('Set up response vectors and K, M, and C matrices: %f sec'%(time()-tic))
    # Option 3: memory efficient way
    P_effs = np.zeros((3*len(DRM_sideNodeLabels), len(timePoints)))
    P_effs[:3*len(DRM_interiorNodeLabels), :] = -np.matmul(Mbe, ae) - np.matmul(Cbe, ve) - np.matmul(Kbe, ue) # P_eff_b
    P_effs[3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels), :] = np.matmul(Meb, ab) + np.matmul(Ceb, vb) + np.matmul(Keb, ub) # P_eff_e
    # Option 3.2: Scipy implementation (similar performance on my machine)
    # P_effs[:3*len(DRM_interiorNodeLabels), :] = -dgemm(1.0, Mbe, ae) - dgemm(1.0, Cbe, ve) - dgemm(1.0, Kbe, ue) # P_eff_b
    # P_effs[3*len(DRM_interiorNodeLabels):3*len(DRM_sideNodeLabels), :] = dgemm(1.0, Meb, ab) + dgemm(1.0, Ceb, vb) + dgemm(1.0, Keb, ub) # P_eff_e
    # Option 2:
    # P_eff_b = -np.matmul(Mbe, ae) - np.matmul(Cbe, ve) - np.matmul(Kbe, ue)
    # P_eff_e = np.matmul(Meb, ab) + np.matmul(Ceb, vb) + np.matmul(Keb, ub)
    # P_effs = np.vstack((P_eff_b, P_eff_e)) # np.vstack creates a new array. Not memory efficient.
    # Option 1: Just left for reference
    # P_eff_bx = P_eff_b[0:3*len(DRM_interiorNodeLabels):3]
    # P_eff_by = P_eff_b[1:3*len(DRM_interiorNodeLabels):3]
    # P_eff_bz = P_eff_b[2:3*len(DRM_interiorNodeLabels):3]
    # P_eff_ex = P_eff_e[0:3*len(DRM_exteriorNodeLabels):3]
    # P_eff_ey = P_eff_e[1:3*len(DRM_exteriorNodeLabels):3]
    # P_eff_ez = P_eff_e[2:3*len(DRM_exteriorNodeLabels):3]
    del ub, ue, vb, ve, ab, ae, Kbe, Keb, Mbe, Meb, Cbe, Ceb # release memory
    startIndices = {'x': 0, 'y': 1, 'z': 2}
    DoF = {'x': 1, 'y': 2, 'z': 3}
    instanceName = getInstanceName(jobName, partName)
    # /// Read reaction forces from the file
    RF = {}
    if RFFileName is not None:
        with open(RFFileName, 'r') as f:
            lines = f.readlines()
        lowerLines = [line.lower() for line in lines]
        CloadLineIndices = [i for i, line in enumerate(lowerLines) if line.startswith('*cload')]
        for i in CloadLineIndices:
            amplitudeName = lines[i].strip('\n').split(',')[1].split('=')[1]
            region, directionIndex, magnitude = lines[i+1].strip('\n').split(',')
            instanceNameInLine, nodeLabel = region.split('.')
            if instanceNameInLine.lower() != instanceName.lower():
                continue
            nodeLabel = int(nodeLabel)
            if nodeLabel not in RF.keys():
                RF[nodeLabel] = {}
            # NOTE: This approach only applies to the RF.txt written by getReactionForce.py.
            staticStepDuration, RF[nodeLabel][int(directionIndex)] = [float(x) for x in lines[i+3].strip('\n').split(',')[-2:]]
            # NOTE: For a general input file that contains *Cload and *Amplitude at uncertain lines, use the following lines.
            # heading = '*Amplitude, name='+amplitudeName
            # index = getNextLineIndexStartsWith(lowerLines, heading=heading, isConversionNeeded=False)
            # staticStepDuration, RF[nodeLabel][int(directionIndex)] = [float(x) for x in lines[index+1].strip('\n').split(',')[-2:]]
    else:
        staticStepDuration = 0
    print('Get reaction force: %f sec'%(time()-tic))
    # /// Build Cload.txt
    lines = ['']*(3+len(timePoints))
    timePoints = [timePoint+staticStepDuration for timePoint in timePoints]
    with open(cLoadFileName, 'w') as f:
        for direction, startIndex in startIndices.items():
            for nodeIndex, nodeLabel in enumerate(DRM_sideNodeLabels):
                lines[0:2] = ['*Cload, amplitude=Amp-%d%s\n'%(nodeLabel, direction),
                    '%s.%d, %s, 1\n'%(instanceName, nodeLabel, DoF[direction])]
                if staticStepDuration != 0:
                    lines[2] = '*Amplitude, name=Amp-%d%s, time=TOTAL TIME\n%.4f,0\n'%(nodeLabel, direction, 0)
                else:
                    lines[2] = '*Amplitude, name=Amp-%d%s\n'%(nodeLabel, direction)
                # NOTE: P_effs[startIndex:3*len(DRM_sideNodeLabels):3, :] returns all effective forces in one direction (x, y, or z)
                # P_eff = P_effs[startIndex:3*len(DRM_sideNodeLabels):3, :][nodeIndex, :] # NOTE: This is (just a tad) less efficient
                P_eff = P_effs[startIndex+3*nodeIndex, :]
                if nodeLabel in RF.keys():
                    P_eff += RF[nodeLabel][DoF[direction]]
                lines[3:] = ['%.4f,%e\n'%(timePoints[timeIndex], P) if P != 0 else '%.4f,0\n'%timePoints[timeIndex] for timeIndex, P in enumerate(P_eff)]
                f.writelines(lines)
    print('Build Cload.txt: %f sec'%(time()-tic))
    return