import numpy as np
import pandas as pd
import sys
import csv
sys.path.append('../HerculesDatabaseInquirySystem')
from hdf5Search import getConvertedGridPointsForAbaqusModel, getInterpolatedHistoryDataForGridPoints, getDistanceFromPlaneOrigin
from getMaterialProperties import getMaterialProperties

def modifyInput(lengths, PML_depth, partName, materialName, alpha, beta, jobName, preInputFileName=None, dummyElementType='C3D8R', cLoadFileName='Cload.txt'):
    if preInputFileName is None:
        preInputFileName = jobName+'_pre.inp'
    density, youngsModulus, poissonsRatio = getMaterialPropertiesFromInputFile(preInputFileName.rstrip('.inp'), materialName)
    # Read the input file and make modifications for UEL use
    with open(preInputFileName, 'r') as f:
        lines = f.readlines()
        dummyElementLine = getLineIndex(lines, '*Element, type=%s\n'%dummyElementType, isConversionNeeded=True)
        lines[dummyElementLine] = '*Element, type=U3\n'
        # [NOTE] Info for parameters under *USER ELEMENT:
        #   TYPE: Must be 'Un' where n is a positive integer less than 10000, and it must be the same as the element type key used to identify this element on the *ELEMENT option.
        #   NODE: Number of nodes associated with an element of this type.
        #   COORDINATES: 3 for 3D, 2 for 2D.
        #   PROPERTY: The number of property values needed as data in UEL to define such an element.
        #   VARIABLES: 360 for 3D problems and 80 for 2D problems. This is determined by Wenyang's UEL subroutine.
        #   [TODO (maybe?)] So how to calculate the number of or variables (VARIABLES) that are needed for *USER ELEMENT?
        #   Numbers in the second line: The active DoF
        lines.insert(dummyElementLine, '*User Element, Type=U3, Nodes=8, Coordinates=3, Properties=12, Variables=360, Unsymm\n')
        lines.insert(dummyElementLine+1, '1, 2, 3, 21, 22, 23, 24, 25, 26\n')
        # NOTE: Since lines are changed above, we cannot create lowerLines like in other functions and reuse it. As a result, we need to create the lowerLines list again (with isConversionNeeded=True).
        endPartLine = getLineIndex(lines, '*End Part\n', isConversionNeeded=True)
        # [NOTE] Info for parameters used for UEL
        #   E: Young's Modulus
        #   xnu: Poisson Ratio
        #   rho: Density
        #   EleType_pos: 0.0 (A trivial property)
        #   PML_L: Depth of PML (distance from the boundary of DRM layer to the exterior of PML region)
        #   afp: Described as `m` and is equal to 2.0 in Wenyang's paper. It is the polynomial degree for alpha and beta functions.
        #   PML_Rcoef: Described as `R` and is equal to 10^(-10) in Wenyang's paper. It is a user-tunable reflection coefficient.
        #   RD_half_width_x: Half width of the domain in x direction without PML layers (i.e., the distance in x direction from the center to the boundary of DRM layer)
        #   RD_half_width_y: Half width of the domain in y direction without PML layers (i.e., the distance in y direction from the center to the boundary of DRM layer)
        #   RD_depth: The depth (in z direction) of the domain without PML (i.e., the depth of the interested region + the thickness of DRM layer)
        #   Damp_alpha, Damp_beta: alpha and beta used in Rayleigh Damping
        lines[endPartLine:endPartLine] = ['*Parameter\n', 'E=%e\n'%youngsModulus, 'xnu=%f\n'%poissonsRatio, 'rho=%f\n'%density, 
            'EleType_pos=0.0\n', 'PML_L=%f\n'%PML_depth, 'afp=2.0\n', 'PML_Rcoef=1e-10\n', 
            'RD_half_width_x=%f\n'%(lengths['x']/2-PML_depth), 'RD_half_width_y=%f\n'%(lengths['y']/2-PML_depth),
            'RD_depth=%f\n'%(lengths['z']-PML_depth), 'Damp_alpha=%f\n'%alpha, 'Damp_beta=%f\n'%beta, 
            '*UEL Property, elset=PML\n', '<E>, <xnu>, <rho>, <EleType_pos>, <PML_L>, <afp>, <PML_Rcoef>, <RD_half_width_x>,\n', 
            '<RD_half_width_y>, <RD_depth>, <Damp_alpha>, <Damp_beta>\n']
        # NOTE: Again, the lines are changed above, so we need to create the lowerLines list again.
        endStepLine = getLineIndex(lines, '*End Step\n', isConversionNeeded=True)
        lines.insert(endStepLine, '*Include, input=%s\n'%cLoadFileName)
    with open(jobName+'.inp', 'w') as f:
        f.writelines(lines)
    return

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

def getNextKeywordLine(lines, startLine=0):
    for line in lines[startLine:]:
        if line.startswith('*'):
            return getLineIndex(lines, line, startLine=startLine)
    return None

def getNodesOnAPart(jobName, partName):
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    partLine = getLineIndex(lowerLines, '*Part, name=%s\n'%partName)
    nodeLine = getLineIndex(lowerLines, '*Node\n', startLine=partLine)
    nextKeywordLine = getNextKeywordLine(lowerLines, nodeLine+1) # NOTE: the +1 is critical!
    nodeLines = lines[nodeLine+1:nextKeywordLine]
    nodes = [line.rstrip(',\n').split(',') for line in nodeLines]
    nodes = {int(node[0]): [float(x) for x in node[1:]] for node in nodes}
    return nodes

def getElementsOnAPart(jobName, partName, elementType):
    # TODO: elementType could be None and collect all elements regardless their element types
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    partLine = getLineIndex(lowerLines, '*Part, name=%s\n'%partName)
    elementLine = getLineIndex(lowerLines, '*Element, type=%s\n'%elementType, startLine=partLine)
    nextKeywordLine = getNextKeywordLine(lowerLines, elementLine+1) # NOTE: the +1 is critical!
    elementLines = lines[elementLine+1:nextKeywordLine]
    elements = [line.rstrip(',\n').split(',') for line in elementLines]
    elements = {int(element[0]): [int(x) for x in element[1:]] for element in elements}
    return elements

def getLabelsInSet(jobName, setName, setType):
    # NOTE: setType can be 'element' or 'node'
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    target = {'element': '*Elset, elset=%s\n'%setName, 
        'node': '*Nset, nset=%s\n'%setName}
    try:
        setLine = getLineIndex(lowerLines, target[setType])
        nextKeywordLine = getNextKeywordLine(lowerLines, setLine+1)
        labelLines = lines[setLine+1:nextKeywordLine]
        # NOTE: iterable unpacking cannot be used in comprehension
        labels = [int(label) for line in labelLines for label in line.rstrip(',\n').split(',')]
    except ValueError: # target is not in the input file
        target = target[setType].rstrip(',\n') + ', generate\n'
        setLine = getLineIndex(lowerLines, target)
        nextKeywordLine = getNextKeywordLine(lowerLines, setLine+1)
        labelLines = lines[setLine+1:nextKeywordLine]
        labelRanges = [[int(label) for label in line.rstrip(',\n').split(',')] for line in labelLines]
        labels = [label for lRange in labelRanges for label in range(lRange[0], lRange[1]+1, lRange[2])]
    return labels

def getShapeFunction(xi, eta, zeta):
    N = 1/8 * np.array([(1-xi)*(1-eta)*(1-zeta), (1+xi)*(1-eta)*(1-zeta), 
        (1+xi)*(1+eta)*(1-zeta), (1-xi)*(1+eta)*(1-zeta),
        (1-xi)*(1-eta)*(1+zeta), (1+xi)*(1-eta)*(1+zeta),
        (1+xi)*(1+eta)*(1+zeta), (1-xi)*(1+eta)*(1+zeta)])
    dN = 1/8 * np.array([[-(1-eta)*(1-zeta), (1-eta)*(1-zeta), 
        (1+eta)*(1-zeta), -(1+eta)*(1-zeta),
        -(1-eta)*(1+zeta), (1-eta)*(1+zeta),
        (1+eta)*(1+zeta), -(1+eta)*(1+zeta)],
        [-(1-xi)*(1-zeta), -(1+xi)*(1-zeta), 
        (1+xi)*(1-zeta), (1-xi)*(1-zeta),
        -(1-xi)*(1+zeta), -(1+xi)*(1+zeta),
        (1+xi)*(1+zeta), (1-xi)*(1+zeta)],
        [-(1-xi)*(1-eta), -(1+xi)*(1-eta), 
        -(1+xi)*(1+eta), -(1-xi)*(1+eta),
        (1-xi)*(1-eta), (1+xi)*(1-eta),
        (1+xi)*(1+eta), (1-xi)*(1+eta)]])
    return N, dN

def getMergeDicts(dict1, dict2):
    mergedDict = dict1.copy()
    mergedDict.update(dict2)
    return mergedDict

def getStiffnessMatrix(E, nu):
    G = E/(2*(1+nu)) # NOTE: G = mu
    C12 = E*nu/((1+nu)*(1-2*nu)) # C12 = lambda
    C11 = E*(1-nu)/((1+nu)*(1-2*nu))
    # NOTE: return the stiffness matrix for isotropic material.
    return(np.array([[C11, C12, C12, 0, 0, 0], 
        [C12, C11, C12, 0, 0, 0], 
        [C12, C12, C11, 0, 0, 0],
        [0, 0, 0, G, 0, 0],
        [0, 0, 0, 0, G, 0],
        [0, 0, 0, 0, 0, G]]))

def getJacobian(nodes, dN):
    """ dN is the derivative of shape function """
    if type(nodes) is dict:
        nodes = np.array(nodes.values())
    # NOTE: the definitions of A and J follow CEE 235B conducted by Prof. ET. 
    # Some others might use J = np.matmul(dN, nodes), and Bi will be changed accordingly.
    A = np.matmul(dN, nodes)
    J = A.transpose()
    if np.linalg.det(J) < 0:
        # NOTE: A nonzero volume element in the real element is mapped 
        # into zero volume in the master element, which is unacceptable.
        raise ValueError('Negative Jacobian determinant!')
    return J
    
def getElementStiffnessMatrix(D, J, dN, wXi, wEta, wZeta):
    dNdx = np.matmul(np.linalg.inv(J.transpose()), dN)
    # NOTE: Hard-coded 6 and 8 is used here. This only applies to 3D problem.
    Bi = [[] for i in range(6)]
    for i in range(8):
        Bi[0] += [dNdx[0, i], 0, 0]
        Bi[1] += [0, dNdx[1, i], 0]
        Bi[2] += [0, 0, dNdx[2, i]]
        # Bi[3] += [0, dNdx[2, i], dNdx[1, i]]
        # Bi[4] += [dNdx[2, i], 0, dNdx[0, i]]
        # Bi[5] += [dNdx[1, i], dNdx[0, i], 0]
        # NOTE: the definition of Bi follows CEE 235B conducted by Prof. ET. 
        # Some others might use the definition above. Then the definition of 
        # J will also be changed accordingly.
        Bi[3] += [dNdx[1, i], dNdx[0, i], 0]
        Bi[4] += [dNdx[2, i], 0, dNdx[0, i]]
        Bi[5] += [0, dNdx[2, i], dNdx[1, i]]
    B = np.array(Bi)
    J_det = np.linalg.det(J)
    # NOTE: @ is the special notation used for numpy to handle matrix multiplication. (same as np.matmul())
    # NOTE 2: Unfortunately, the version of numpy used in Abaqus might be too old to support it.
    eleK = np.matmul(np.matmul(B.transpose(), D), B)*J_det*wXi*wEta*wZeta
    return eleK

def getElementMassMatrix(density, N, J, wXi, wEta, wZeta):
    Ni = [[] for i in range(3)]
    for i in range(len(N)):
        Ni[0] += [N[i], 0, 0]
        Ni[1] += [0, N[i], 0]
        Ni[2] += [0, 0, N[i]]
    N = np.array(Ni)
    J_det = np.linalg.det(J)
    eleM = density*np.matmul(N.transpose(), N)*J_det*wXi*wEta*wZeta
    return eleM

def getGlobalMatrices(D, density, integrationPoints, weighting, DRM_elements, DRM_sideNodeLabels, nodes, alpha, beta):
    """ Returns K, M, C matrices for DRM elements """
    nDOF = 3*len(DRM_sideNodeLabels)
    drmK = np.zeros([nDOF, nDOF])
    drmM = np.zeros([nDOF, nDOF])
    drmC = np.zeros([nDOF, nDOF])
    for i, intP in enumerate(integrationPoints):
        N, dN = getShapeFunction(*intP)
        for eleLabel, nodeLabels in DRM_elements.items():
            J = getJacobian(np.array([nodes[label] for label in nodeLabels]), dN)
            eleK = getElementStiffnessMatrix(D, J, dN, *weighting[i])
            eleM = getElementMassMatrix(density, N, J, *weighting[i])
            eleC = alpha*eleK + beta*eleM
            nodeIndices = [DRM_sideNodeLabels.index(label) for label in nodeLabels]
            globalDOF = [label for nodeIndex in nodeIndices for label in range(3*nodeIndex, 3*nodeIndex+3)]
            for row, globalRow in enumerate(globalDOF):
                for col, globalCol in enumerate(globalDOF):
                    drmK[globalRow, globalCol] += eleK[row, col]
                    drmM[globalRow, globalCol] += eleM[row, col]
                    drmC[globalRow, globalCol] += eleC[row, col]
    return drmK, drmM, drmC

def getAndWriteDisplacementHistoryForDRM(jobName, partName, targetOrigin, dispHistoryFileName='DispHistory.csv', dbPath='../HerculesDatabaseInquirySystem/database/planedisplacements.hdf5'):
    # NOTE: mpirun works here since we use `getInterpolatedHistoryDataForGridPoints()` function
    nodes = getNodesOnAPart(jobName, partName)
    DRM_interiorNodeLabels = getLabelsInSet(jobName, setName='inDRM', setType='node')
    DRM_exteriorNodeLabels = getLabelsInSet(jobName, setName='sideDRM', setType='node')
    DRM_sideNodeLabels = DRM_interiorNodeLabels + DRM_exteriorNodeLabels
    DRM_nodes = [nodes[label] for label in DRM_sideNodeLabels]
    gridPoints = getConvertedGridPointsForAbaqusModel(dbPath, DRM_nodes, origin=targetOrigin, gridPointsInMeter=True)
    df = getInterpolatedHistoryDataForGridPoints(gridPoints, dbPath, pointLabelList=DRM_sideNodeLabels, gridPointsInMeter=True)
    df.to_csv(dispHistoryFileName)
    return df

def getHistoryOutputforDRM(dispHistoryFileName='DispHistory.csv'):
    """ This function reads the displacement history and compute velocity and 
    acceleration histories. Finally, it returns a Pandas DataFrame that contains 
    all 3 histories. """
    df = pd.read_csv(dispHistoryFileName, index_col=0)
    pointLabelList = df['pointLabel'].drop_duplicates().to_list()
    histories = {}
    for pointLabel in pointLabelList:
        history = df[df['pointLabel'] == pointLabel]
        history.columns = list(history.columns[:-3]) + ['ux', 'uy', 'uz']
        dt = history['time'].iloc[1] - history['time'].iloc[0]
        histories[pointLabel] = {}
        for direction in ['x', 'y', 'z']:
            history.insert(len(history.columns), 'v'+direction, history['u'+direction].diff()/dt)
            history.loc[history.index[0], 'v'+direction] = 0
            history.insert(len(history.columns), 'a'+direction, history['v'+direction].diff()/dt)
            history.loc[history.index[0], 'a'+direction] = 0
            for quantity in ['u', 'v', 'a']:
                histories[pointLabel][quantity+direction] = history[quantity+direction].to_numpy()
    return history['time'].to_list(), histories

def getQuantityMatrixAtPoints(histories, pointLabelList, quantity):
    matrix = [histories[pointLabel][quantity+direction] for pointLabel in pointLabelList for direction in ['x', 'y', 'z']]
    return np.array(matrix)

# NOTE: This function assume there is only one instance created from the given part
def getInstanceName(jobName, partName):
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    for line in lowerLines:
        if line.startswith('*instance') and line.endswith('part=%s\n'%partName):
            columns = [x.strip() for x in line.split(',')]
            for column in columns:
                if column.startswith('name='):
                    instanceName = column[5:]
                    return instanceName
    return None

def getEquivalentForces(jobName, partName, materialName, alpha, beta, cLoadFileName='Cload.txt'):
    density, youngsModulus, poissonsRatio = getMaterialPropertiesFromInputFile(jobName, materialName)
    # NOTE: The variable D below is usually denoted as C, but we also use C to denote the damping matrix.
    # To avoid the confusion, here we use D to denote the stiffness matrix for isotropic material.
    D = getStiffnessMatrix(youngsModulus, poissonsRatio)
    nodes = getNodesOnAPart(jobName, partName)
    elements = getElementsOnAPart(jobName, partName, elementType='C3D8')
    DRM_elementLabels = getLabelsInSet(jobName, setName='DRM', setType='element')
    DRM_elements = {label: elements[label] for label in DRM_elementLabels}
    DRM_interiorNodeLabels = getLabelsInSet(jobName, setName='inDRM', setType='node')
    DRM_exteriorNodeLabels = getLabelsInSet(jobName, setName='sideDRM', setType='node')
    DRM_sideNodeLabels = DRM_interiorNodeLabels + DRM_exteriorNodeLabels
    # NOTE: The simplest Gaussian quadrature is used here. Should be able to change it if needed.
    # For more information: https://en.wikipedia.org/wiki/Gaussian_quadrature
    integrationPoints = 1/np.sqrt(3) * np.array([[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1]])
    weighting = [[1, 1, 1] for i in range(8)]
    drmK, drmM, drmC = getGlobalMatrices(D, density, integrationPoints, weighting, DRM_elements, DRM_sideNodeLabels, nodes, alpha, beta)
    timePoints, histories = getHistoryOutputforDRM()
    ub = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'u')
    ue = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'u')
    vb = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'v')
    ve = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'v')
    ab = getQuantityMatrixAtPoints(histories, DRM_interiorNodeLabels, 'a')
    ae = getQuantityMatrixAtPoints(histories, DRM_exteriorNodeLabels, 'a')
    Kbe = drmK[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Keb = drmK[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    Mbe = drmM[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Meb = drmM[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    Cbe = drmC[range(3*len(DRM_interiorNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels))]
    Ceb = drmC[range(3*len(DRM_interiorNodeLabels), 3*len(DRM_sideNodeLabels)), :][:, range(3*len(DRM_interiorNodeLabels))]
    P_eff_b = -np.matmul(Mbe, ae) - np.matmul(Cbe, ve) - np.matmul(Kbe, ue)
    P_eff_e = np.matmul(Meb, ab) + np.matmul(Ceb, vb) - np.matmul(Keb, ub)
    # P_eff_bx = P_eff_b[range(0, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_by = P_eff_b[range(1, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_bz = P_eff_b[range(2, 3*len(DRM_interiorNodeLabels), 3)]
    # P_eff_ex = P_eff_e[range(0, 3*len(DRM_exteriorNodeLabels), 3)]
    # P_eff_ey = P_eff_e[range(1, 3*len(DRM_exteriorNodeLabels), 3)]
    # P_eff_ez = P_eff_e[range(2, 3*len(DRM_exteriorNodeLabels), 3)]
    startIndices = {'x': 0, 'y': 1, 'z': 2}
    DoF = {'x': 1, 'y': 2, 'z': 3}
    lines = []
    P_effs = np.vstack((P_eff_b, P_eff_e))
    instanceName = getInstanceName(jobName, partName)
    for direction, startIndex in startIndices.items():
        for nodeIndex, nodeLabel in enumerate(DRM_sideNodeLabels):
            lines += ['*Cload, amplitude=Amp-%d%s\n'%(nodeLabel, direction),
                '%s.%d, %s, 1\n'%(instanceName, nodeLabel, DoF[direction]),
                '*Amplitude, name=Amp-%d%s\n'%(nodeLabel, direction)]
            for timeIndex, P_eff in enumerate(P_effs[range(startIndex, 3*len(DRM_sideNodeLabels), 3)][nodeIndex]):
                lines.append('%.4f, %e\n'%(timePoints[timeIndex], P_eff))
    with open(cLoadFileName, 'w') as f:
        f.writelines(lines)
    return

# NOTE: This function is used to get and centroid of the part and later use it 
# to get the material properties at this centroid. However, it should be done 
# when building the model. So this function might be moved later.
# def getCentroid(jobName, partName):
#     nodes = getNodesOnAPart(jobName, partName)
#     nodeCoordinates = {'x': [], 'y': [], 'z': []}
#     for node in nodes.values():
#         for i, direction in enumerate(nodeCoordinates.keys()):
#             nodeCoordinates[direction].append(node[i])
#     # NOTE: Because the domain should be a simple hexagon, we can get the centroid 
#     # of it by the following easy way. Be aware if the domain is more complicated.
#     centroid = [max(nodeCoordinates[direction])+min(nodeCoordinates[direction]) for direction in nodeCoordinates.keys()]
#     return centroid

def getMaterialPropertiesAtCentroid(centroid, isAdjustmentNeeded=True, origin=None, dbPath='../HerculesDatabaseInquirySystem/database/planedisplacements.hdf5'):
    if isAdjustmentNeeded and origin is not None:
        planesData = pd.read_hdf(dbPath, 'planesData')
        distance, x_dis, y_dis = getDistanceFromPlaneOrigin(origin, planesData.loc[0])
        centroid = [centroid[0]+x_dis, centroid[1]+y_dis, centroid[2]]
    # NOTE: The input of `getMaterialProperties` should be a dict. The key is designed 
    # to be an element's tag (label). Here we just use a dummy 1 to make it works.
    materialProperties = getMaterialProperties({1: centroid})
    Vs, Vp, rho = materialProperties[1]
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

if __name__ == '__main__':
    # NOTE: Change working diretory to the folder of this script
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # NOTE: `lengths` is a dict with keys ('x', 'y', 'z') and the values are the 
    # total length of the part (including the interested domain, DRM layer, and PML layer).
    # lengths = {'x': 365, 'y': 190, 'z': 150}
    lengths = {'x': 100, 'y': 100, 'z': 100}
    DRM_depth = 5 # The thickness of the DRM layer should be the mesh size used in the Abaqus model
    PML_depth = DRM_depth*5 # Could be ranged from 5 to 10 times of the thickness of DRM. But 5 is good enough.
    partName = 'domain'
    materialName = 'soil'
    # ===== Rayleigh Damping =====
    alpha = 0
    beta = 0
    # =====
    jobName = 'DRM_PML_Analysis'
    subroutineFileName = 'PML3dUEL912.f'
    timeIncrement = 0.05 # unit: sec. This should be read from the site response file
    duration = 29.95 # unit: sec. This should be read from the site response file
    targetOrigin = [41.022024713821054, 28.88575077152881]

    # ===== Get material properties =====
    centroid = [x/2 for x in lengths.values()]
    youngsModulus, poissonsRatio, density = getMaterialPropertiesAtCentroid(centroid, origin=targetOrigin)

    # ===== Write data for creating Abaqus model =====
    # NOTE: All the lengths are the total length of the part (including the interested domain, DRM layer, and PML layer).
    dataForAbaqusModel = {'length_x': lengths['x'], 'length_y': lengths['y'], 'length_z': lengths['z'], 
        'DRM_depth': DRM_depth, 'PML_depth': PML_depth,
        'partName': partName, 'materialName': materialName, 
        'youngsModulus': youngsModulus, 'poissonsRatio': poissonsRatio, 'density': density, 
        'alpha': alpha, 'beta': beta, 
        'jobName': jobName, 'subroutineFileName': subroutineFileName,
        'timeIncrement': timeIncrement, 'duration': duration}
    # writeDataForCreatingAbaqusModel(dataForAbaqusModel)
    # NOTE: At this step, move 'dataForAbaqusModel.csv' to the Abaqus working directory and prepare the model.
    # After the _pre.inp file generated, move the _pre.inp file back.

    # ===== Modify the preliminary Abaqus input file =====
    # modifyInput(lengths, PML_depth, partName, materialName, alpha, beta, jobName)

    # ===== Getting the displacement histories for DRM nodes =====
    # NOTE: mpirun works here
    # import timeit
    # numExec = 1
    # print('Averaged Elapsed Time: %.2f secs' % (timeit.timeit(lambda: getAndWriteDisplacementHistoryForDRM(jobName, partName, targetOrigin), number=numExec)/numExec))
    # df = getAndWriteDisplacementHistoryForDRM(jobName, partName, targetOrigin)
    
    # ===== Compute the equivalent forces =====
    getEquivalentForces(jobName, partName, materialName, alpha, beta)
    # NOTE: After this step, move the .inp file and Cload.txt to the Abaqus working directory and run the .inp file