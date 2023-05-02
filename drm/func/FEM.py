import numpy as np
from getDataFromInputFile import getMaterialPropertiesFromInputFile
# from get_MPI_data import get_MPI_data

def getShapeFunction(xi: float, eta: float, zeta: float) -> tuple[np.array, np.array]:
    ''' getShapeFunction returns shape functions, N (1-by-8 array), and its derivative, dN. 
    dN is a 3-by-8 matrix, three rows from the top are dN_i/dxi, dN_i/deta, dN_i/dzeta. '''
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

def getJacobian(nodes: np.array, dN: np.array) -> np.array:
    """ dN is the derivative of shape function """
    # NOTE: Forget why nodes might be a dict. Comment it out for now.
    # if type(nodes) is dict:
    #     nodes = np.array(nodes.values())
    # NOTE: the definitions of A and J follow CEE 235B conducted by Prof. ET. 
    # Some others might use J = np.matmul(dN, nodes), and Bi will be changed accordingly.
    A = np.matmul(dN, nodes)
    J = A.transpose()
    if np.linalg.det(J) < 0:
        # NOTE: A nonzero volume element in the real element is mapped 
        # into zero volume in the master element, which is unacceptable.
        raise ValueError('Negative Jacobian determinant!')
    return J

def getStiffnessMatrix(E: float, nu: float) -> np.array:
    ''' getStiffnessMatrix returns the stiffness matrix, D, for isotropic material. 
    Inputs: E is the Young's modulus, nu is the poisson's ratio. '''
    G = E/(2*(1+nu)) # NOTE: G = mu
    C12 = E*nu/((1+nu)*(1-2*nu)) # C12 = lambda
    C11 = E*(1-nu)/((1+nu)*(1-2*nu))
    return(np.array([[C11, C12, C12, 0, 0, 0], 
        [C12, C11, C12, 0, 0, 0], 
        [C12, C12, C11, 0, 0, 0],
        [0, 0, 0, G, 0, 0],
        [0, 0, 0, 0, G, 0],
        [0, 0, 0, 0, 0, G]]))

def getElementStiffnessMatrix(D: np.array, J: np.array, dN: np.array, wXi: float, wEta: float, wZeta: float) -> np.array:
    ''' getElementStiffnessMatrix returns 3D element stiffness matrix. '''
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
    # eleK = B.transpose()@D@B *J_det*wXi*wEta*wZeta
    # NOTE 2: The version of numpy used in Abaqus is too old to support it. 
    # If this is going to be integrated into Abaqus, use the following version.
    eleK = np.matmul(np.matmul(B.transpose(), D), B)*J_det*wXi*wEta*wZeta
    return eleK

def getElementMassMatrix(density: float, N: np.array, J: np.array, wXi: float, wEta: float, wZeta: float) -> np.array:
    ''' getElementMassMatrix returns 3D element mass (inertia) matrix. '''
    Ni = [[] for i in range(3)]
    for i in range(len(N)):
        Ni[0] += [N[i], 0, 0]
        Ni[1] += [0, N[i], 0]
        Ni[2] += [0, 0, N[i]]
    N = np.array(Ni)
    J_det = np.linalg.det(J)
    eleM = density*np.matmul(N.transpose(), N)*J_det*wXi*wEta*wZeta
    return eleM

def get3DGaussQuadPoints(order: int) -> tuple[list, list]:
    ''' getGaussQuadPoints returns the Gauss quadrature points and weights. '''
    if order == 1:
        integrationPoints = [0]
        weights = [2]
    elif order == 2:
        integrationPoints = [-1/np.sqrt(3), 1/np.sqrt(3)]
        weights = [1, 1]
    elif order == 3:
        integrationPoints = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
        weights = [5/9, 8/9, 5/9]
    return integrationPoints, weights

# def AssembleMatrixFromMPI(matrix: np.array, numElements: int, comm, size, rank, isLast=False) -> np.array:
#     if rank > 0:
#         comm.send(matrix, dest=0, tag=41) # tag is arbitrary
#         if isLast:
#             exit() # only keeping the first processor after sending the collected data to it.
#     else:
#         for i in range(1, size):
#             if i >= numElements:
#                 break
#             tmp = comm.recv(source=i, tag=41)
#             matrix += tmp
#         return matrix

def getGlobalMatrices(GaussQuadOrder: int, DRM_elements: dict[int, list[int]], DRM_sideNodeLabels: list[int], nodes:  dict[int, list[float]], 
    materialFileName='ElementMaterial.txt', isHomogeneous=False, jobName=None, materialName=None) -> tuple[np.array, np.array, np.array]:
    """ Returns K, M, C matrices for DRM elements """
    nDOF = 3*len(DRM_sideNodeLabels)
    drmK = np.zeros([nDOF, nDOF])
    drmM = np.zeros([nDOF, nDOF])
    drmC = np.zeros([nDOF, nDOF])
    if isHomogeneous:
        density, youngsModulus, poissonsRatio, alpha, beta = getMaterialPropertiesFromInputFile(materialName, fileName=jobName+'.inp')
        # NOTE: The variable D below is usually denoted as C, but we also use C to denote the damping matrix.
        # To avoid the confusion, here we use D to denote the stiffness matrix for isotropic material.
        D = getStiffnessMatrix(youngsModulus, poissonsRatio)
    else: # heterogeneous
        # NOTE: To avoid keeping converting letter case, create lines with all lowercase letters from the material file first and reuse them
        with open(materialFileName, 'r') as f:
            lines = f.readlines()
        materialFileLowerLines = [line.lower() for line in lines] # For case insensitive search
    integrationPoints, weights = get3DGaussQuadPoints(GaussQuadOrder)
    # NOTE: Although it's possible to use MPI to distribute the computation of each element, the improvement is not significant (since this is not the major bottleneck) and makes the code more complicated.
    # eleLabels = list(DRM_elements.keys())
    # MPI_enabled, comm, size, rank = get_MPI_data()
    # if MPI_enabled:
    #     numElements = len(eleLabels)
    #     if rank > numElements:
    #         exit() # The requested number of processors exceeds the number of elements (although this is very unlikely to happen)
    #     numDistributedElements = int(numElements/size)
    #     if rank+1 == size: # The last thread
    #         eleLabels = eleLabels[rank*numDistributedElements:numElements]
    #     else:
    #         eleLabels = eleLabels[rank*numDistributedElements:(rank+1)*numDistributedElements]
    # for eleLabel in eleLabels:
    #     nodeLabels = DRM_elements[eleLabel]
    for eleLabel, nodeLabels in DRM_elements.items():
        if not isHomogeneous: # heterogeneous
            materialName = 'MATERIAL-'+str(eleLabel)
            density, youngsModulus, poissonsRatio, alpha, beta = getMaterialPropertiesFromInputFile(materialName, lowerLines=materialFileLowerLines)
            D = getStiffnessMatrix(youngsModulus, poissonsRatio)
        for i, xi in enumerate(integrationPoints):
            wXi = weights[i]
            for j, eta in enumerate(integrationPoints):
                wEta = weights[j]
                for k, zeta in enumerate(integrationPoints):
                    wZeta = weights[k]
                    N, dN = getShapeFunction(xi, eta, zeta)
                    J = getJacobian(np.array([nodes[label] for label in nodeLabels]), dN)
                    eleK = getElementStiffnessMatrix(D, J, dN, wXi, wEta, wZeta)
                    eleM = getElementMassMatrix(density, N, J, wXi, wEta, wZeta)
                    eleC = alpha*eleK + beta*eleM
                    nodeIndices = [DRM_sideNodeLabels.index(label) for label in nodeLabels]
                    globalDOF = [label for nodeIndex in nodeIndices for label in range(3*nodeIndex, 3*nodeIndex+3)]
                    # NOTE: Although this approach utilizes Advanced indexing (not Basic indexing, thus creates a copy every time), 
                    # it's still faster than a double for loop with single element indexing (one kind of Basic indexing).
                    drmK[np.ix_(globalDOF, globalDOF)] += eleK
                    drmM[np.ix_(globalDOF, globalDOF)] += eleM
                    drmC[np.ix_(globalDOF, globalDOF)] += eleC
    # if MPI_enabled:
    #     drmK = AssembleMatrixFromMPI(drmK, numElements, comm, size, rank)
    #     drmM = AssembleMatrixFromMPI(drmM, numElements, comm, size, rank)
    #     drmC = AssembleMatrixFromMPI(drmC, numElements, comm, size, rank, isLast=True)
    return drmK, drmM, drmC