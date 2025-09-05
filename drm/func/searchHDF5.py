import pandas as pd
# import pyproj
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import Point, Polygon
from getDistance import getDistanceBetweenTwoCoordinates
from get_MPI_data import get_MPI_data
from time import time
from tqdm import tqdm

def getDistanceFromPlaneOrigin(targetPoint, planeData):
    origin = list(planeData[:2])
    return getDistanceBetweenTwoCoordinates(origin, targetPoint)

def getProperFilter(planesData: pd.DataFrame, gridPoints: list[list[float]], 
    outputFormat: dict|str = str):
    # NOTE: gridPoints should be in meters
    if type(gridPoints[0]) is not list: # If there is only one point
        gridPoints = [gridPoints]
    planeData = planesData.iloc[0]
    x_dis = [p[0] for p in gridPoints]
    y_dis = [p[1] for p in gridPoints]
    depth = [p[2] for p in gridPoints]
    zList = planesData['z']
    stepAlongStrike = planeData.iloc[3] # Distance between nodes in X direction (dstrk)
    stepDownDip = planeData.iloc[5] # Distance between nodes in Y direction
    filters = {'x': {'upper': max(x_dis)+stepAlongStrike, 'lower': min(x_dis)-stepAlongStrike}, 
        'y': {'upper': max(y_dis)+stepDownDip, 'lower': min(y_dis)-stepDownDip}, 
        'z': {'upper': zList[zList>=max(depth)].min(), 'lower': zList[zList<=min(depth)].max()}}
    if outputFormat is dict:
        return filters
    elif outputFormat is str:
        whereStatement = ""
        criteria = []
        for key, value in filters.items():
            if len(value) == 1 and list(value.keys()) == ['equal']:
                criteria.append(key+" = "+str(value['equal']))
            elif all([v is not None for v in value.values()]):
                criteria.append(key+" >= "+str(value['lower'])+" and "+key+" <= "+str(value['upper']))
            elif value['lower'] is not None:
                criteria.append(key+" >= "+str(value['lower']))
            elif value['upper'] is not None:
                criteria.append(key+" <= "+str(value['upper']))
        whereStatement += ' and '.join(criteria)
        return whereStatement
    else:
        raise ValueError('The output format should be either dict or str.')

def getSortedDisp(xyz, uvw):
    dimension = len(xyz[0])
    xyz_comp = [[] for i in range(dimension)]
    for i in range(dimension):
        xyz_comp[i] = list(set([cord[i] for cord in xyz]))
    [cord.sort() for cord in xyz_comp]
    sortedDisp = []
    if dimension == 3:
        for x in xyz_comp[0]:
            for y in xyz_comp[1]:
                for z in xyz_comp[2]:
                    idx = xyz.index([x, y, z])
                    sortedDisp.append(uvw[idx])
    elif dimension == 2:
        for x in xyz_comp[0]:
            for y in xyz_comp[1]:
                idx = xyz.index([x, y])
                sortedDisp.append(uvw[idx])
    length_xyz = [len(comp) for comp in xyz_comp]
    disp = np.array(sortedDisp).reshape([*length_xyz, 3])
    return xyz_comp, disp

def getInterpolatedDisp(targetPointInMeter, xyz, uvw):
    # NOTE: If any dimension of disp used in RegularGridInterpolator is 1, any return value will be nan.
    # NOTE: To correctly apply trilinear interpolation, use `RegularGridInterpolator` instead of `LinearNDInterpolator`.
    xyz_comp, disp = getSortedDisp(xyz, uvw)
    interpolatedFunc = RegularGridInterpolator(xyz_comp, disp)
    return interpolatedFunc(targetPointInMeter)[0]

def getCollectedResultFromMPI(listData, numDomains, comm, size, rank):
    if rank > 0:
        comm.send(listData, dest=0, tag=41) # tag is arbitrary
        exit() # only keeping the first processor after sending the collected data to it.
    else:
        for i in range(1, size):
            if i >= numDomains:
                break
            tmp = comm.recv(source=i, tag=41)
            listData += tmp
        return listData

# TODO: domainSurfaceCorners is no longer written in the HDF5 file.
# Maybe compute it back or find another way to define the domain surface.
def isPointWithinDomain(targetPoint, dbPath):
    planesData = pd.read_hdf(dbPath, 'planesData')
    domainSurfaceCorners = pd.read_hdf(dbPath, 'domainSurfaceCorners')
    domainSurface = Polygon(domainSurfaceCorners.to_numpy())
    surfacePoint = Point(targetPoint[1], targetPoint[0])
    if domainSurface.contains(surfacePoint) and targetPoint[2]<planesData['z'].max():
        return True
    else:
        return False

def getClosestPlaneKeys(planesData, targetPoint):
    zList = planesData['z']
    depth = targetPoint[-1]
    return ['plane%i'%i for i in planesData.index if planesData.loc[i]['z'] in [zList[zList>=depth].min(), zList[zList<=depth].max()]]

def getAllPossiblePlaneKeys(planesData, gridPoints):
    zList = planesData['z']
    depthList = [targetPoint[-1] for targetPoint in gridPoints]
    maxZ = zList[zList>=max(depthList)].min()
    minZ = zList[zList<=min(depthList)].max()
    possibleZ = zList[(zList>=minZ) & (zList<=maxZ)]
    return ['plane%i'%i for i in planesData.index if planesData.loc[i]['z'] in possibleZ.values]

def getCheckedAndConvertedGridPoints(dbPath, gridPoints, gridPointsInMeter=False):
    planesData = pd.read_hdf(dbPath, 'planesData')
    convertedGridPoints = []
    for targetPoint in gridPoints:
        # TODO: Capability to check whether a point is within the domain if gridPointsInMeter=True
        # if not isPointWithinDomain(targetPoint, dbPath) and gridPointsInMeter is False:
        #     return 'The inquiry point is not located within the domain.'
        if gridPointsInMeter:
            targetPointInMeter = targetPoint
        else:
            distance, x_dis, y_dis = getDistanceFromPlaneOrigin(targetPoint, planesData.loc[0])
            targetPointInMeter = [x_dis, y_dis, targetPoint[-1]]
        convertedGridPoints.append(targetPointInMeter)
    return convertedGridPoints

def getInterpolatedHistoryDataForGridPoints(gridPoints: list[list[float]]|list[float], 
        dbPath: str, pointLabelList: list[int]|None = None, gridPointsInMeter: bool = False,
        planeIndices: list[int]|None = None, **kwargs):
    tic = time()
    if type(gridPoints[0]) is not list: # If there is only one point
        gridPoints = [gridPoints]
    numDomains = len(gridPoints)
    if pointLabelList is None:
        pointLabelList = list(range(1, numDomains+1))
    MPI_enabled, comm, size, rank = get_MPI_data()
    if not MPI_enabled or (MPI_enabled and rank == 0):
        # Construct a minimum Pandas DataFrame (`df`) that contains the whole interested domain (`gridPoints`)
        planesData = pd.read_hdf(dbPath, 'planesData')
        if planeIndices is not None:
            planesData = planesData.loc[planeIndices]
        planeKeys = getAllPossiblePlaneKeys(planesData, gridPoints)
        gridPoints = getCheckedAndConvertedGridPoints(dbPath, gridPoints, gridPointsInMeter)
        filters = getProperFilter(planesData, gridPoints)
        df = pd.read_hdf(dbPath, key=planeKeys[0], where=filters) # get filtered dataset from the HDF5 file
        for planeKey in planeKeys[1:]:
            newData = pd.read_hdf(dbPath, key=planeKey, where=filters)
            df = pd.concat([df, newData])
        print('Preliminaries: %.2f secs' % (time() - tic))
    else:
        df = None
        planesData = None
    if MPI_enabled:
        gridPoints = comm.bcast(gridPoints, root=0)
        df = comm.bcast(df, root=0)
        planesData = comm.bcast(planesData, root=0)
        if rank >= numDomains:
            exit() # The requested number of processors exceeds the number of domains
        numDistributedDomains = int(numDomains/size)
        if rank+1 == size: # The last thread
            gridPoints = gridPoints[rank*numDistributedDomains:numDomains]
            pointLabelList = pointLabelList[rank*numDistributedDomains:numDomains]
        else:
            gridPoints = gridPoints[rank*numDistributedDomains:(rank+1)*numDistributedDomains]
            pointLabelList = pointLabelList[rank*numDistributedDomains:(rank+1)*numDistributedDomains]
    interpolatedData = []
    for i, targetPointInMeter in tqdm(enumerate(gridPoints), total=len(gridPoints), 
        position=rank, desc=f'Processor {rank} - Interpolating history data'):
        newInterpolatedData, columns = getInterpolatedHistoryData(targetPointInMeter, df, planesData, returnType=list)
        interpolatedData.extend([[pointLabelList[i]]+line for line in newInterpolatedData])
    if MPI_enabled:
        interpolatedData = getCollectedResultFromMPI(interpolatedData, numDomains, comm, size, rank)
        print('Combining results from all MPI processes: %.2f secs' % (time() - tic))
    df = pd.DataFrame(interpolatedData, columns=columns.insert(0, 'pointLabel'))
    return df

def getInterpolatedHistoryData(targetPointInMeter, df, planesData, returnType=pd.DataFrame):
    interpolatedData = []
    filters = getProperFilter(planesData, targetPointInMeter)
    df = df.query(filters)
    for timeStep in df['timeStep'].unique():
        df_at_t = df.loc[df['timeStep']==timeStep]
        uvw = df_at_t[['u', 'v', 'w']].values.tolist()
        t = df_at_t['time'].iloc[0]
        # NOTE: If there is only one element in any dimension in the data set, `getInterpolatedDisp` will not work.
        if len(planesData['z']) == 1 or targetPointInMeter[-1] in list(planesData['z']):
            xy = df_at_t[['x', 'y']].values.tolist()
            interpolatedData.append([timeStep, t, *targetPointInMeter, *getInterpolatedDisp(targetPointInMeter[:2], xy, uvw)])
        else:
            xyz = df_at_t[['x', 'y', 'z']].values.tolist()
            interpolatedData.append([timeStep, t, *targetPointInMeter, *getInterpolatedDisp(targetPointInMeter, xyz, uvw)])
    if returnType is list:
        return interpolatedData, df.columns
    else: # returnType is pd.DataFrame
        df = pd.DataFrame(interpolatedData, columns=df.columns)
        return df

def alignGridPoints(gridPoints: pd.DataFrame, plane: pd.DataFrame, 
        allowance: float = 0.1, **kwargs) -> None:
    '''
    Align the x and y coordinates of gridPoints to the closest grid point in the plane.
    This function modifies the gridPoints DataFrame in place.
    '''
    for comp in ['x', 'y']:
        completedPointLabels = []
        for pointLabel, point in tqdm(gridPoints.iterrows(), desc='Aligning grid points in %s-axis' % comp):
            if pointLabel in completedPointLabels:
                continue
            closestPoints = plane[(plane[comp] >= point[comp] - allowance) & 
                                  (plane[comp] <= point[comp] + allowance)]
            if closestPoints.empty:
                raise ValueError("No data found for pointLabel %i. " \
                    "Consider increasing the allowance." % (pointLabel))
            closestCompValue = closestPoints[comp].unique()
            if closestCompValue.size > 1:
                raise ValueError("Multiple closest %s components found for pointLabel %i. " \
                    "Consider decreasing the allowance." % (comp, pointLabel))
            subset = gridPoints.loc[(gridPoints[comp] == point[comp])]
            gridPoints.loc[subset.index, comp] = closestCompValue
            completedPointLabels.extend(subset.index.tolist())
    return None

def getNearestHistoryDataForGridPoints(gridPoints: list[list[float]]|list[float], 
    dbPath: str, pointLabelList: list[int]|None = None, planeIndices: list[int]|None = None, 
    **kwargs):
    if type(gridPoints[0]) is not list: # If there is only one point
        gridPoints = [gridPoints]
    if pointLabelList is None:
        pointLabelList = list(range(1, len(gridPoints)+1))
    gridPoints = pd.DataFrame(gridPoints, columns=['x', 'y', 'z'], index=pointLabelList)
    # NOTE: Users should ensure planes used have the same horizontal dimensions, 
    # same origin, but unique z.
    # TODO: Implement MPI
    planesData = pd.read_hdf(dbPath, 'planesData')
    if planeIndices is not None:
        planesData = planesData.loc[planeIndices]
    alignGridPoints(gridPoints, pd.read_hdf(dbPath, key='plane%i'%planesData.index[0]), **kwargs)
    df = pd.DataFrame()
    for i, planeData in planesData.iterrows():
        print("Processing plane %i: " % i)
        plane = pd.read_hdf(dbPath, key='plane%i'%i)
        plane.insert(0, 'pointLabel', 0)
        gP = gridPoints[gridPoints['z'] == planeData['z']]
        if gP.empty:
            print("No grid points found for plane %i at z = %f. Skipping..." % (i, planeData['z']))
            continue
        newDataIndices = []
        for pointLabel, point in tqdm(gP.iterrows()):
            newData = plane[(plane['x'] == point['x']) & (plane['y'] == point['y'])]
            if newData.empty:
                raise ValueError("No data found for pointLabel %i at plane %i. " \
                    "Making sure all planes have the same grid except z value." % (pointLabel, i))
            plane.loc[newData.index, 'pointLabel'] = pointLabel
            newDataIndices.extend(newData.index.to_list())
        df = pd.concat([df, plane.loc[newDataIndices]])
    df.reset_index(drop=True, inplace=True)
    return df

def getGrid(cornerPoints, depth, numNodes=3):
    gridPoints = []
    # /// For points on four walls
    for i in range(len(cornerPoints)):
        prevPoint = cornerPoints[i-1]
        point = cornerPoints[i]
        latNodes = np.linspace(prevPoint[0], point[0], numNodes)
        lonNodes = np.linspace(prevPoint[1], point[1], numNodes)
        linspacedPoints = list(zip(latNodes[1:], lonNodes[1:])) # skip the first (lat, lon) since it will be the last combination for the next iteration.
        if depth != 0:
            depthNodes = np.linspace(0, depth, numNodes)
        else:
            depthNodes = [depth]
        gridPoints += [[linspacedPoint[0], linspacedPoint[1], depthNode] for depthNode in depthNodes for linspacedPoint in linspacedPoints] 
    # /// For points on the bottom surface
    numSpaces = numNodes - 1
    pointsOnLines = []
    for i in range(2):
        primaryVector = cornerPoints[2*i+1] - cornerPoints[2*i]
        # NOTE: Excluding the first and the last points on the line since it has been included when we generate points on the four walls
        pointsOnLines.append([ cornerPoints[2*i] + primaryVector*(k/numSpaces) for k in range(1, numSpaces) ])
    for i in range(numNodes-2):
        secondaryVector = pointsOnLines[1][-(i+1)] - pointsOnLines[0][i]
        pointsOnSecondaryLine = [ pointsOnLines[0][i] + secondaryVector*(k/numSpaces) for k in range(1, numSpaces) ]
        gridPoints += [list(points)+[depth] for points in pointsOnSecondaryLine]
    return gridPoints
    
def getParallelogramVertices(cornerPoints):
    # NOTE: cornerPoints should be a list containing three points (lat, lon), and this function will find the fourth vertex on the parallelogram.
    if len(cornerPoints) != 3:
        raise ValueError("There should be only 3 points in the input variable 'cornerPoints'.")
    # NOTE: For easier vector calculation, the cornerPoints will be transformed to np.array instead of list
    cornerPoints = np.array(cornerPoints)
    vector = cornerPoints[2] - cornerPoints[1]
    cornerPoints = np.vstack([cornerPoints, cornerPoints[0] + vector])
    return cornerPoints

def visualizeGirdPoints(gridPoints): # Used only for debugging
    import matplotlib.pyplot as plt
    fig = plt.figure()
    axis = fig.add_subplot(projection='3d')
    axis.scatter(*zip(*gridPoints))
    plt.show()
    return

def getPointsFromCSV(f):
    df = pd.read_csv(f)
    return df.values.tolist()

# TODO: Adjustment if the Abaqus model origin is not same as Hercules' database (practically 
# every case) and axes are not aligned to latitude and longitude direction
def getConvertedGridPointsForAbaqusModel(dbPath: str, nodes: list[list[float]], 
        origin: tuple[float]|None = None, gridPointsInMeter: bool = False, 
        isCoordinateConverted: bool = False, planeIndices: list[int] = [0], 
        **kwargs):
    # istanbulProj = pyproj.Proj('+proj=utm +ellps=WGS84 +units=m +lat_0=%f +lon_0=%f'%(origin[0], origin[1]))
    if gridPointsInMeter is False:
        raise ValueError("Grid Point Conversion with Latitude and Longitude is not yet supported.")
    planesData = pd.read_hdf(dbPath, 'planesData')
    if origin is not None:
        distance, x_dis, y_dis = getDistanceFromPlaneOrigin(origin, planesData.loc[planeIndices[0]])
        if isCoordinateConverted:
            # NOTE: If the directions used in Abaqus model are EW, NS, and pointing up from the earth,
            # then isCoordinateConverted can be set to True to correct the differences between the directions used in Hercules and Abaqus.
            gridPoints = [[node[1]+x_dis, node[0]+y_dis, -node[2]] for node in nodes]
        else:
            gridPoints = [[node[0]+x_dis, node[1]+y_dis, node[2]] for node in nodes]
    elif isCoordinateConverted: # origin is None
        gridPoints = [[node[1], node[0], -node[2]] for node in nodes]
    return gridPoints
