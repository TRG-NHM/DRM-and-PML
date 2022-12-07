import csv
import pandas as pd
from math import sin, cos, asin, atan2, sqrt, degrees as rad2deg, radians as deg2rad
from getDistance import getDistanceBetweenTwoCoordinates
from getLine import getLineIndex, getNextKeywordLine
from getInstanceName import getInstanceName
from getLabelsInSet import getLabelsInSet

def getInstanceTranslation(inpLines: list[str], partName: str, instanceName=None) -> list[list[float]]:
    # NOTE: Not consider the rotation at this moment
    if instanceName is None:
        instanceName = partName+'-1'
    target = '*Instance, name=%s, part=%s\n'%(instanceName, partName)
    index = getLineIndex(inpLines, target)
    line = inpLines[index+1]
    if line.lower() == '*end instance\n':
        return [0, 0, 0]
    else:
        translation = [float(x) for x in line.split(',')]
        return translation

def getNodeCoordinates(jobName: str, partName: str) -> dict[int, list[float]]:
    ''' getNodeCoordinates returns all the node coordinates on the given part in the format of {nodeLabel: [x, y, z]}. '''
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    translation = getInstanceTranslation(lowerLines, partName, instanceName=getInstanceName(jobName, partName))
    partLine = getLineIndex(lowerLines, '*Part, name=%s\n'%partName, isConversionNeeded=False)
    nodeStartLine = getLineIndex(lowerLines, '*Node\n', startLine=partLine+1, isConversionNeeded=False)
    nodeEndLine = getNextKeywordLine(lowerLines, startLine=nodeStartLine+1)
    nodes = {}
    for line in lines[nodeStartLine+1:nodeEndLine]:
        node = line.rstrip(',\n').split(',')
        nodes[int(node[0].strip())] = [float(x)+translation[i] for i, x in enumerate(node[1:])]
    return nodes

def getDestination(origin: tuple, bearing: float, d: float) -> tuple[float]:
    ''' 
    origin = (lat, lon)
    bearing is the initial bearing, which sometimes referred to as forward azimuth, in radian from north (counted clockwise)
    d is distance traveled in meters
    '''
    R = 6371E3 # Earth's mean radius in meters
    lat1, lon1 = origin
    phi1 = deg2rad(lat1)
    lambda1 = deg2rad(lon1)
    phi2 = asin( sin(phi1)*cos(d/R) + cos(phi1)*sin(d/R)*cos(bearing) )
    lambda2 = lambda1 + atan2( sin(bearing)*sin(d/R)*cos(phi1), cos(d/R) - sin(phi1)*sin(phi2) )
    lat2 = rad2deg(phi2)
    lon2 = rad2deg(lambda2)
    return (lat2, lon2)

def getCorrectedCoordinatesForHercules(nodes: dict[int | str, list[float]], 
    HerculesModelOrigin=None, AbaqusModelOrigin=None, result='LatLon', 
    isCoordinateConverted=False) -> dict[int | str, list[float]]:
    if result=='Meter' and HerculesModelOrigin is not None and AbaqusModelOrigin is not None:
        distance, x_dis, y_dis = getDistanceBetweenTwoCoordinates(HerculesModelOrigin, AbaqusModelOrigin)
    else:
        x_dis, y_dis = 0, 0
    # NOTE: After this for loop, the order of coordinates for nodes becomes [north, east, point down to the earth] in meters
    for label, node in nodes.items():
        if isCoordinateConverted:
            nodes[label] = [node[1]+x_dis, node[0]+y_dis, -node[2]]
        else:
            nodes[label] = [node[0]+x_dis, node[1]+y_dis, node[2]]
        if node[2] == 0:
            nodes[label][2] = 0
    if result == 'LatLon' and AbaqusModelOrigin is not None:
        for label, node in nodes.items():
            d = sqrt(node[0]**2 + node[1]**2)
            bearing = atan2(node[0], node[1])
            lat, lon = getDestination(AbaqusModelOrigin, bearing, d)
            nodes[label] = [lat, lon, node[2]]
    return nodes

# def getNodesAndWriteNodeTable(fileName: str, jobName: str, partName: str, coordinateCorrection=False, 
#     HerculesModelOrigin=None, AbaqusModelOrigin=None, isCoordinateConverted=False) -> None:
#     with open(fileName, 'r') as f:
#         lines = f.readlines()
#     nodeCoordinates = getNodeCoordinates(jobName, partName)
#     if coordinateCorrection:
#         nodeCoordinates = getCorrectedCoordinatesForHercules(nodeCoordinates, HerculesModelOrigin, AbaqusModelOrigin, isCoordinateConverted=isCoordinateConverted)
#     startLine = 0
#     index = 0
#     numNode = 1
#     nodeTable = []
#     target = '*Cload, amplitude=AMP-Node'+str(numNode)+'_x\n'
#     while index is not None:
#         nodeLabel = int(lines[index+1].split(',')[0].split('.')[1])
#         nodeTable.append([numNode, nodeLabel]+nodeCoordinates[nodeLabel])
#         numNode = numNode + 1
#         startLine = index + 1
#         target = '*Cload, amplitude=AMP-Node'+str(numNode)+'_x\n'
#         index = getLineIndex(lines, target, startLine)
#     with open('nodeTable.csv', 'w', newline='') as f:
#         writer = csv.writer(f)
#         writer.writerow(['nodeNum', 'nodelabel', 'x', 'y', 'z'])
#         writer.writerows(nodeTable)
#     return

def writeNodeTable(jobName: str, partName: str, nodeSetName='DRM', coordinateCorrection=False, 
    AbaqusModelOrigin=None, isCoordinateConverted=False) -> None:
    nodeLabelsOnDRM = getLabelsInSet(jobName, nodeSetName, 'node', partName=partName)
    nodeCoordinates = getNodeCoordinates(jobName, partName)
    if coordinateCorrection:
        nodeCoordinates = getCorrectedCoordinatesForHercules(nodeCoordinates, AbaqusModelOrigin=AbaqusModelOrigin, isCoordinateConverted=isCoordinateConverted)
    nodeTable = []
    for i, nodeLabel in enumerate(nodeLabelsOnDRM):
        nodeTable.append([i+1, nodeLabel]+nodeCoordinates[nodeLabel])
    with open('nodeTable.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['nodeNum', 'nodelabel', 'x', 'y', 'z'])
        writer.writerows(nodeTable)
    return

def fitOldNodeTable(newTableFileName, oldTableFileName):
    newTable = pd.read_csv(newTableFileName)
    oldTable = pd.read_csv(oldTableFileName)
    for i in range(len(newTable)):
        oldRow = oldTable[(oldTable['x'] == newTable.iloc[i]['x']) & (oldTable['y'] == newTable.iloc[i]['y']) & (oldTable['z'] == newTable.iloc[i]['z'])]
        # newTable.iloc[i]['nodeNum'] = oldRow['nodeNum'] # NOTE: this way won't change the value and trigger SettingWithCopyWarning
        newTable['nodeNum'][i] = oldRow['nodeNum'] # NOTE: although this way still triggers SettingWithCopyWarning, the value will be changed
    newTable.to_csv(newTableFileName, index=False)
    return

if __name__ == '__main__':
    jobName = 'Istanbul_model_complete'
    partName = 'Part-Soil'
    newTableFileName = 'nodeTable.csv'
    oldTableFileName = '../Abaqus Steel Building Model from Bulent/nodeTable.csv'
    AbaqusModelOrigin = (41.0318, 28.9417)
    # HerculesModelOrigin=(40.58711947, 28.58207271)
    # writeNodeTable(jobName, partName, coordinateCorrection=True, 
    #     AbaqusModelOrigin=AbaqusModelOrigin, isCoordinateConverted=True)
    fitOldNodeTable(newTableFileName, oldTableFileName)