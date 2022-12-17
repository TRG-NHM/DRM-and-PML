import pandas as pd
from getNodes import getNodeCoordinates
from getLabelsInSet import getLabelsInSet
from searchHDF5 import getConvertedGridPointsForAbaqusModel, getInterpolatedHistoryDataForGridPoints

def writeDisplacementHistoryForDRM(dbPath: str, jobName: str, partName: str, 
    targetOrigin: list[float], dispHistoryFileName='DispHistory.csv') -> pd.DataFrame:
    ''' # NOTE: mpirun works here since we use `getInterpolatedHistoryDataForGridPoints()` function '''
    nodes = getNodeCoordinates(jobName, partName)
    DRM_interiorNodeLabels = getLabelsInSet(jobName, setName='inDRM', setType='node')
    DRM_exteriorNodeLabels = getLabelsInSet(jobName, setName='outDRM', setType='node')
    DRM_sideNodeLabels = DRM_interiorNodeLabels + DRM_exteriorNodeLabels
    DRM_nodes = [nodes[label] for label in DRM_sideNodeLabels]
    gridPoints = getConvertedGridPointsForAbaqusModel(dbPath, DRM_nodes, origin=targetOrigin, gridPointsInMeter=True)
    df = getInterpolatedHistoryDataForGridPoints(gridPoints, dbPath, pointLabelList=DRM_sideNodeLabels, gridPointsInMeter=True)
    df.to_csv(dispHistoryFileName)
    # return df
    return