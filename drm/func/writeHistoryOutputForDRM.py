import pandas as pd
from getNodes import getNodeCoordinates
from getLabelsInSet import getLabelsInSet
from searchHDF5 import getConvertedGridPointsForAbaqusModel, getInterpolatedHistoryDataForGridPoints, getNearestHistoryDataForGridPoints

def writeDisplacementHistoryForDRM(jobName: str, partName: str, dispHistoryFileName='DispHistory.csv', 
    isCoordinateConverted: bool = False, matchMethod: str = 'nearest', 
    **kwargs) -> pd.DataFrame:
    ''' # NOTE: mpirun works here since we use `getInterpolatedHistoryDataForGridPoints()` function '''
    nodes = getNodeCoordinates(jobName, partName)
    DRM_interiorNodeLabels = getLabelsInSet(jobName, setName='inDRM', setType='node')
    DRM_exteriorNodeLabels = getLabelsInSet(jobName, setName='outDRM', setType='node')
    DRM_sideNodeLabels = DRM_interiorNodeLabels + DRM_exteriorNodeLabels
    DRM_nodes = [nodes[label] for label in DRM_sideNodeLabels]
    kwargs.update({'nodes': DRM_nodes, 'pointLabelList': DRM_sideNodeLabels,
        'isCoordinateConverted': isCoordinateConverted, 'gridPointsInMeter': True})
    kwargs['gridPoints'] = getConvertedGridPointsForAbaqusModel(**kwargs)
    if matchMethod == 'interpolated':
        df = getInterpolatedHistoryDataForGridPoints(**kwargs)
    elif matchMethod == 'nearest':
        df = getNearestHistoryDataForGridPoints(**kwargs)
    if isCoordinateConverted:
        # NOTE: If the directions used in Abaqus model are EW, NS, and pointing up from the earth,
        # then isCoordinateConverted can be set to True to correct the differences between the directions used in Hercules and Abaqus.
        # NOTE: It might be confusing, but no matter whether the coordinate is converted, 
        # the coordinates (x, y, z) are always in Hercules' notation. On the other 
        # hand, the displacements will be rearranged to fit Abaqus' notation.
        tmp = df['u'].copy()
        df['u'] = df['v']
        df['v'] = tmp
        df.loc[df['w']!=0, 'w'] = -df.loc[df['w']!=0, 'w']
    # Round the floating number to 6 decimal places
    for col in ['time', 'x', 'y', 'z', 'u', 'v', 'w']:
        df[col] = df[col].round(6)
    if dispHistoryFileName.endswith('.csv'):
        df.to_csv(dispHistoryFileName, index=False)
    else:
        df.to_hdf(dispHistoryFileName, key='dispHistory', mode='w', format='table', data_columns=True, index=False, complevel=9, complib="blosc:lz4")
    return