import pandas as pd
import numpy as np
import os, io

def getHistoryOutputForDRMFromDispHistoryFile(dispHistoryFileName='DispHistory.csv') -> tuple[list[float], dict[int, dict[str, np.array]]]:
    """ This function reads the displacement history and compute velocity and 
    acceleration histories. Finally, it returns a Dict that contains all 3 
    histories. """
    if dispHistoryFileName.endswith('.csv'):
        df = pd.read_csv(dispHistoryFileName)
    else:
        df = pd.read_hdf(dispHistoryFileName)
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

def getFileWithoutUnnecessaryHeading(filePath: str) -> io.StringIO:
    ''' getFileWithoutUnnecessaryHeading returns the file in string IO format so that it keeps in memory without writing/overwriting to a file. '''
    with open(filePath, 'r') as f:
        lines = f.readlines()
        lines[0] = lines[0].lstrip('#')
    # with open(filePath, 'w') as f:
    #     f.writelines(lines)
    return io.StringIO(''.join(lines))

def getHistoryOutputForDRMFromStationFiles(stationFolder: str, nodeTableFileName='nodeTable.csv', isCoordinateConverted=False, 
    nodeLabels=None, truncateTime=None, isResponseRecalculationNeeded=False)  -> tuple[list[float], dict[int, dict[str, np.array]]]:
    """ This function reads the displacement history and compute velocity and 
    acceleration histories. Finally, it returns a Dict that contains all 
    3 histories. """
    # NOTE: If the response is computed from Hercules, the directions are NS, EW, and pointing down toward the earth.
    # If this is not the direction combination used in the Abaqus model, isCoordinateConverted can be set to True to correct it.
    nodeTable = pd.read_csv(nodeTableFileName, index_col=0)
    # TODO: MPI can be implemented here
    histories = {}
    if isCoordinateConverted:
        columns = [quantity+direction for quantity in ['u', 'v', 'a'] for direction in ['y', 'x', 'z']]
    else:
        columns = [quantity+direction for quantity in ['u', 'v', 'a'] for direction in ['x', 'y', 'z']]
    for nodeNum in nodeTable.index:
        nodeLabel = int(nodeTable.loc[nodeNum, 'nodelabel'])
        if nodeLabels is not None and nodeLabel not in nodeLabels:
            continue
        stationFileName = 'station.'+str(nodeNum)
        stationFilePath = os.path.join(stationFolder, stationFileName)
        stationFile = getFileWithoutUnnecessaryHeading(stationFilePath)
        df = pd.read_csv(stationFile, delim_whitespace=True, index_col='Time(s)')
        if truncateTime is not None:
            df = df.loc[truncateTime[0]:truncateTime[1]]
            df.index = df.index - truncateTime[0]
        df.columns = columns
        # NOTE: We may need to compute the velocity and acceleration by ourselves because sometimes Hercules generate weirdly large acceleration.
        if isResponseRecalculationNeeded:
            for direction in ['x', 'y', 'z']:
                df['v'+direction] = df['u'+direction].diff()/df.index.to_series().diff()
                df['v'+direction].iloc[0] = 0
                df['a'+direction] = df['v'+direction].diff()/df.index.to_series().diff()
                df['a'+direction].iloc[0] = 0
        histories[nodeLabel] = {column: df[column].to_numpy() for column in columns}
        if isCoordinateConverted:
            for quantity in ['u', 'v', 'a']:
                histories[nodeLabel][quantity+'z'] = -histories[nodeLabel][quantity+'z']
    return df.index.to_list(), histories