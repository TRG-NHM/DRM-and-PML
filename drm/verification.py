import matplotlib.pyplot as plt
import pandas as pd
import os
from func.getHistoryOutputForDRM import getFileWithoutUnnecessaryHeading

def getStationList(folderPath: str, neglect_subfolder=True) -> list[str]:
    stationList = []
    for dirpath, dirNames, fileNames in os.walk(folderPath):
        for fileName in fileNames:
            prefix, suffix = os.path.splitext(fileName)
            if prefix == 'station' or suffix == '.rpt':
                stationList.append(fileName)
        if neglect_subfolder:
            dirNames[:] = [] # clear subdirectories. This act would neglect subdirectories after searching at current directory
    return stationList

def getTimeHistoryFromStationAndAbaqusResult(stationFolder: str, AbaqusResultFolder: str, truncateTime=None) -> list[pd.DataFrame]:
    # NOTE: The sequence would be [1, 0]. So the location sequence would be ['Center', 'Top Center']
    stationList = getStationList(stationFolder)
    AbaqusResultList = getStationList(AbaqusResultFolder)
    stations = {}
    AbaqusResults = {}
    for stationFileName in stationList:
        stationNum = int(stationFileName.split('.')[1])
        stationFile = getFileWithoutUnnecessaryHeading(os.path.join(stationFolder, stationFileName))
        df = pd.read_csv(stationFile, delim_whitespace=True, index_col='Time(s)')
        stations[stationNum] = df
    for AbaqusResultFileName in AbaqusResultList:
        nodeLabel = int(AbaqusResultFileName.split('.')[1])
        df = pd.read_csv(os.path.join(AbaqusResultFolder, AbaqusResultFileName), sep='\s{2,}', engine='python', skiprows=1, index_col='X')
        df.index.name = 'Time(s)'
        df.columns = ['X|(m)', 'Y-(m)', 'Z.(m)', 'X|(m/s)', 'Y-(m/s)', 'Z.(m/s)', 'X|(m/s2)', 'Y-(m/s2)', 'Z.(m/s2)']
        if truncateTime is not None:
            df = df.loc[truncateTime[0]:truncateTime[1]]
            df.index = df.index - truncateTime[0]
            df.index = df.index.to_series().round(4) # Pandas Float64Index doesn't define round method
        AbaqusResults[nodeLabel] = df
    return stations, AbaqusResults

def plotAndSaveResults(locations: dict[str, dict[str, int]], stations: dict[int, pd.DataFrame], AbaqusResults: dict[int, pd.DataFrame], outputFolder='', isCoordinateConverted=False) -> None:
    columnSeries = {'Displacement (m)': {'$u$': 'X|(m)', '$v$': 'Y-(m)', '$w$': 'Z.(m)'},
        'Velocity (m/s)': {'$\dot{u}$': 'X|(m/s)', '$\dot{v}$': 'Y-(m/s)', '$\dot{w}$': 'Z.(m/s)'},
        'Acceleration (m/s$^2$)': {'$\ddot{u}$': 'X|(m/s2)', '$\ddot{v}$': 'Y-(m/s2)', '$\ddot{w}$': 'Z.(m/s2)'}}
    plt.rcParams.update({'font.size': 10, 'font.weight': 'regular'})
    for locationName, resultNum in locations.items():
        stationNum = resultNum['station']
        nodeLabel = resultNum['nodeLabel']
        fig, axes = plt.subplots(3, 3, sharex='col', sharey='row')
        fig.set_size_inches([8.5, 8.5])
        times = AbaqusResults[nodeLabel].index
        if isCoordinateConverted:
            AbaqusResults[nodeLabel].columns = ['Y-(m)', 'X|(m)', 'Z.(m)', 'Y-(m/s)', 'X|(m/s)', 'Z.(m/s)', 'Y-(m/s2)', 'X|(m/s2)', 'Z.(m/s2)']
            for column in ['Z.(m)', 'Z.(m/s)', 'Z.(m/s2)']:
                AbaqusResults[nodeLabel][column] = -AbaqusResults[nodeLabel][column]
        for i, (quantity, columns) in enumerate(columnSeries.items()):
            for j, (title, column) in enumerate(columns.items()):
                axes[i, j].plot(times, AbaqusResults[nodeLabel][column], '-', linewidth=0.5, label='Abaqus')
                axes[i, j].plot(times, stations[stationNum].loc[times, column], '-', linewidth=0.5, label='Hercules')
                axes[i, j].legend(frameon=False)
                axes[i, j].set(title=title)
                axes[-1, j].set(xlabel='Time (s)')
            axes[i, 0].set(ylabel=quantity)
        fig.savefig(os.path.join(outputFolder, locationName+'.pdf'))

if __name__ == '__main__':
    # /// CASE 1: prospectusModel
    # stationFolder = './Verification/stations'
    # AbaqusResultFolder = './Verification/AbaqusResults'
    # locations = {'Soil box top center': {'station': 0, 'nodeLabel': 8302}, 'Soil box center': {'station': 1, 'nodeLabel': 8922}}
    # stations, AbaqusResults = getTimeHistoryFromStationAndAbaqusResult(stationFolder, AbaqusResultFolder)
    # plotAndSaveResults(locations, stations, AbaqusResults, outputFolder='./Verification')
    # /// CASE 2: IstanbulModel
    stationFolder = '../Cases/Abaqus Steel Building Model from Bulent/Verification/stations'
    AbaqusResultFolder = '../Cases/Abaqus Steel Building Model from Bulent/Verification/AbaqusResults'
    locations = {'Istanbul soil box top center': {'station': 0, 'nodeLabel': 29406}}
    truncateTime = [1.0, 40.99]
    stations, AbaqusResults = getTimeHistoryFromStationAndAbaqusResult(stationFolder, AbaqusResultFolder, truncateTime=truncateTime)
    plotAndSaveResults(locations, stations, AbaqusResults, outputFolder='./Verification', isCoordinateConverted=True)