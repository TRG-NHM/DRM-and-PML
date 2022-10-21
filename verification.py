import matplotlib.pyplot as plt
import pandas as pd
import os
import io

def getStationList(folderPath: str, neglect_subfolder=True) -> list[str]:
    stationList = []
    for dirpath, dirNames, fileNames in os.walk(folderPath):
        for fileName in fileNames:
            prefix, suffix = os.path.splitext(fileName)
            if prefix == 'station' or suffix == '.csv':
                stationList.append(fileName)
        if neglect_subfolder:
            dirNames[:] = [] # clear subdirectories. This act would neglect subdirectories after searching at current directory
    return stationList

def getFileWithoutUnnecessaryHeading(filePath: str) -> io.StringIO:
    with open(filePath, 'r') as f:
        lines = f.readlines()
        lines[0] = lines[0].lstrip('#')
    # with open(filePath, 'w') as f:
    #     f.writelines(lines)
    return io.StringIO(''.join(lines))

def getTimeHistoryFromStationAndAbaqusResult(stationFolder: str, AbaqusResultFolder: str) -> list[pd.DataFrame]:
    # NOTE: The sequence would be [1, 0]. So the location sequence would be ['Center', 'Top Center']
    stationList = getStationList(stationFolder)
    AbaqusResultList = getStationList(AbaqusResultFolder)
    stations = []
    AbaqusResults = []
    for stationFileName in stationList:
        prefix, suffix = os.path.splitext(stationFileName)
        stationFilePath = os.path.join(stationFolder, stationFileName)
        stationFile = getFileWithoutUnnecessaryHeading(stationFilePath)
        df = pd.read_csv(stationFile, delim_whitespace=True, index_col='Time(s)')
        stations.append(df)
        AbaqusResultFileName = [fileName for fileName in AbaqusResultList if suffix+'.csv' in fileName][0]
        df = pd.read_csv(os.path.join(AbaqusResultFolder, AbaqusResultFileName), index_col='Time(s)')
        AbaqusResults.append(df)
    return stations, AbaqusResults

def plotAndSaveResults(locationNames: list[str], stations: list[pd.DataFrame], AbaqusResults: list[pd.DataFrame], outputFolder='') -> None:
    columnSeries = {'Displacement (m)': {'$u$': 'X|(m)', '$v$': 'Y-(m)', '$w$': 'Z.(m)'},
        'Velocity (m/s)': {'$\dot{u}$': 'X|(m/s)', '$\dot{v}$': 'Y-(m/s)', '$\dot{w}$': 'Z.(m/s)'},
        'Accerleration (m/s$^2$)': {'$\ddot{u}$': 'X|(m/s2)', '$\ddot{v}$': 'Y-(m/s2)', '$\ddot{w}$': 'Z.(m/s2)'}}
    plt.rcParams.update({'font.size': 10, 'font.weight': 'regular'})
    for resultNum, locationName in enumerate(locationNames):
        fig, axes = plt.subplots(3, 3, sharex='col', sharey='row')
        fig.set_size_inches([8.5, 8.5])
        times = AbaqusResults[resultNum].index
        for i, (quantity, columns) in enumerate(columnSeries.items()):
            for j, (title, column) in enumerate(columns.items()):
                axes[i, j].plot(times, AbaqusResults[resultNum][column], '-', linewidth=0.5, label='Abaqus')
                axes[i, j].plot(times, stations[resultNum].loc[times, column], '-', linewidth=0.5, label='Hercules')
                axes[i, j].legend(frameon=False)
                axes[i, j].set(title=title)
                axes[-1, j].set(xlabel='Time (s)')
            axes[i, 0].set(ylabel=quantity)
        fig.savefig(os.path.join(outputFolder, locationName+'.pdf'))

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    stationFolder = './Verification/stations'
    AbaqusResultFolder = './Verification/AbaqusResults'
    stations, AbaqusResults = getTimeHistoryFromStationAndAbaqusResult(stationFolder, AbaqusResultFolder)
    locationNames = ['Center', 'Top Center']
    plotAndSaveResults(locationNames, stations, AbaqusResults, outputFolder='./Verification')