import pandas as pd
import os

def print_material_statistics(filePaths: list[str]|str) -> None:
    ''' Print basic statistics of the material properties stored in the given file(s) '''
    if isinstance(filePaths, str):
        filePaths = [filePaths]
    for filePath in filePaths:
        if 'df' not in locals():
            df = pd.read_csv(filePath, index_col=0)
        else:
            df_new = pd.read_csv(filePath, index_col=0)
            df = pd.concat([df, df_new], axis=0)
    # if there's column named 'poissonsRatio_corrected', remove any rows with negative Poisson's Ratio
    if 'poissonsRatio_corrected' in df.columns:
        df = df[df['poissonsRatio_corrected'] >= 0]
    df.describe().to_csv('material_statistics.csv')
    return

if __name__ == '__main__':
    # Change the fileFolder to the folder where your material CSV files are located
    fileFolder = '/Volumes/CORSAIR/Abaqus/Cases/1Hz_v3/site4_soil_only_8m_element'
    fileNames = ['rawMaterialDRM.csv', 'rawMaterialPML.csv']
    filePaths = [os.path.join(fileFolder, fileName) for fileName in fileNames]
    print_material_statistics(filePaths)