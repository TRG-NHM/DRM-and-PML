from getLine import getLineIndex, getNextLineIndexStartsWith

def getDataFromHerculesInputFile(inputFilePath: str) -> dict:
    ''' getDataFromHerculesInputFile returns a dict with all 
    parameters as the dict keys and values as the dict values. '''
    with open(inputFilePath) as f:
        lines = f.readlines()
        lines = [line.rstrip().split('=') for line in lines if (line[0] != '#' and line.rstrip() != '')]
    inputData = {}
    heldKey = None
    for line in lines:
        if len(line) == 2 and line[1] != '':
            inputData[line[0].rstrip()] = line[1].lstrip()
            heldKey = None
        elif len(line) == 2 and line[1] == '':
            heldKey = line[0].rstrip()
            inputData[line[0].rstrip()] = []
        elif len(line) == 1:
            inputData[heldKey].append(line[0].split())
        else:
            raise ValueError('Unexpected input file format.')
    return inputData

# def getMaterialPropertiesFromInputFile(jobName: str, materialName: str) -> tuple[float, float, float]:
#     with open(jobName+'.inp', 'r') as f:
#         lines = f.readlines()
#     lowerLines = [line.lower() for line in lines] # For case insensitive search
#     materialLine = getLineIndex(lowerLines, '*Material, name=%s\n'%materialName, isConversionNeeded=False)
#     densityLine = getLineIndex(lowerLines, '*Density\n', startLine=materialLine, isConversionNeeded=False)
#     density = float(lines[densityLine+1].split(',')[0])
#     elasticLine = getLineIndex(lowerLines, '*Elastic\n', startLine=materialLine, isConversionNeeded=False)
#     elasticProperties = lines[elasticLine+1].split(',')
#     youngsModulus = float(elasticProperties[0])
#     poissonsRatio = float(elasticProperties[1])
#     return density, youngsModulus, poissonsRatio

def getMaterialPropertiesFromInputFile(materialName: str, fileName=None, lowerLines=None) -> tuple[float, float, float, float, float]:
    ''' getMaterialPropertiesFromInputFile returns density, youngsModulus, poissonsRatio, alpha, and beta of the given material. 
    If the file is an Abaqus input file, the full name of it should be passed (i.e., 'example.inp' for fileName). 
    If lowerLines is provided, fileName is not necessary. '''
    if fileName is not None and lowerLines is None:
        with open(fileName, 'r') as f:
            lines = f.readlines()
        lowerLines = [line.lower() for line in lines] # For case insensitive search
    materialLine = getLineIndex(lowerLines, '*Material, name=%s\n'%materialName, isConversionNeeded=False)
    densityLine = getLineIndex(lowerLines, '*Density\n', startLine=materialLine, isConversionNeeded=False)
    density = float(lowerLines[densityLine+1].split(',')[0])
    elasticLine = getLineIndex(lowerLines, '*Elastic\n', startLine=materialLine, isConversionNeeded=False)
    elasticProperties = lowerLines[elasticLine+1].split(',')
    youngsModulus = float(elasticProperties[0])
    poissonsRatio = float(elasticProperties[1])
    nextMaterialLine = getNextLineIndexStartsWith(lowerLines, '*Material,', startLine=materialLine+1, isConversionNeeded=False)
    if nextMaterialLine is None:
        nextMaterialLine = -1
    dampingLine = getNextLineIndexStartsWith(lowerLines[:nextMaterialLine], '*Damping,', startLine=materialLine+1, isConversionNeeded=False)
    if dampingLine is None:
        alpha = 0
        beta = 0
    else:
        dampingProperties = lowerLines[dampingLine].rstrip('\n').split(',')[1:]
        properties = {}
        for prop in dampingProperties:
            x = prop.split('=')
            properties[x[0].strip()] = float(x[1])
        alpha = properties['alpha']
        beta = properties['beta']
    return density, youngsModulus, poissonsRatio, alpha, beta