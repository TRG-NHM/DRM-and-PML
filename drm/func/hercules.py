"""
hercules.py includes tools related to Hercules.
"""

def getInputData(inputFilePath: str) -> dict:
    ''' getInputData returns a dict with all Hercules'
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
