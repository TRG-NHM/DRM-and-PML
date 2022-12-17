from getLine import getLineIndex, getNextKeywordLine

def getLabelsInSet(jobName: str, setName: str, setType:str, partName=None) -> list[int]:
    ''' getLabelsInSet returns a list containing node/element labels in a node/element set.
    NOTE: setType can be 'element' or 'node' '''
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    # NOTE: if partName is given, it should be faster and more accurate when searching through the inp file. 
    # If partName is not given, search through the whole inp file.
    partLine = 0
    if partName is not None:
        partLine = getLineIndex(lowerLines, '*Part, name=%s\n'%partName, isConversionNeeded=False)
    target = {'element': '*Elset, elset=%s\n'%setName, 
        'node': '*Nset, nset=%s\n'%setName}
    setLine = getLineIndex(lowerLines, target[setType], startLine=partLine, isConversionNeeded=False)
    if setLine is not None:
        nextKeywordLine = getNextKeywordLine(lowerLines, setLine+1)
        labelLines = lines[setLine+1:nextKeywordLine]
        # NOTE: iterable unpacking cannot be used in comprehension
        labels = [int(label) for line in labelLines for label in line.rstrip(',\n').split(',')]
    else: # target is not in the input file
        target = target[setType].rstrip(',\n') + ', generate\n'
        setLine = getLineIndex(lowerLines, target, startLine=partLine, isConversionNeeded=False)
        nextKeywordLine = getNextKeywordLine(lowerLines, setLine+1)
        labelLines = lines[setLine+1:nextKeywordLine]
        labelRanges = [[int(label) for label in line.rstrip(',\n').split(',')] for line in labelLines]
        labels = [label for lRange in labelRanges for label in range(lRange[0], lRange[1]+1, lRange[2])]
    return labels