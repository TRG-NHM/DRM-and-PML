from getLine import getLineIndex, getNextKeywordLine

def getElements(jobName: str, partName: str, elementType: str) -> dict[int, list[int]]:
    ''' getElements returns the element connectivity on the given part in the format of {elementLabel: [nodeLabels]}. '''
    # TODO: elementType could be None and collect all elements regardless their element types
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    partLine = getLineIndex(lowerLines, '*Part, name=%s\n'%partName, isConversionNeeded=False)
    elementLine = getLineIndex(lowerLines, '*Element, type=%s\n'%elementType, startLine=partLine, isConversionNeeded=False)
    nextKeywordLine = getNextKeywordLine(lowerLines, elementLine+1) # NOTE: the +1 is critical!
    elementLines = lines[elementLine+1:nextKeywordLine]
    elements = [line.rstrip(',\n').split(',') for line in elementLines]
    elements = {int(element[0]): [int(x) for x in element[1:]] for element in elements}
    return elements
