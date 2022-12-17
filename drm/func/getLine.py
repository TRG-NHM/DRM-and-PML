def getLineIndex(lines: list[str], target: str, startLine=0, isConversionNeeded=True) -> int | None:
    ''' To make it foolproof, isConversionNeeded is set to True. For better 
    performance, you may set it to be False. '''
    if isConversionNeeded: # NOTE: make the search case insensitive due to different input file writing preferences
        lowerLines = [line.lower() for line in lines[startLine:]]
    else:
        lowerLines = lines[startLine:]
    try:
        return startLine + lowerLines.index(target.lower())
    except ValueError: # .index() cannot find the target
        return None

def getNextKeywordLine(lines: list[str], startLine=0) -> int | None:
    for i, line in enumerate(lines[startLine:]):
        if line.startswith('*'):
            return startLine + i
    return None

def getNextLineIndexStartsWith(lines: list[str], heading='*', startLine=0, isConversionNeeded=True) -> int | None:
    if isConversionNeeded and heading != '*': # NOTE: This is used to make the search case insensitive due to different input file writing preferences
        lowerLines = [line.lower() for line in lines[startLine:]]
    else:
        lowerLines = lines[startLine:]
    lowerHeading = heading.lower()
    for lineIndex, line in enumerate(lowerLines):
        if line.startswith(lowerHeading):
            return startLine+lineIndex
    return None

def getAllLinesIndicesStartsWith(lines: list[str], heading='*', startLine=0, isConversionNeeded=True) -> list[int]:
    ''' NOTE: If there is no line starts with the given heading, the returned list will be empty. '''
    if isConversionNeeded and heading != '*': # NOTE: This is used to make the search case insensitive due to different input file writing preferences
        lines[startLine:] = [line.lower() for line in lines[startLine:]]
    return [startLine+lineIndex for lineIndex, line in enumerate(lines[startLine:]) if line.startswith(heading.lower())]
