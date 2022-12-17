def getInstanceName(jobName: str, partName: str) -> str | None:
    ''' NOTE: This function assume there is only one instance created from the given part '''
    with open(jobName+'.inp', 'r') as f:
        lines = f.readlines()
    lowerLines = [line.lower() for line in lines] # For case insensitive search
    for index, line in enumerate(lowerLines):
        if line.startswith('*instance') and line.endswith('part=%s\n'%partName.lower()):
            columns = [x.strip() for x in lines[index].split(',')]
            for column in columns:
                if column.startswith('name='):
                    instanceName = column[5:]
                    return instanceName
    return None