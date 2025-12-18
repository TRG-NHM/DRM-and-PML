import utm

def planeCreator(origin: list[float], dimensions: list[float], 
        elementSize: float, includePreAndPostLines: bool = True, **kwargs):
    '''
    Create strings for Hercules input file to generate plane outputs.

    Args:
        origin (list[float]): [latitude, longitude] of the top-center point used in local (Abaqus) model.
        dimensions (list[float]): [EW (m), NS (m), Depth (m)] of the local (Abaqus) model.
        elementSize (float): Size of the elements used in the local (Abaqus) model.
        includePreAndPostLines (bool): Whether to include pre and post lines in the output.

    Returns:
        str: Formatted string for Hercules input file.
    '''
    # Check the element size
    if any(dim%elementSize != 0 for dim in dimensions):
        raise ValueError("Dimensions must be multiples of elementSize.")
    # Convert latitude and longitude to UTM
    utm_easting, utm_northing, zone_number, zone_letter = utm.from_latlon(origin[0], origin[1])
    southwest = (utm_easting - dimensions[0]/2, utm_northing - dimensions[1]/2)
    lat_lon_sw = utm.to_latlon(southwest[0], southwest[1], zone_number, zone_letter)
    numPlanes = int(dimensions[2] // elementSize) + 1
    output = ['number_output_planes = ' + str(numPlanes), 'output_planes ='] if includePreAndPostLines else []
    for i in range(numPlanes):
        # NOTE: x in Hercules is NS, y is EW
        output.append('   '.join(['  %.6f'%lat_lon_sw[0], '%.6f'%lat_lon_sw[1], 
            '%7d'%(i*elementSize), 
            '%10d'%elementSize, '%4d'%(dimensions[1]//elementSize+1),
            '%10d'%elementSize, '%4d'%(dimensions[0]//elementSize+1),
            '   0', '   0']))
    if includePreAndPostLines:
        output += ['# ---------   ---------   -------   ----------   ----   ----------   ----   ----   ----',
                   '#   x|lat       y-lon        z          dx        nx        dy        ny    strk    dp ']
    return '\n'.join(output)

if __name__ == "__main__":
    # ===== Example 1 =====
    # kwargs = {
    #     "origin": [41.031450, 28.784960],
    #     "dimensions": [80.0, 80.0, 20.0],
    #     "elementSize": 2
    # }
    # print(planeCreator(**kwargs))
    # ===== Example 2 =====
    kwargs = {
        "origin": [41.031450, 28.784960],
        'includePreAndPostLines': True
    }
    numElements = [40, 40, 10] # number of elements in each direction
    for elementSize in [16, 8, 4, 2]:
        kwargs["elementSize"] = elementSize
        kwargs["dimensions"] = [num*elementSize for num in numElements]
        print(planeCreator(**kwargs))