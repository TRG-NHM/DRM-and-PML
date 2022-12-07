import pyproj
import numpy as np

def getDistanceBetweenTwoCoordinates(point1, point2):
    ''' The format of point1 and point 2 should be (latitude, longitude)'''
    geod = pyproj.Geod(ellps='WGS84')
    # NOTE: distance is in meter, azimuths are in degree (forward azimuth is counted clockwise from the north)
    # NOTE 2: the input order of geod.inv is longitude, latitude
    forwardAzimuth, backAzimuth, distance = geod.inv(point1[1], point1[0], point2[1], point2[0])
    # NOTE: x_dis is the distance component along North-South direction (North is positive)
    x_dis = distance*np.cos(np.deg2rad(forwardAzimuth))
    # NOTE: y_dis is the distance component along East-West direction (East is positive)
    y_dis = distance*np.sin(np.deg2rad(forwardAzimuth))
    return distance, x_dis, y_dis