import utm
import pyproj
import numpy as np

def getDistanceBetweenTwoCoordinatesPyproj(point1, point2):
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

def getDistanceBetweenTwoCoordinates(point1, point2):
    ''' The format of point1 and point 2 should be (latitude, longitude)'''
    utm_easting1, utm_northing1, zone_number1, zone_letter1 = utm.from_latlon(point1[0], point1[1])
    utm_easting2, utm_northing2, zone_number2, zone_letter2 = utm.from_latlon(point2[0], point2[1])
    if zone_number1 != zone_number2 or zone_letter1 != zone_letter2:
        return getDistanceBetweenTwoCoordinatesPyproj(point1, point2)
    else:
        x_dis = utm_northing2 - utm_northing1
        y_dis = utm_easting2 - utm_easting1
        distance = np.sqrt(x_dis**2 + y_dis**2)
        return distance, x_dis, y_dis