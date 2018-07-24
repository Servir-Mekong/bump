#******************************************************************************
# FILE: utils.py
# AUTHOR: Kel Markert
# EMAIL: kel.markert@nasa.gov
# ORGANIZATION: NASA-SERVIR, UAH/ESSC
# MODIFIED BY: n/a
# CREATION DATE: 25 Jan. 201y
# LAST MOD DATE: 03 Apr. 2017
# PURPOSE: This script converts a length unit and coverts it from meters to
#          decimal degrees or vise versa
# DEPENDENCIES: numpy
#******************************************************************************

# import dependency
import math

def meters2dd(inPt,scale=30):
    """Function to convert meters to decimal degrees based on the approximation
    given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system

    Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.

    Keywords:
        scale (int): Resolution of the raster value to covert into decimal
        degrees, must be in meters.

    Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees

    """

    lat = inPt[0] # get latitude value

    radLat = math.radians(lat) # convert degree latitude to radians

    a = 6378137 # radius of Earth in meters

    ba = 0.99664719 # constant of b/a

    ss = math.atan(ba*math.tan(radLat)) # calculate the reduced latitude

    # factor to convert meters to decimal degrees for X axis
    xfct = (math.pi/180)*a*math.cos(ss)

    # factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*math.cos(2*radLat)+1.175*math.cos(4*radLat)-
              0.0023*math.cos(6*radLat))

    # get decimal degree resolution
    ydd = scale / yfct
    xdd = scale / xfct

    # return list of converted resolution values
    return [ydd,xdd]

def dd2meters(inPt,scale=0.1):
    """Function to convert decimal degrees to meters based on the approximation
        given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system

        Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.

        Keywords:
        scale (int): Resolution of the raster value to covert into meters,
        must be in decimal degrees.

        Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees

        """

    lat = inPt[0] # get latitude value

    radLat = math.radians(lat) # convert degree latitude to radians

    a = 6378137 # radius of Earth in meters

    ba = 0.99664719 # constant of b/a

    ss = math.atan(ba*math.tan(radLat)) # calculate the reduced latitude

    # factor to convert meters to decimal degrees for X axis
    xfct = (math.pi/180)*a*math.cos(ss)

    # factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*math.cos(2*radLat)+1.175*math.cos(4*radLat)-
            0.0023*math.cos(6*radLat))

    # get meter resolution
    y_meters = scale * yfct
    x_meters = scale * xfct

    # return list of converted resolution values
    return [y_meters,x_meters]

def test():
    # main level program for testing
    inPt = [-1,0] # -1 degrees latitude, 0 degrees longitue
    outMeters = meters2dd(inPt,100)
    outDD = dd2meters(inPt,1)

    # print results
    print 'Output decimal degrees: {0}\nOutput meters: {1}'.format(outMeters,outDD)

    return

# Execute the main level program if run as standalone
if __name__ == "__main__":
    test()
