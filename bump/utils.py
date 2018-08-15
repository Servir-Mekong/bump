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



# gee helpers

import ee
import math

UPPER_LEFT = 0
LOWER_LEFT = 1
LOWER_RIGHT = 2
UPPER_RIGHT = 3
PI = lambda: ee.Number(math.pi)
MAX_SATELLITE_ZENITH = 7.5

def line_from_coords(coordinates, fromIndex, toIndex):
    return ee.Geometry.LineString(ee.List([
        coordinates.get(fromIndex),
        coordinates.get(toIndex)]))


def line(start, end):
    return ee.Geometry.LineString(ee.List([start, end]))


def degToRad(deg):
    return deg.multiply(PI().divide(180))


def value(list, index):
    return ee.Number(list.get(index))


def radToDeg(rad):
    return rad.multiply(180).divide(PI())


def where(condition, trueValue, falseValue):
    trueMasked = trueValue.mask(condition)
    falseMasked = falseValue.mask(invertMask(condition))
    return trueMasked.unmask(falseMasked)


def invertMask(mask):
    return mask.multiply(-1).add(1)


def x(point):
    return ee.Number(ee.List(point).get(0))


def y(point):
    return ee.Number(ee.List(point).get(1))


def determine_footprint(image):
    footprint = ee.Geometry(image.get('system:footprint'))
    bounds = ee.List(footprint.bounds().coordinates().get(0))
    coords = footprint.coordinates()

    xs = coords.map(lambda item: x(item))
    ys = coords.map(lambda item: y(item))

    def findCorner(targetValue, values):
        diff = values.map(lambda value: ee.Number(value).subtract(targetValue).abs())
        minValue = diff.reduce(ee.Reducer.min())
        idx = diff.indexOf(minValue)
        return coords.get(idx)

    lowerLeft = findCorner(x(bounds.get(0)), xs)
    lowerRight = findCorner(y(bounds.get(1)), ys)
    upperRight = findCorner(x(bounds.get(2)), xs)
    upperLeft = findCorner(y(bounds.get(3)), ys)

    return ee.List([upperLeft, lowerLeft, lowerRight, upperRight, upperLeft])


def replace_bands(image, bands):
    result = image
    for band in bands:
        result = result.addBands(band, overwrite=True)
    return result



# Execute the main level program if run as standalone
if __name__ == "__main__":
    test()
