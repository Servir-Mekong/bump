# base functionality package

from __future__ import division, print_function

import glob
import copy

import numpy as np
import xarray as xr

import env
import utils


bumper = env.environment()


class grid(object):
    def __init__(self,resolution):

        if '4326' in bumper.crs:
            midPoint = [np.mean(bumper.north,bumper.south),
                        np.mean(bumper.east,bumper.west)]
            ySpacing, xSpacing = utils.meters2dd(midPoint,resolution)

        elif type(resolution) == list:
            xSpacing = resolution[0]
            ySpacing = resolution[1]

        else:
            ySpacing,xSpacing= resolution,resolution

        lons = np.arange(bumper.west,bumper.east,xSpacing)
        lats = np.arange(bumper.south,bumper.north,ySpacing)

        self.xx,self.yy = np.meshgrid(lons,lats)

        self.nominalResolution = [xSpacing,ySpacing]

        return

    def raster2grid(self,raster):

        return


class raster(object):
    def __init__(self,path,sensor,crs='epsg:4326'):
        self.src = path
        self.sensor = sensor
        self.crs = {'init':bumper.crs}

        return


    def _copy(self):
        return copy.deepcopy(self)


    def _extractBits(self,image,start,end):
        """Helper function to convert Quality Assurance band bit information to flag values

        Args:
            image (ndarray): Quality assurance image as a numpy array
            start (int): Bit position to start value conversion
            end (int): Bit position to end value conversion

        Returns:
            out (ndarray): Output quality assurance in values from bit range
        """

        pattern = 0;
        for i in range(start,end+1):
            pattern += math.pow(2, i)

        bits = image.astype(np.uint16) & int(pattern)
        out = bits >> start

        return out


class collection(object):
    def __init__(self,folder,grid):


        return
