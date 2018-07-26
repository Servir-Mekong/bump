# main environment package for bump

from __future__ import division, print_function

import ee
import netrc

class environment(object):

    def __init__(self):
        """Initialize the environment."""
        # Initialize the Earth Engine object, using the authentication credentials.
        ee.Initialize()
        
        self.acct = netrc.netrc('../mycredentials.netrc')
        
        self.crs = 'epsg:4326'
        
        # set dates
        self.startYear = 2018
        self.endYear = 2018
        self.startJulian = 10
        self.endJulian = 50
        
        self.startDate='2018-07-22'
        self.endDate='2018-07-25'
        
        self.north = 17.500
        self.south = 15.000
        self.east = 105.000
        self.west = 102.500
        
        # define study area
        self.studyArea = [[self.west,self.south],[self.east,self.south],
        				   [self.east,self.north],[self.west,self.north],
   					   [self.west,self.south]]
        
        
        self.viirs = True
        self.atms = False
        self.landsat = False
        self.sentinel1 = False
        self.sentinel2 = False
        
        return

