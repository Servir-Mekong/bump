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
        self.startYear = 2015
        self.endYear = 2015
        self.startJulian = 10
        self.endJulian = 50
        
        self.north = 19.999
        self.south = 18.552
        self.east = 105.806
        self.west = 103.876
        
        # define study area
        self.studyArea = [[self.west,self.south],[self.east,self.south],
        				   [self.east,self.north],[self.west,self.north],
   					   [self.west,self.south]]
        return

