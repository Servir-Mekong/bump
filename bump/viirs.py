# VIIRS packge

import datetime
import numpy as np
from osgeo import gdal
from scipy import ndimage


import core
import env

bumper = env.environment()


class viirs(core.raster):

    def __init__(self,infile):
        core.raster.__init__(self,infile,'viirs')

        tree = '//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/'
        field = 'SurfReflect_{0}{1}_1'
        base = 'HDF5:"{0}":{1}{2}'

        m = [i for i in range(12) if i not in [0,6,9]]
        i = [i for i in range(1,4)]
        bands = [m,i]

        res = ['1km','500m']
        mode = ['M','I']

        band = gdal.Open(base.format(self.src,tree.format('1km'),field.format('QF',1)))
        self.metadata = band.GetMetadata()
        cloudQA = self._extractBits(band.ReadAsArray(),2,3)
        hiresCloudQA = ndimage.zoom(cloudQA,2,order=0)
        band = None

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',2)))
        shadowQA = self._extractBits(band.ReadAsArray(),3,3)
        hiresShadowQA = ndimage.zoom(shadowQA,2,order=0)

        # qa = (cloudQA>0)&(shadowQA<1)
        mask = (hiresCloudQA>0)&(hiresShadowQA<1)

        east,west = float(self.metadata['EastBoundingCoord']), float(self.metadata['WestBoundingCoord'])
        north,south = float(self.metadata['NorthBoundingCoord']), float(self.metadata['SouthBoundingCoord'])

        self.extent = [west,south,east,north]

        databands = {'QA':mask}

        bandNames = ['QA']

        for i in range(2):
            for j in range(len(bands[i])):

                subdataset = base.format(infile,tree.format(res[i]),field.format(mode[i],bands[i][j]))

                band = gdal.Open(subdataset)
                if i == 0:
                    data = ndimage.zoom(band.ReadAsArray(),2,order=0)
                else:
                    data = band.ReadAsArray()

                data = np.ma.masked_where(data<0,data)
                data = np.ma.masked_where(data>10000,data)
                data = np.ma.masked_where(mask!=0,data)
                bName = '{0}{1}'.format(mode[i],bands[i][j])
                databands[bName] = data.astype(np.int16)
                bandNames.append(bName)

                band = None
                data = None

        self.bands = databands
        self.bandNames = bandNames

        coords = {}


        self.crs = {'init':'epsg:6974'}
        self.proj = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
        
        coords['Lon'],coords['Lat'] = self._geoGrid(self.extent,self.bands['I1'].shape,self.proj,wgsBounds=False)

        self.coords = coords

        self.gt = self._getGt(north,west,self.proj)
#
#        transform = {}
#        transform['Mtransform'] = Affine.from_gdal(*mgt)
#        transform['Itransform'] = Affine.from_gdal(*igt)
#        self.transform = transform

        date = '{0}{1}{2}'.format(self.metadata['RangeBeginningDate'],self.metadata['RangeBeginningTime'],' UTC')

        self.date = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')

        return
