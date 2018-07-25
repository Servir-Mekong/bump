# VIIRS packge

from __future__ import division, print_function

import ftplib
import datetime
import numpy as np
from osgeo import gdal
from scipy import ndimage


import core
import env

bumper = env.environment()


class viirs(core.raster):

    def __init__(self):
        core.raster.__init__(self,'viirs')
        
        return
        
        
    def read(self,infile):
        
        out = self._copy()

        tree = '//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/'
        field = 'SurfReflect_{0}{1}_1'
        base = 'HDF5:"{0}":{1}{2}'

        m = [i for i in range(12) if i not in [0,6,9]]
        i = [i for i in range(1,4)]
        bands = [m,i]

        res = ['1km','500m']
        mode = ['M','I']

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',1)))
        out.metadata = band.GetMetadata()
        cloudQA = self._extractBits(band.ReadAsArray(),2,3)
        hiresCloudQA = ndimage.zoom(cloudQA,2,order=0)
        band = None

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',2)))
        shadowQA = self._extractBits(band.ReadAsArray(),3,3)
        hiresShadowQA = ndimage.zoom(shadowQA,2,order=0)

        # qa = (cloudQA>0)&(shadowQA<1)
        mask = ~(hiresCloudQA>0)&(hiresShadowQA<1)

        east,west = float(out.metadata['EastBoundingCoord']), float(out.metadata['WestBoundingCoord'])
        north,south = float(out.metadata['NorthBoundingCoord']), float(out.metadata['SouthBoundingCoord'])

        out.extent = [west,south,east,north]

        databands = {'mask':mask}

        bandNames = ['mask']

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
                
                bName = '{0}{1}'.format(mode[i],bands[i][j])
                databands[bName] = data.astype(np.int16)
                bandNames.append(bName)

                band = None
                data = None
                
        out.bands = databands
        out.bandNames = bandNames
        
        out.updateMask()

        coords = {}


        out.crs = {'init':'epsg:6974'}
        out.proj = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
        
        coords['lon'],coords['lat'] = self._geoGrid(out.extent,out.bands['I1'].shape,out.proj,wgsBounds=False)

        out.coords = coords

        out.gt = self._getGt(north,west,500,out.proj)

        date = '{0}{1}{2}'.format(out.metadata['RangeBeginningDate'],out.metadata['RangeBeginningTime'],' UTC')

        out.coords['date'] = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')

        return out
    
    def fetch(self,date,h,v,outdir='./'):
        """Function to download VIIRS NRT data for specified time and tile

        Args:
            date (datetime.datetime): Datetime object specifying which date the data of interest was acquired.
            h (int): horizontal tile grid to fetch
            v (int): vertical tile grid to fetch
            outdir (str, optional): out directory to dump retrieved data to
            default = './' or current working directory
            default = None
    
        Returns:
            None
        """
    
        if outdir[-1] != '/':
            outdir = outdir+'/'
    
        usr,_,pswrd = bumper.acct.hosts['https://urs.earthdata.nasa.gov']
    
        basename = 'VNP09GA_NRT.A{0}{1:03d}.h{2:02d}v{3:02d}.001.h5'
    
        yr = date.year
        dt = (date-datetime.datetime(yr,1,1)).days + 1
    
        url = 'nrt3.modaps.eosdis.nasa.gov'
        directory = '/allData/5000/VNP09GA_NRT/Recent/'
    
        ftp = ftplib.FTP(url)
        ftp.login(usr,pswrd)
    
        ftp.cwd(directory)
    
        filename = basename.format(yr,dt,h,v)
        with open(outdir + filename, 'wb') as f:
           ftp.retrbinary('RETR ' + filename, f.write)
    
        ftp.close()
    
        return outdir+filename
        
