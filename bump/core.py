# base functionality package

from __future__ import division, print_function

#import glob
import copy
import math
from itertools import groupby

import numpy as np
import xarray as xr
from osgeo import gdal,osr
from pyproj import Proj, transform
from scipy import interpolate 

import env
import utils
#import viirs


bumper = env.environment()


class grid(object):
    def __init__(self,resolution):

        if '4326' in bumper.crs:
            midPoint = [np.mean([bumper.north,bumper.south]),
                        np.mean([bumper.east,bumper.west])]
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
        self.dims = self.xx.shape

        return
    
    def mapRaster(self,raster,interpMethod='linear'):
        
        rasterLons = raster.coords['lon']
        rasterLats = raster.coords['lat']
        
        
        xSelect = (rasterLons>bumper.west) & (rasterLons<bumper.east)
        ySelect = (rasterLats>bumper.south) & (rasterLats<bumper.north)
        spatialSelect = ySelect & xSelect
        
        idx = np.where(spatialSelect == True)
        
        # Format geolocation coordinates for gridding
        pts = np.zeros((idx[0].size,2))
        pts[:,0] = rasterLons[idx].ravel()
        pts[:,1] = rasterLats[idx].ravel()
        
        bNames = raster.bands.keys()
        
        qualityBands = np.zeros([self.dims[0],self.dims[1],len(bNames)])
                
        out = raster._copy()
        
        for i in range(len(bNames)):
            if 'mask' in bNames[i]:
                iMethod = 'nearest'
            else:
                iMethod = interpMethod
            
            # Regrid data to common grid 
            out.bands[bNames[i]] = interpolate.griddata(pts, 
                                     raster.bands[bNames[i]][idx].ravel(),
                                     (self.xx,self.yy), method=iMethod)
            
            qualityBands[:,:,i] = (out.bands[bNames[i]] < 16000) & \
                                  (out.bands[bNames[i]] > -1)
            if i == 0:
                interpMask = np.isnan(out.bands[bNames[i]])
                
        qualityMask = np.min(qualityBands,axis=2).astype(np.bool)
        out.bands['mask'] = out.bands['mask'] & ~interpMask & qualityMask
        out.updateMask()
            
        out.coords['Lon'],out.coords['Lat'] = self.xx,self.yy
        out.extent = bumper.west,bumper.south,bumper.east,bumper.north
        
        out.gt = out._getGt(bumper.north,bumper.west,self.nominalResolution)

        return out


class raster(object):
    def __init__(self,sensor):
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
    

    def _geoGrid(self,extent,dims,nativeProj=None,wgsBounds=True):

        west, south, east, north = extent

        gcsProj = Proj(init=bumper.crs)
        native = Proj(nativeProj)

        if wgsBounds == False:
            native = Proj(nativeProj)

            llx,lly = transform(gcsProj,native,west,south)
            urx,ury = transform(gcsProj,native,east,north)

            yCoords = np.linspace(lly,ury,dims[0],endpoint=False)[::-1]
            xCoords = np.linspace(llx,urx,dims[1],endpoint=False)

            xx,yy = np.meshgrid(xCoords,yCoords)

            lons,lats = transform(native,gcsProj,xx,yy)

        elif wgsBounds == True:
            llx,lly = west,south
            urx,ury = east,north

            yCoords = np.linspace(lly,ury,dims[0],endpoint=False)[::-1]
            xCoords = np.linspace(llx,urx,dims[1],endpoint=False)

            lons,lats = np.meshgrid(xCoords,yCoords)

        else:
            raise ValueError('argument "wgsBounds" needs to be binary')

        return lons,lats
    
    
    def _getGt(self,north,west,gridSize,projStr=None):
        if projStr:
            outProj = Proj(projStr)
            inProj = Proj(init=bumper.crs)

            ulx,uly = transform(inProj,outProj,west,north)
            
        else:
            ulx,uly = bumper.west,bumper.north
            
        if type(gridSize) == list:
            pxSize = gridSize[0]
            pySize = gridSize[1]
        else:
            pxSize,pySize = gridSize, gridSize

        #(originX, pixelWidth, 0, originY, 0, pixelHeight)

        gt = (ulx,pxSize,0,uly,0,-pySize)

        return gt
    
    
    def updateMask(self):
        bNames = [i for i in self.bands.keys() if i != 'mask']
        mask = self.bands['mask']
                
        for i in bNames:
            data = self.bands[i]
            self.bands[i] = np.ma.masked_where(mask==0,data)
            
        return
    
    
    def unmask(self,value=None):
        bNames = [i for i in self.bands.keys() if i != 'mask']
        mask = self.bands['mask']
        
        out = self._copy()

        if value:
            for i in bNames:
                out.bands[i][np.where(mask==0)] = value
                
        else:
            for i in bNames:
                out.bands[i] = self.bands[i].data
                
        return out
    
    
    def writeGeotiff(self,fileName):
        
            return
        
        
    def normalizedDifference(self,band1,band2):
        nd = (self.bands[band1] - self.bands[band2]) / \
             (self.bands[band1] + self.bands[band2])
        
        self.bands['nd'] = nd
        
        clone = self._copy()

        for n in clone.bandNames:
            if n != 'mask':
                del clone.bands[n]
            
        clone.bandNames = ['mask','nd']
        
        return clone
    


class collection(object):
    def __init__(self,rasterList,gr,tile=False):
        gt = rasterList[0].gt
        sensor =rasterList[0].sensor
        bandList = rasterList[0].bandNames
        collDates = []
        for i in rasterList:
            collDates.append([i,i.coords['date']])
        
        outArr = []
        outDates = []
        masks = []
        for k,g in groupby(collDates,lambda x: x[1]):
            this = list(g)
            outDates.append(this[0][1])
            newList = []
            ms = []
            for i in this:
                newList.append(i[0].bands.values())
                ms.append(i[0].bands['mask'])
                
            outArr.append(np.mean(newList,axis=0))
            masks.append(np.min(ms))
        
        x = gr.xx[0,:]
        y = gr.yy[:,0]
        
        data = {'x':x,'y':y}
        
        for i in range(len(bandList)):
            tData = np.zeros([gr.dims[0],gr.dims[1],len(outArr)])
            for j in range(len(outArr)):
                if rasterList[j].bandNames != bandList:
                    raise AttributeError('All band names for rasters must be \
                                         the same for collection')
                        
                else:
                    temp = outArr[j][i,:,:]
                    temp[np.where(masks[j]==0)] = 0
                    if sensor == 'viirs':
                        tData[:,:,j] = np.flipud(temp)
                    else:
                        tData[:,:,j] = temp
                    
                data[bandList[i]] = (['y','x','time'],tData)
        
            
        ds = xr.Dataset(data,
                        coords={'lon': (['y', 'x'], gr.xx),
                                'lat': (['y', 'x'], gr.yy),
                                'time': outDates},
                        attrs={'gt':gt,}
                        )
        
        if tile:
            ds = ds.chunk({'x': 100, 'y': 100})
        
        self.data = ds.sortby('time').where(ds!=0)
        
        return 
    
    def writeNetCDF(self,fileName):   
        
        return
    
    def writeGeotiffs(self,folder,prefix='bump_out'):
        
        nFiles,y,x = self.data.dims.values()
        dates = self.data.time.values
        
        drv = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(bumper.crs.split(':')[1]))
                
        
        for i in range(nFiles):
            t = str(dates[i]).split('.')[0].replace(':','')
            print(t)
            tValues = self.data.isel(time=i).to_array().values
            
            bands = tValues.shape[0]
            
            outDs = drv.Create(prefix + '_' + t + '.tif',
                               y,
                               x,
                               bands,
                               gdal.GDT_Int16
                               )
            
            for b in range(bands):
                band = outDs.GetRasterBand(b+1)
                band.WriteArray(tValues[b,:,:])
                band.SetNoDataValue(0)
                band = None
            
            outDs.SetGeoTransform(self.data.attrs['gt'])
            
            outDs.SetProjection(srs.ExportToWkt())
            
            outDs.FlushCache()
        
        return
