# ATMS package

from __future__ import division, print_function

import datetime
import numpy as np
from osgeo import gdal
#from scipy import ndimage

import core
import env


bumper = env.environment()


class atms(core.raster):
    def __init__(self):
        core.raster.__init__(self,'atms')
        
        df = np.load('atms_mlc_coeffs.npy')
        self.calFactors = dict(df.item())
        
        return
    
    
    def read(self,infile,ingeo):
        
        out = self._copy()

        trees = ['//All_Data/ATMS-SDR_All/',
                 '//All_Data/ATMS-SDR-GEO_All/'
                ]

        fields = [['GainCalibration','BrightnessTemperature',],
                  ['Latitude','Longitude','SatelliteZenithAngle']
                 ]

        base = 'HDF5:"{0}":{1}{2}'    

        band = gdal.Open(base.format(infile,trees[0],fields[0][0]))
        out.gainCal = band.ReadAsArray()
        
        band = None
        
        band = gdal.Open(base.format(infile,trees[0],fields[0][1]))
        data = band.ReadAsArray()
        out.metadata = band.GetMetadata()
        
        btrBands = ['C1','C2','C3','C4','C5','C16','C17','C18','C19','C20','C21','C22']
        toKeep = [0,1,2,3,4,15,16,17,18,19,20,21]
        
        databands = {'mask':np.ones([data.shape[0],data.shape[1]])}

        bandNames = ['mask']
        
        for b in range(len(toKeep)):
            calData = data[:,:,toKeep[b]]

            BTr = np.zeros_like(calData).astype(np.uint16)
            for i in range(BTr.shape[0]):
                BTr[i,:] = (calData[i,:].astype(int) * out.gainCal[i,toKeep[b]]) * 0.0001
        
            databands[btrBands[b]] = BTr
            bandNames.append(btrBands[b])
            
        band = None
        data = None
        
        out.bands = databands
        out.bandNames = bandNames
        
        out.updateMask()
        
        out.nativeCRS = {'init':'epsg:4326'}
        out.proj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

        geoBands = ['lat','lon','view']
        
        coords = {}
        
        for c in range(len(geoBands)):
            band = gdal.Open(base.format(ingeo,trees[1],fields[1][c]))
            coords[geoBands[c]] = band.ReadAsArray()

            band = None
            data = None
            
        out.coords = coords
                
        out.bands['mask'] = out.coords['view']<50 | out.bands['mask'].astype(np.bool)
        
        out.updateMask()
        
        out.gt = None
        
        fileComps = infile.split('_')
        timeComps = fileComps[2][1:] + fileComps[3][1:] + 'UTC'
        
        out.date = datetime.datetime.strptime(timeComps,'%Y%m%d%H%M%S%f%Z')
        
        return out
    
    def _mlc(self):
        
        keys = self.bands.keys()
#        keys = list(self.calFactors.keys())
        
        Btr = np.moveaxis(np.array(self.bands.values()),0,-1)
        
        cal_mlc = np.zeros_like(Btr)
        
        for c in range(len(keys)):
            chnl = keys[c]
            if chnl not in "mask":
                tbi = self.calFactors[chnl]['Tbi']
                preds = self.calFactors[chnl]['predictors']
                popt = self.calFactors[chnl]['popt']
                for i in range(cal_mlc.shape[0]):
                    factor = np.zeros(cal_mlc.shape[1])
                    for j in range(popt.size):
                        Tbij = self.calFactors[keys[j]]['Tbij']
                        factor = factor + (popt[j]*(Btr[i,:,preds[j]]-Tbij))
                    cal_mlc[i,:] = tbi + factor
                
        return cal_mlc