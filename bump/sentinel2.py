# Sentinel-2 package

import ee
from Py6S import *
import math
import datetime
import os, sys
from utils import *
sys.path.append("/gee-atmcorr-S2/bin/")
from atmospheric import Atmospheric
import sun_angles
import view_angles

 
#import ee.mapclient

class env(object):

	def __init__(self):
		"""Initialize the environment."""

		# Initialize the Earth Engine object, using the authentication credentials.
		ee.Initialize()
		self.startDate = "2016-01-15"
		self.endDate = "2016-01-26"
		self.location = ee.Geometry.Point([-80.72,-1.34])
		countries = ee.FeatureCollection('ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw')
		self.location  = countries.filter(ee.Filter.inList('Country', ['Ecuador'])).geometry();
		
		self.metadataCloudCoverMax = 80
		self.cloudThreshold = 10
		self.hazeThresh = 200
		self.exportBands = ee.List(['blue','green','red'])
        
		self.s2BandsIn = ee.List(['QA60','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12'])
		self.s2BandsOut = ee.List(['QA60','cb','blue','green','red','re1','re2','re3','nir1','nir2','waterVapor','cirrus','swir1','swir2'])
		self.divideBands = ee.List(['cb','blue','green','red','re1','re2','re3','nir1','nir2','waterVapor','cirrus','swir1','swir2'])
		
		self.feature = 0
		
		# 9. Cloud and cloud shadow masking parameters.
		# If cloudScoreTDOM is chosen
		# cloudScoreThresh: If using the cloudScoreTDOMShift method-Threshold for cloud 
		#    masking (lower number masks more clouds.  Between 10 and 30 generally works best)
		self.cloudScoreThresh = 10;
		
		# Percentile of cloud score to pull from time series to represent a minimum for 
		# the cloud score over time for a given pixel. Reduces commission errors over 
		# cool bright surfaces. Generally between 5 and 10 works well. 0 generally is a bit noisy	
		self.cloudScorePctl = 10; 
		
		# zScoreThresh: Threshold for cloud shadow masking- lower number masks out 
		# less. Between -0.8 and -1.2 generally works well
		self.zScoreThresh = -1.15

		# shadowSumThresh: Sum of IR bands to include as shadows within TDOM and the 
		# shadow shift method (lower number masks out less)
		self.shadowSumThresh = 3000;
		
		# contractPixels: The radius of the number of pixels to contract (negative buffer) clouds and cloud shadows by. Intended to eliminate smaller cloud 
		#    patches that are likely errors (1.5 results in a -1 pixel buffer)(0.5 results in a -0 pixel buffer)
		# (1.5 or 2.5 generally is sufficient)
		self.contractPixels = 2.5; 
		
		# dilatePixels: The radius of the number of pixels to dilate (buffer) clouds 
		# and cloud shadows by. Intended to include edges of clouds/cloud shadows 
		# that are often missed (1.5 results in a 1 pixel buffer)(0.5 results in a 0 pixel buffer)
		# (2.5 or 3.5 generally is sufficient)
		self.dilatePixels = 2.5;
		
		self.calcSR = True      
		self.brdf = True
		self.QAcloudMask = False
		self.cloudMask = False
		self.shadowMask = True
		self.terrainCorrection = True


class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
	
	def TOAtoSR(self,img):
		
		#TDOMMask = img.select(['TDOMMask'])
		info = self.collectionMeta[self.env.feature]['properties']
		scene_date = datetime.datetime.utcfromtimestamp(info['system:time_start']/1000)# i.e. Python uses seconds, EE uses milliseconds
		solar_z = info['MEAN_SOLAR_ZENITH_ANGLE']
        
		geom = ee.Geometry.Point([info['centroid']['coordinates'][0],info['centroid']['coordinates'][1]])
		date = ee.Date.fromYMD(scene_date.year,scene_date.month,scene_date.day)
		
		h2o = Atmospheric.water(geom,date).getInfo()
		o3 = Atmospheric.ozone(geom,date).getInfo()
		aot = Atmospheric.aerosol(geom,date).getInfo()

		SRTM = ee.Image('CGIAR/SRTM90_V4')# Shuttle Radar Topography mission covers *most* of the Earth
		alt = SRTM.reduceRegion(reducer = ee.Reducer.mean(),geometry = geom).get('elevation').getInfo()
		
		if alt:
			km = alt/1000 # i.e. Py6S uses units of kilometers
		
		else:
			km = 0
		# Instantiate
		s = SixS()
		
		# Atmospheric constituents
		s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o,o3)
		s.aero_profile = AeroProfile.Continental
		s.aot550 = aot
		
		# Earth-Sun-satellite geometry
		s.geometry = Geometry.User()
		s.geometry.view_z = 0               # always NADIR (I think..)
		s.geometry.solar_z = solar_z        # solar zenith angle
		s.geometry.month = scene_date.month # month and day used for Earth-Sun distance
		s.geometry.day = scene_date.day     # month and day used for Earth-Sun distance
		s.altitudes.set_sensor_satellite_level()
		s.altitudes.set_target_custom_altitude(km)
		
		def spectralResponseFunction(bandname):
			"""
			Extract spectral response function for given band name
			"""
            
			bandSelect = {
				'B1':PredefinedWavelengths.S2A_MSI_01,
				'B2':PredefinedWavelengths.S2A_MSI_02,
				'B3':PredefinedWavelengths.S2A_MSI_03,
				'B4':PredefinedWavelengths.S2A_MSI_04,
				'B5':PredefinedWavelengths.S2A_MSI_05,
				'B6':PredefinedWavelengths.S2A_MSI_06,
				'B7':PredefinedWavelengths.S2A_MSI_07,
				'B8':PredefinedWavelengths.S2A_MSI_08,
				'B8A':PredefinedWavelengths.S2A_MSI_09,
				'B9':PredefinedWavelengths.S2A_MSI_10,
				'B10':PredefinedWavelengths.S2A_MSI_11,
				'B11':PredefinedWavelengths.S2A_MSI_12,
				'B12':PredefinedWavelengths.S2A_MSI_13}
								    
			return Wavelength(bandSelect[bandname])

		def toa_to_rad(bandname):
			"""
			Converts top of atmosphere reflectance to at-sensor radiance
			"""
			
			# solar exoatmospheric spectral irradiance
			ESUN = info['SOLAR_IRRADIANCE_'+bandname]
			solar_angle_correction = math.cos(math.radians(solar_z))
			
			# Earth-Sun distance (from day of year)
			doy = scene_date.timetuple().tm_yday
			d = 1 - 0.01672 * math.cos(0.9856 * (doy-4))# http://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year
			
			# conversion factor
			multiplier = ESUN*solar_angle_correction/(math.pi*d**2)
			
			# at-sensor radiance
			rad = img.select(bandname).multiply(multiplier)
			return rad
			
		
		def surface_reflectance(bandname):
			"""
			Calculate surface reflectance from at-sensor radiance given waveband name
			"""
			
			# run 6S for this waveband
			s.wavelength = spectralResponseFunction(bandname)
			s.run()
			
			# extract 6S outputs
			Edir = s.outputs.direct_solar_irradiance             #direct solar irradiance
			Edif = s.outputs.diffuse_solar_irradiance            #diffuse solar irradiance
			Lp   = s.outputs.atmospheric_intrinsic_radiance      #path radiance
			absorb  = s.outputs.trans['global_gas'].upward       #absorption transmissivity
			scatter = s.outputs.trans['total_scattering'].upward #scattering transmissivity
			tau2 = absorb*scatter                                #total transmissivity
			
			# radiance to surface reflectance
			rad = toa_to_rad(bandname)
			
			ref = rad.subtract(Lp).multiply(math.pi).divide(tau2*(Edir+Edif))
			
			return ref
		
		# all wavebands
		output = img.select('QA60')
		for band in ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']:
			output = output.addBands(surface_reflectance(band))
			
		self.env.feature += 1
		
		return output #.addBands(TDOMMask)


	def getSentinel2(self):
	
		s2s = ee.ImageCollection('COPERNICUS/S2').filterDate(self.env.startDate,self.env.endDate) \
	                                                 .filterBounds(self.env.location) \
							 .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT',self.env.metadataCloudCoverMax)) \
						
#		2s = ee.ImageCollection([ee.Image("COPERNICUS/S2/20170606T153621_20170606T154217_T17MRT"),ee.Image("COPERNICUS/S2/20170606T153621_20170606T154217_T17MRT")])		
		s2sAll = ee.ImageCollection('COPERNICUS/S2').filterBounds(self.env.location) \
                                                 .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT',self.env.metadataCloudCoverMax)) \
#						 .map(self.QAMaskCloud)

		print("got" + str(s2s.size().getInfo()) + " images")

		if self.env.shadowMask == True:
			print("applying shadow mask..")
			s2s = self.maskShadows(s2s,s2sAll)

		print("scaling bands..")
		s2s = s2s.map(self.scaleS2) #.select(self.env.s2BandsIn,self.env.s2BandsOut)
 
		self.collectionMeta = s2s.getInfo()['features']
		
		if self.env.calcSR == True:
			print("calculate surface reflectance using 6s..")
			s2s = s2s.map(self.TOAtoSR).select(self.env.s2BandsIn,self.env.s2BandsOut)	
		print(ee.Image(s2s.first()).bandNames().getInfo())


		if self.env.QAcloudMask == True:
			print("use QA band for cloud Masking")
			s2s = s2s.map(self.QAMaskCloud)

		print(ee.Image(s2s.first()).bandNames().getInfo())
			
		if self.env.cloudMask == True:
			print("sentinel cloud score...")
			s2s = s2s.map(self.sentinelCloudScore)
			s2s = self.cloudMasking(s2s)

		print(ee.Image(s2s.first()).bandNames().getInfo())
		img = ee.Image(s2s.first())
#		print(img.geometry().bounds(1).getInfo())


		if self.env.brdf == True:
			print("apply brdf correction..")
			s2s = s2s.map(self.brdf)
			##print(s2s.first().getInfo())
			#img = self.brdf(s2s.first())
				
	
	
		if self.env.terrainCorrection == True:
			print("apply terrain correction..")
			s2s = s2s.map(self.terrain)
		
		print(ee.Image(s2s.first()).bandNames().getInfo())	
		print("calculating medoid")
		img = self.medoidMosaic(s2s)
		
		print(img.bandNames().getInfo())
		print("rescaling..")
		img = self.reScaleS2(img)
		
		return img
		

	def scaleS2(self,img):
		t = ee.Image(img.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']));
		t = t.divide(10000)
		t = t.addBands(img.select(['QA60']));
		out = t.copyProperties(img).copyProperties(img,['system:time_start','system:footprint']).set("centroid",img.geometry().centroid())

		return out;

	def reScaleS2(self,img):
		t = ee.Image(img.select(['red','green','blue']));
		t = t.multiply(10000)
		#t = t.addBands(img.select(['QA60']));
		out = ee.Image(t.copyProperties(img).copyProperties(img,['system:time_start'])).int16()

		return out;

		

	# Function to mask clouds using the Sentinel-2 QA band.
	def QAMaskCloud(self,image):
		qa = image.select('QA60').int16();
		
		# Bits 10 and 11 are clouds and cirrus, respectively.
		cloudBitMask = int(math.pow(2, 10));
		cirrusBitMask = int(math.pow(2, 11));
		
		# Both flags should be set to zero, indicating clear conditions.
		mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0));
		
		# Return the masked and scaled data.
		return image.updateMask(mask);

	def sentinelCloudScore(self,img):
		"""
		Computes spectral indices of cloudyness and take the minimum of them.
    
		Each spectral index is fairly lenient because the group minimum 
		is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
		
		originally written by Matt Hancher for Landsat imagery
		adapted to Sentinel by Chris Hewig and Ian Housman
		"""

        
		def rescale(img, thresholds):
			"""
			Linear stretch of image between two threshold values.
			"""
			return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])
        
		# cloud until proven otherwise
		score = ee.Image(1)
		blueCirrusScore = ee.Image(0)
		
		# clouds are reasonably bright
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['blue']), [0.1, 0.5]))
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cb']), [0.1, 0.5]))
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cirrus']), [0.1, 0.3]))
		score = score.min(blueCirrusScore)
		
		score = score.min(rescale(img.select(['red']).add(img.select(['green'])).add(img.select('blue')), [0.2, 0.8]))
		score = score.min(rescale(img.select(['nir1']).add(img.select(['swir1'])).add(img.select('swir2')), [0.3, 0.8]))

		# clouds are moist
		ndsi = img.normalizedDifference(['green','swir1'])
		score=score.min(rescale(ndsi, [0.8, 0.6]))
		score = score.multiply(100).byte();
		score = score.clamp(0,100);
		
		return img.addBands(score.rename(['cloudScore']))         
	
	
	def cloudMasking(self,collection):

		
		def maskClouds(img):
			cloudMask = img.select(['cloudScore']).lt(self.env.cloudScoreThresh)\
											      .focal_min(self.env.dilatePixels) \
											      .focal_max(self.env.contractPixels) \
											      .rename(['cloudMask'])    
            
			return img.updateMask(cloudMask).addBands(cloudMask);

		
		# Find low cloud score pctl for each pixel to avoid comission errors
		#minCloudScore = collection.select(['cloudScore']).reduce(ee.Reducer.percentile([self.env.cloudScorePctl]));

		collection = collection.map(maskClouds)
		
		return collection
		

	def maskShadows(self,collection,allCollection):

		def TDOM(image):
			zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
			irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
			TDOMMask = zScore.lt(self.env.zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
							 .And(irSum.lt(self.env.shadowSumThresh)).Not()
			TDOMMask = TDOMMask.focal_min(self.env.dilatePixels)
			return image.addBands(TDOMMask.rename(['TDOMMask']))

		def mask(image):
			outimg = image.updateMask(image.select(['TDOMMask']))
			return outimg

		shadowSumBands = ['B8','B11']

		# Get some pixel-wise stats for the time series
		irStdDev = allCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = allCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

		# Mask out dark dark outliers
		collection_tdom = collection.map(TDOM)

		return collection_tdom.map(mask)

	def addAllTasselCapIndices(self,img): 
		""" Function to get all tasselCap indices """
		
		def getTasseledCap(img):
			"""Function to compute the Tasseled Cap transformation and return an image"""
			
			coefficients = ee.Array([
				[0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
				[-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
				[0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
				[-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
				[-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
				[0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
			]);
		
			bands=ee.List(['blue','green','red','nir1','swir1','swir2'])
			
			# Make an Array Image, with a 1-D Array per pixel.
			arrayImage1D = img.select(bands).toArray()
		
			# Make an Array Image with a 2-D Array per pixel, 6x1.
			arrayImage2D = arrayImage1D.toArray(1)
		
			componentsImage = ee.Image(coefficients).matrixMultiply(arrayImage2D).arrayProject([0]).arrayFlatten([['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]).float();
	  
			# Get a multi-band image with TC-named bands.
			return img.addBands(componentsImage);	
			
			
		def addTCAngles(img):

			""" Function to add Tasseled Cap angles and distances to an image. Assumes image has bands: 'brightness', 'greenness', and 'wetness'."""
				
			# Select brightness, greenness, and wetness bands	
			brightness = img.select('brightness');
			greenness = img.select('greenness');
			wetness = img.select('wetness');
	  
			# Calculate Tasseled Cap angles and distances
			tcAngleBG = brightness.atan2(greenness).divide(math.pi).rename(['tcAngleBG']);
			tcAngleGW = greenness.atan2(wetness).divide(math.pi).rename(['tcAngleGW']);
			tcAngleBW = brightness.atan2(wetness).divide(math.pi).rename(['tcAngleBW']);
			tcDistBG = brightness.hypot(greenness).rename(['tcDistBG']);
			tcDistGW = greenness.hypot(wetness).rename(['tcDistGW']);
			tcDistBW = brightness.hypot(wetness).rename(['tcDistBW']);
			img = img.addBands(tcAngleBG).addBands(tcAngleGW).addBands(tcAngleBW).addBands(tcDistBG).addBands(tcDistGW).addBands(tcDistBW);
			
			return img;
	
	
		img = getTasseledCap(img)
		img = addTCAngles(img)
		return img

 
	def tcbwi(self,img):   
		
		out = img.expression('(-Abg*Dbg)+(Agw*Dgw**2)/(Awb*Dwb)',{
						'Agw':img.select('tcAngleGW'),
						'Awb':img.select('tcAngleBW'),
						'Abg':img.select('tcAngleBG'),
						'Dgw':img.select('tcDistGW'),
						'Dwb':img.select('tcDistBW'),
						'Dbg':img.select('tcDistBG'),
						}).rename(['tcbwi']);
		
		return img.addBands(out).addBands(out);

	def terrain(self,img):
		degree2radian = 0.01745;
 
		def topoCorr_IC(img):
			
			dem = ee.Image("USGS/SRTMGL1_003").unmask(0)
						
			# Extract image metadata about solar position
			SZ_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'))).multiply(3.14159265359).divide(180).clip(img.geometry().buffer(10000)); 
			SA_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).multiply(3.14159265359).divide(180)).clip(img.geometry().buffer(10000)); 
			
			# Creat terrain layers
			slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
			slp_rad = ee.Terrain.slope(dem).multiply(degree2radian).clip(img.geometry().buffer(10000));
			asp_rad = ee.Terrain.aspect(dem).multiply(degree2radian).clip(img.geometry().buffer(10000));
  			
			# Calculate the Illumination Condition (IC)
			# slope part of the illumination condition
			cosZ = SZ_rad.cos();
			cosS = slp_rad.cos();
			slope_illumination = cosS.expression("cosZ * cosS", \
												{'cosZ': cosZ, 'cosS': cosS.select('slope')});
	
			# aspect part of the illumination condition
			sinZ = SZ_rad.sin(); 
			sinS = slp_rad.sin();
			cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
			aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", \
                                           {'sinZ': sinZ, \
                                            'sinS': sinS, \
                                            'cosAziDiff': cosAziDiff});
	
			# full illumination condition (IC)
			ic = slope_illumination.add(aspect_illumination);
			
			# Add IC to original image
			img_plus_ic = ee.Image(img.addBands(ic.rename(['IC'])).addBands(cosZ.rename(['cosZ'])).addBands(cosS.rename(['cosS'])).addBands(slp.rename(['slope'])));
			
			return ee.Image(img_plus_ic);
 
			# Calculate the Illumination Condition (IC)
			# slope part of the illumination condition
			cosZ = SZ_rad.cos();
			cosS = slp_rad.cos();
			slope_illumination = cosS.expression("cosZ * cosS", \
												{'cosZ': cosZ, 'cosS': cosS.select('slope')});
			
			
			# aspect part of the illumination condition
			sinZ = SZ_rad.sin(); 
			sinS = slp_rad.sin();
			cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
			aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", \
                                           {'sinZ': sinZ, \
                                            'sinS': sinS, \
                                            'cosAziDiff': cosAziDiff});
			
			# full illumination condition (IC)
			ic = slope_illumination.add(aspect_illumination);
			
			

			# Add IC to original image
			img_plus_ic = ee.Image(img.addBands(ic.rename(['IC'])).addBands(cosZ.rename(['cosZ'])).addBands(cosS.rename(['cosS'])).addBands(slp.rename(['slope'])));
			
			return ee.Image(img_plus_ic);
 
		def topoCorr_SCSc(img):
			img_plus_ic = img;
			mask1 = img_plus_ic.select('nir1').gt(-0.1);
			mask2 = img_plus_ic.select('slope').gte(5) \
                            .And(img_plus_ic.select('IC').gte(0)) \
                            .And(img_plus_ic.select('nir1').gt(-0.1));

			img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));

			bandList = self.env.exportBands   			
 
			def applyBands(image):
				#cb = apply_SCSccorr('cb').select(['cb'])
				blue = apply_SCSccorr('blue').select(['blue'])
				green = apply_SCSccorr('green').select(['green'])
				red = apply_SCSccorr('red').select(['red'])
				re1 = apply_SCSccorr('re1').select(['re1'])
				re2 = apply_SCSccorr('re2').select(['re2'])
				re3 = apply_SCSccorr('re3').select(['re3'])
				nir1 = apply_SCSccorr('nir1').select(['nir1'])
				nir2 = apply_SCSccorr('nir2').select(['nir2'])
				swir1 = apply_SCSccorr('swir1').select(['swir1'])
				swir2 = apply_SCSccorr('swir2').select(['swir2'])

				return image.select([]).addBands([blue,green,red,re1,re2,re3,nir1,nir2,swir1,swir2])


			def apply_SCSccorr(band):
				method = 'SCSc';
		
				out = img_plus_ic_mask2.select('IC', band).reduceRegion(reducer= ee.Reducer.linearFit(), \
											geometry= ee.Geometry(img.geometry().buffer(-5000)), \
											scale= 100, \
											maxPixels = 1e13); 

				#out_a = ee.Number(out.get('scale'));
				#out_b = ee.Number(out.get('offset'));
				
				#out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')));
				out_c = 1	
				
				# apply the SCSc correction
				SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
															'image': img_plus_ic_mask2.select([band]),
															'ic': img_plus_ic_mask2.select('IC'),
															'cosB': img_plus_ic_mask2.select('cosS'),
															'cosZ': img_plus_ic_mask2.select('cosZ'),
															'cvalue': out_c });
      
				return ee.Image(SCSc_output);
																  
			#img_SCSccorr = ee.Image([apply_SCSccorr(band) for band in bandList]).addBands(img_plus_ic.select('IC'));
			img_SCSccorr = applyBands(img).select(bandList).addBands(img_plus_ic.select('IC'))
		
			bandList_IC = ee.List([bandList, 'IC']).flatten();
			
			img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
  			
			return img_SCSccorr.unmask(img_plus_ic.select(bandList)) 
	
		
		
		img = topoCorr_IC(img)
		img = topoCorr_SCSc(img)
		
		return img
  	

	def medoidMosaic(self,collection):
		""" medoid composite with equal weight among indices """
                    
		collection = collection.select(self.env.exportBands)

		bandNames = self.env.exportBands;
		bandNumbers = ee.List.sequence(1,bandNames.length());
        
		# calculate medion
		median = ee.ImageCollection(collection).median()
        
		def subtractmedian(img):
			diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
			return diff.reduce('sum').addBands(img);
        
		medoid = collection.map(subtractmedian)
  
		medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(bandNames.length().add(1))).select(bandNumbers,bandNames);
  
		return medoid;

 
	def brdf(self,img):   
		

	
		def _apply(image, kvol, kvol0):
			blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
			green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
			red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
			re1 = _correct_band(image, 're1', kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845)
			re2 = _correct_band(image, 're2', kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003)
			re3 = _correct_band(image, 're3', kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197)
			nir1 = _correct_band(image, 'nir1', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
			nir2 = _correct_band(image, 'nir2', kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611)
			swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
			swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
			return replace_bands(image, [blue, green, red,re1,re2,re3, nir1,nir2, swir1, swir2])


		def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
			"""fiso + fvol * kvol + fgeo * kgeo"""
			iso = ee.Image(f_iso)
			geo = ee.Image(f_geo)
			vol = ee.Image(f_vol)
			pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
			pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
			cfac = pred0.divide(pred).rename(['cfac'])
			corr = image.select(band_name).multiply(cfac).rename([band_name])
			return corr


		def _kvol(sunAz, sunZen, viewAz, viewZen):
			"""Calculate kvol kernel.
			From Lucht et al. 2000
			Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""
			
			relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
			pa1 = viewZen.cos() \
				.multiply(sunZen.cos())
			pa2 = viewZen.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle1 = pa1.add(pa2)
			phase_angle = phase_angle1.acos()
			p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
			p2 = p1.multiply(phase_angle1)
			p3 = p2.add(phase_angle.sin())
			p4 = sunZen.cos().add(viewZen.cos())
			p5 = ee.Image(PI().divide(4))

			kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

			viewZen0 = ee.Image(0)
			pa10 = viewZen0.cos() \
				.multiply(sunZen.cos())
			pa20 = viewZen0.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle10 = pa10.add(pa20)
			phase_angle0 = phase_angle10.acos()
			p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
			p20 = p10.multiply(phase_angle10)
			p30 = p20.add(phase_angle0.sin())
			p40 = sunZen.cos().add(viewZen0.cos())
			p50 = ee.Image(PI().divide(4))

			kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

			return (kvol, kvol0)
         
		date = img.date()
		footprint = determine_footprint(img)
		(sunAz, sunZen) = sun_angles.create(date, footprint)
		(viewAz, viewZen) = view_angles.create(footprint)
		(kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)
		img = ee.Image(_apply(img, kvol.multiply(PI()), kvol0.multiply(PI())))
		return img



     		
if __name__ == "__main__":        
	
	img = functions().getSentinel2().select(["red","green","blue"]) #,"cloudScore","cloudMask","TDOMMask"]).multiply(10000)
	#img = functions().getSentinel2().select(["TDOMMask"])	
	print(img.bandNames().getInfo())	
	
	geom = ee.Image(img).geometry().getInfo()
	countries = ee.FeatureCollection('ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw');
	geom  = countries.filter(ee.Filter.inList('Country', ['Ecuador'])).geometry().bounds().getInfo();

	import time
	t= time.strftime("%Y%m%d_%H%M%S")
	task_ordered= ee.batch.Export.image.toAsset(image=img, 
								  description="tempwater", 
								  assetId="users/servirmekong/temp/077ters209"+t ,
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  scale=20)
	
	
	task_ordered.start() 
