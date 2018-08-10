# Sentinel-2 package

import ee
import math 
import ee.mapclient

class env(object):

    def __init__(self):
        """Initialize the environment."""   
         
        # Initialize the Earth Engine object, using the authentication credentials.
                
        ee.Initialize()
        
        self.startDate = "2017-02-01"
        self.endDate = "2017-09-01"
        self.location = ee.Geometry.Point([105.216064453125,19.041348796589016])
        
        self.metadataCloudCoverMax = 30
        self.cloudThreshold = 10
        self.hazeThresh = 200
        
        self.s2BandsIn = ee.List(['QA60','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12'])
        self.s2BandsOut = ee.List(['QA60','cb','blue','green','red','re1','re2','re3','nir','nir2','waterVapor','cirrus','swir1','swir2'])
        self.divideBands = ee.List(['cb','blue','green','red','re1','re2','re3','nir','nir2','waterVapor','cirrus','swir1','swir2'])
      
        
         # 9. Cloud and cloud shadow masking parameters.
        # If cloudScoreTDOM is chosen
        # cloudScoreThresh: If using the cloudScoreTDOMShift method-Threshold for cloud 
        #    masking (lower number masks more clouds.  Between 10 and 30 generally works best)
        self.cloudScoreThresh = 0.3;

        # Percentile of cloud score to pull from time series to represent a minimum for 
        # the cloud score over time for a given pixel. Reduces commission errors over 
        # cool bright surfaces. Generally between 5 and 10 works well. 0 generally is a bit noisy
        self.cloudScorePctl = 0; 

        # zScoreThresh: Threshold for cloud shadow masking- lower number masks out 
        # less. Between -0.8 and -1.2 generally works well
        self.zScoreThresh = -0.8;

        # shadowSumThresh: Sum of IR bands to include as shadows within TDOM and the 
        # shadow shift method (lower number masks out less)
        self.shadowSumThresh = 0.30;

        # contractPixels: The radius of the number of pixels to contract (negative buffer) clouds and cloud shadows by. Intended to eliminate smaller cloud 
        #    patches that are likely errors (1.5 results in a -1 pixel buffer)(0.5 results in a -0 pixel buffer)
        # (1.5 or 2.5 generally is sufficient)
        self.contractPixels = 1.5; 

        # dilatePixels: The radius of the number of pixels to dilate (buffer) clouds 
        # and cloud shadows by. Intended to include edges of clouds/cloud shadows 
        # that are often missed (1.5 results in a 1 pixel buffer)(0.5 results in a 0 pixel buffer)
        # (2.5 or 3.5 generally is sufficient)
        self.dilatePixels = 3.5;

              
        self.maskSR = True
        self.cloudMask = True
        self.hazeMask = True
        self.shadowMask = True
        self.terrainCorrection = True


class surfaceReflectance():
	def __init__(self):
		self.env = environment()
    
	def TAOtoSR(self,img):
		
		cloudMask = img.select(['cloudScore'])
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
		km = alt/1000 # i.e. Py6S uses units of kilometers
		
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
		
		return output.addBands(cloudMask)


class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
		
	def getSentinel2(self):
		s2s = ee.ImageCollection('COPERNICUS/S2').filterDate(self.env.startDate,self.env.endDate) \
                                                 .filterBounds(self.env.location) \
												 .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT',self.env.metadataCloudCoverMax)) \
												 .select(self.env.s2BandsIn,self.env.s2BandsOut)
		
		s2s = s2s.map(self.scaleS2)
		s2s = s2s.map(self.QAMaskCloud)
		s2s = s2s.map(self.sentinelCloudScore)
			
		s2s = self.cloudMasking(s2s)
		s2s = self.maskShadows(s2s)
		
		s2s = s2s.map(self.addAllTasselCapIndices)
		s2s = s2s.map(self.tcbwi)
		s2s = s2s.map(self.terrain)
		
		img = ee.Image(s2s.first())
		
				
		return img
		

	def scaleS2(self,img):
		t = img.select(self.env.divideBands).divide(10000);
		t = t.addBands(img.select(['QA60']));
		return t.copyProperties(img).copyProperties(img,['system:time_start']).set("centroid",img.geometry().centroid())
		

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
		
		# clouds are reasonably bright
		score = score.min(rescale(img.select(['blue']), [0.1, 0.5]))
		score = score.min(rescale(img.select(['cb']), [0.1, 0.3]))
		score = score.min(rescale(img.select(['cb']).add(img.select(['cirrus'])), [0.15, 0.2]))
		score = score.min(rescale(img.select(['red']).add(img.select(['green'])).add(img.select('blue')), [0.2, 0.8]))
		
		# clouds are moist
		ndmi = img.normalizedDifference(['re3','swir1'])
		score=score.min(rescale(ndmi, [-0.1, 0.1]))
		
		# clouds are not snow.
		ndsi = img.normalizedDifference(['green', 'swir1'])
		score=score.min(rescale(ndsi, [0.8, 0.6])).rename(['cloudScore'])
		
		return img.addBands(score)         
	
	
	def cloudMasking(self,collection):

		
		def maskClouds(img):
			cloudMask = img.select(['cloudScore']).subtract(minCloudScore) \
											      .focal_min(self.env.dilatePixels) \
											      .focal_max(self.env.contractPixels) \
											      .lt(self.env.cloudScoreThresh).rename(['cloudMask'])    
            
			return img.updateMask(cloudMask).addBands(cloudMask);

		
		# Find low cloud score pctl for each pixel to avoid comission errors
		minCloudScore = collection.select(['cloudScore']).reduce(ee.Reducer.percentile([self.env.cloudScorePctl]));

		collection = collection.map(maskClouds)
		
		return collection
		

	def maskShadows(self,collection):

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

		shadowSumBands = ['nir','swir1']

		# Get some pixel-wise stats for the time series
		irStdDev = collection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = collection.select(shadowSumBands).reduce(ee.Reducer.mean())

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
		
			bands=ee.List(['blue','green','red','nir','swir1','swir2'])
			
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
			
			dem = ee.Image("USGS/SRTMGL1_003")
						
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
 
		def topoCorr_SCSc(img):
			img_plus_ic = img;
			mask1 = img_plus_ic.select('nir').gt(-0.1);
			mask2 = img_plus_ic.select('slope').gte(5) \
                            .And(img_plus_ic.select('IC').gte(0)) \
                            .And(img_plus_ic.select('nir').gt(-0.1));

			img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));

			bandList = ee.List(['cb','blue','green','red']) #,'re1','re2','re3','nir','nir2','waterVapor','cirrus','swir1','swir2']); # Specify Bands to topographically correct
    
			def apply_SCSccorr(band):
				method = 'SCSc';
			
				out = img_plus_ic_mask2.select('IC', band).reduceRegion(reducer= ee.Reducer.linearFit(), \
																			geometry= img.geometry().buffer(-5000), \
																			scale= 10, \
																			maxPixels = 1e13); 

				out_a = ee.Number(out.get('scale'));
				out_b = ee.Number(out.get('offset'));
				out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')));
				
				# apply the SCSc correction
				SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
															'image': img_plus_ic_mask2.select([band]),
															'ic': img_plus_ic_mask2.select('IC'),
															'cosB': img_plus_ic_mask2.select('cosS'),
															'cosZ': img_plus_ic_mask2.select('cosZ'),
															'cvalue': out_c });
      
				return ee.Image(SCSc_output);
			
				
			# need to fix this in to map.. 
			img_SCSccorr = img.select([]).addBands(apply_SCSccorr("cb")) \
										 .addBands(apply_SCSccorr("blue")) \
										 .addBands(apply_SCSccorr("green")) \
										 .addBands(apply_SCSccorr("red"))\
								#		 .addBands(apply_SCSccorr("re1")) \
								#		 .addBands(apply_SCSccorr("re2"))\
								#		 .addBands(apply_SCSccorr("re3"))\
								#		 .addBands(apply_SCSccorr("nir"))\
								#		 .addBands(apply_SCSccorr("nir2"))\
								#		 .addBands(apply_SCSccorr("waterVapor"))\
								#		 .addBands(apply_SCSccorr("cirrus"))\
								#		 .addBands(apply_SCSccorr("swir1"))\
								#		 .addBands(apply_SCSccorr("swir2"))\
	
	
			return img_SCSccorr.unmask(img_plus_ic.select(bandList)) 
	
		
		
		img = topoCorr_IC(img)
		img = topoCorr_SCSc(img)
		
		return img

     		
if __name__ == "__main__":        
	
	img = functions().getSentinel2()
	
	#print s2Images.bandNames().getInfo()

	#img = ee.Image(s2Images.first()) #.select(["tcbwi"])
	#print img.getInfo()
	geom = ee.Image(img).geometry().getInfo()

	#img = img.mask(img.gt(0))
	# create the vizualization parameters
	#viz = {'min':-1, 'max':1, 'bands':"tcbwi",'palette':"red,green,blue"};
 
	#ee.mapclient.centerMap(105.216064453125,19.0413,10)
	#ee.mapclient.addToMap(img,viz, "mymap")
	
	task_ordered= ee.batch.Export.image.toAsset(image=img, 
								  description="tempwater", 
								  assetId="users/servirmekong/temp/0waters209" ,
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  scale=150)
	
	
	task_ordered.start() 
