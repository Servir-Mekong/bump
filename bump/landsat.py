# Landsat package


import ee
import math 

class env(object):

    def __init__(self):
        """Initialize the environment."""   
         
        # Initialize the Earth Engine object, using the authentication credentials.
                
        ee.Initialize()
        
        self.startDate = "2016-01-01"
        self.endDate = "2016-05-01"
        self.location = ee.Geometry.Point([105.216064453125,19.041348796589016])
        
        self.metadataCloudCoverMax = 30
        self.cloudThreshold = 10
        self.hazeThresh = 200
              
        self.maskSR = True
        self.cloudMask = True
        self.hazeMask = True
        self.shadowMask = True
        
        self.divideBands = ee.List(['blue','green','red','nir','swir1','swir2'])
        self.bandNamesLandsat = ee.List(['blue','green','red','nir','swir1','thermal','swir2','sr_atmos_opacity','pixel_qa','radsat_qa'])
        self.sensorBandDictLandsatSR = ee.Dictionary({'L8' : ee.List([1,2,3,4,5,7,6,9,10,11]),\
                                                      'L7' : ee.List([0,1,2,3,4,5,6,7,9,10]),\
                                                      'L5' : ee.List([0,1,2,3,4,5,6,7,9,10]),\
                                                      'L4' : ee.List([0,1,2,3,4,5,6,7,9,10])})
                                                      

class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 	
		
	def getLandsat(self):
		landsat8 =  ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(self.env.startDate,self.env.endDate).filterBounds(self.env.location)
		landsat8 = landsat8.filterMetadata('CLOUD_COVER','less_than',self.env.metadataCloudCoverMax)
		landsat8 = landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)
		
		print ee.Image(landsat8.first()).bandNames().getInfo()                 
		
		if landsat8.size().getInfo() > 0:
			
			# mask clouds using the QA band
			#if self.env.maskSR == True:
			#	print "removing clouds" 
			#	landsat8 = landsat8.map(self.CloudMaskSRL8)    
					
			# mask clouds using cloud mask function
			#if self.env.cloudMask == True:
			#	print "removing some more clouds"
			#	landsat8 = landsat8.map(self.maskClouds)

				# mask clouds using cloud mask function
			#if self.env.hazeMask == True:
			#	print "removing haze"
			#	landsat8 = landsat8.map(self.maskHaze)

			# mask clouds using cloud mask function
			#if self.env.shadowMask == True:
			#	print "shadow masking"
			#	self.fullCollection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterBounds(self.env.location).select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)  
			#	landsat8 = self.maskShadows(landsat8)		
				
			landsat8 = landsat8.map(self.scaleLandsat)
			#landsat8 = landsat8.map(self.getTasseledCap)
			#landsat8 = landsat8.map(self.addTCAngles)
			#slandsat8 = landsat8.map(self.reScaleLandsat)   
	

		return landsat8
       
	
	def CloudMaskSRL8(self,img):
		"""apply cf-mask Landsat""" 
		QA = img.select("pixel_qa")
		
		shadow = QA.bitwiseAnd(8).neq(0);
		cloud =  QA.bitwiseAnd(32).neq(0);
		return img.updateMask(shadow.Not()).updateMask(cloud.Not()).copyProperties(img)		
         
	def scaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		#thermal = ee.Image(img).select(ee.List(['thermal'])).multiply(0.1)
		scaled = ee.Image(img).select(self.env.divideBands).multiply(0.0001)
		#image = ee.Image(scaled).addBands(thermal)	
		
		#return ee.Image(image.copyProperties(img))
		return ee.Image(scaled)

	def reScaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
        
		thermalBand = ee.List(['thermal'])
		thermal = ee.Image(img).select(thermalBand).multiply(10)
                
		otherBands = ee.Image(img).bandNames().removeAll(thermalBand)
		scaled = ee.Image(img).select(otherBands).divide(0.0001)
        
		image = ee.Image(scaled.addBands(thermal)).int16()
        
		return image.copyProperties(img)

	def maskHaze(self,img):
		""" mask haze """
		opa = ee.Image(img.select(['sr_atmos_opacity']).multiply(0.001))
		haze = opa.gt(self.env.hazeThresh)
		return img.updateMask(haze.Not())
 

	def maskClouds(self,img):
		"""
		Computes spectral indices of cloudyness and take the minimum of them.
		
		Each spectral index is fairly lenient because the group minimum 
		is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
		originally written by Matt Hancher for Landsat imageryadapted to Sentinel by Chris Hewig and Ian Housman
		"""
		
		score = ee.Image(1.0);
		# Clouds are reasonably bright in the blue band.
		blue_rescale = img.select('blue').subtract(ee.Number(0.1)).divide(ee.Number(0.3).subtract(ee.Number(0.1)))
		score = score.min(blue_rescale);

		# Clouds are reasonably bright in all visible bands.
		visible = img.select('red').add(img.select('green')).add(img.select('blue'))
		visible_rescale = visible.subtract(ee.Number(0.2)).divide(ee.Number(0.8).subtract(ee.Number(0.2)))
		score = score.min(visible_rescale);

		# Clouds are reasonably bright in all infrared bands.
		infrared = img.select('nir').add(img.select('swir1')).add(img.select('swir2'))
		infrared_rescale = infrared.subtract(ee.Number(0.3)).divide(ee.Number(0.8).subtract(ee.Number(0.3)))
		score = score.min(infrared_rescale);

		# Clouds are reasonably cool in temperature.
		temp_rescale = img.select('thermal').subtract(ee.Number(300)).divide(ee.Number(290).subtract(ee.Number(300)))
		score = score.min(temp_rescale);

		# However, clouds are not snow.
		ndsi = img.normalizedDifference(['green', 'swir1']);
		ndsi_rescale = ndsi.subtract(ee.Number(0.8)).divide(ee.Number(0.6).subtract(ee.Number(0.8)))
		score =  score.min(ndsi_rescale).multiply(100).byte();
		mask = score.lt(self.env.cloudThreshold).rename(['cloudMask']);
		img = img.updateMask(mask);
        
		return img;
        
	def maskShadows(self,collection,zScoreThresh=-0.8,shadowSumThresh=0.35,dilatePixels=2):

		def TDOM(image):
			zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
			irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
			TDOMMask = zScore.lt(zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
				.And(irSum.lt(shadowSumThresh)).Not()
			TDOMMask = TDOMMask.focal_min(dilatePixels)
			
			return image.updateMask(TDOMMask)
			
		shadowSumBands = ['nir','swir1']

		# Get some pixel-wise stats for the time series
		irStdDev = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

		# Mask out dark dark outliers
		collection_tdom = collection.map(TDOM)

		return collection_tdom



	def getTasseledCap(self,img):
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


	def addTCAngles(self,img):

		""" Function to add Tasseled Cap angles and distances to an image. Assumes image has bands: 'brightness', 'greenness', and 'wetness'. """
		
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


		img = img.addBands(img.expression('(-Abg*Dbg)+(Agw*Dgw**2)/(Awb*Dwb)',{
								'Agw':img.select(['tcAngleGW']),
								'Awb':img.select(['tcAngleBW']),
								'Abg':img.select(['tcAngleBG']),
								'Dgw':img.select(['tcDistGW']),
								'Dwb':img.select(['tcDistBW']),
								'Dbg':img.select(['tcDistBG'])}).rename(['tcbwi']));
			
		return img;
	

            
 
        
if __name__ == "__main__":        
	
	landsatImages = functions().getLandsat()
	
	print landsatImages.size().getInfo()
	
	img = ee.Image(landsatImages.first())
	geom = ee.Image(landsatImages.first()).geometry().getInfo()
	print img.bandNames().getInfo()
	
	task_ordered= ee.batch.Export.image.toAsset(image=img, 
								  description="tempwater", 
								  assetId="users/servirmekong/temp/tempwwater10" ,
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  scale=30)
	
	
	task_ordered.start() 
