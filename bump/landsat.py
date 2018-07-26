# Landsat package


import ee
import math 
from utils import *

class env(object):

    def __init__(self):
        """Initialize the environment."""   
         
        # Initialize the Earth Engine object, using the authentication credentials.
                
        ee.Initialize()
        
        self.startDate = "2016-01-01"
        self.endDate = "2016-09-01"
        self.location = ee.Geometry.Point([105.216064453125,19.041348796589016])
        
        self.metadataCloudCoverMax = 30
        self.cloudThreshold = 10
        self.hazeThresh = 200
              
        self.maskSR = True
        self.cloudMask = True
        self.hazeMask = True
        self.shadowMask = True
        self.brdfCorrect = True
        self.terrainCorrection = True
        
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
			if self.env.maskSR == True:
				print "removing clouds" 
				landsat8 = landsat8.map(self.CloudMaskSRL8)    
					
			# mask clouds using cloud mask function
			if self.env.hazeMask == True:
				print "removing haze"
				landsat8 = landsat8.map(self.maskHaze)

			# mask clouds using cloud mask function
			if self.env.shadowMask == True:
				print "shadow masking"
				self.fullCollection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterBounds(self.env.location).select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)  
				landsat8 = self.maskShadows(landsat8)		
			
			landsat8 = landsat8.map(self.scaleLandsat)
			
			# mask clouds using cloud mask function
			if self.env.cloudMask == True:
				print "removing some more clouds"
				landsat8 = landsat8.map(self.maskClouds)
					
			if self.env.brdfCorrect == True:
				landsat8 = landsat8.map(self.brdf)
			
			img = ee.Image(landsat8.first())
				
			if self.env.terrainCorrection == True:
				print "terrain correction"
				landsat8 = landsat8.map(self.terrain)
			
			#img = self.terrain(img)
			#print img.getInfo()
			#
			#landsat8 = landsat8.map(self.addAllTasselCapIndices)


		return landsat8
       
	
	def CloudMaskSRL8(self,img):
		"""apply cf-mask Landsat""" 
		QA = img.select("pixel_qa")
		
		shadow = QA.bitwiseAnd(8).neq(0);
		cloud =  QA.bitwiseAnd(32).neq(0);
		return img.updateMask(shadow.Not()).updateMask(cloud.Not()).copyProperties(img)		
         
	def scaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		thermal = img.select(ee.List(['thermal'])).multiply(0.1)
		scaled = ee.Image(img).select(self.env.divideBands).multiply(ee.Number(0.0001))
		
		return img.select([]).addBands(scaled).addBands(thermal)
		
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

 
 	def terrain(self,img):   
		degree2radian = 0.01745;
 
		# Extract solar zenith and calculate incidence angle (i)
		# Load USGS/SRTMGL1_003 DEM
		terrain = ee.call('Terrain', ee.Image('USGS/SRTMGL1_003'));
		# Extract slope in radians for each pixel in the image  
		p = terrain.select(['slope']).multiply(degree2radian);
		# Extract solar zenith angle from the image
		z = ee.Number(img.get('SOLAR_ZENITH_ANGLE')).multiply(degree2radian);
		# Extract solar azimuth from the image
		az = ee.Image(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE')).multiply(degree2radian));
		# Extract aspect in radians for each pixel in the image
		o = terrain.select(['aspect']).multiply(degree2radian);
		cosao = (az.subtract(o)).cos()
		# Calculate the cosine of the local solar incidence for every pixel in the image in radians 
		cosi = img.expression('((cosp*cosz) + ((sinp*sinz)*(cosao)))', 
		{
		  'cosp': p.cos(),
		  'cosz': z.cos(),
		  'sinp': p.sin(),
		  'sinz': z.sin(),
		  'az' : az,
		  'o' : o,
		  'cosao': cosao
		});

		# Create the image to apply the linear regression.The first band
		# is the cosi and the second band is the response variable, the reflectance (the bands).
		# L (y) = a + b*cosi(x); a = intercept, b = slope

		# Dependent: Reflectance
		y = img.select(['red', 'nir', 'blue']);
		# Independent: (cosi)
		x = cosi;
		# Intercept: a
		a = ee.Image(1).rename(['a']);
		# create an image collection with the three variables by concatenating them
		reg_img = ee.Image.cat(a,x,y);
		# specify the linear regression reducer
		lr_reducer = ee.Reducer.linearRegression(2,3);

		# fit the model
		fit = reg_img.reduceRegion(
				reducer= lr_reducer,
				geometry= ee.Image(img).geometry(),
				scale= 30,
				maxPixels= 1e9);
				
		fit = fit.combine({"coefficients": ee.Array([[1],[1]])}, False);

		# Get the coefficients as a nested list, 
		# cast it to an array, and get just the selected column
		slo = (ee.Array(fit.get('coefficients')).get([1,0]));
		intercept = (ee.Array(fit.get('coefficients')).get([0,0]));

		# Calculate C parameter C= a/b
		C = intercept.divide(slo);

		image = img.expression(
					'((img * ((cosp*cosz) + C))/(cosi + C))',
					{
					  'img': img,
					  'cosp': p.cos(),
					  'cosz': z.cos(),
					  'cosi': cosi,
					  'C': C
				  });

		return img.select([]).addBands(image).copyProperties(img)
 
 
	def brdf(self,img):   
		
		import sun_angles
		import view_angles

	
		def _apply(image, kvol, kvol0):
			blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
			green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
			red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
			nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
			swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
			swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
			return replace_bands(image, [blue, green, red, nir, swir1, swir2])


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
		return _apply(img, kvol.multiply(PI()), kvol0.multiply(PI()))

        
if __name__ == "__main__":        
	
	landsatImages = functions().getLandsat()
	
	print landsatImages.size().getInfo()
	
	img = ee.Image(landsatImages.first())
	geom = ee.Image(landsatImages.first()).geometry().getInfo()
	
	print img.bandNames().getInfo()
	
	task_ordered= ee.batch.Export.image.toAsset(image=img, 
								  description="tempwater", 
								  assetId="users/servirmekong/temp/tempwwater22" ,
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  scale=100)
	
	
	task_ordered.start() 
