
import ee

class indices():

	def __init__(self):
		
		# list with functions to call for each index
		self.functionList = {"ND_blue_green" : self.ND_blue_green, \
							 "ND_blue_red" : self.ND_blue_red, \
							 "ND_blue_nir" : self.ND_blue_nir, \
							 "ND_blue_swir1" : self.ND_blue_swir1, \
							 "ND_blue_swir2" : self.ND_blue_swir2, \
							 "ND_green_red" : self.ND_green_red, \
							 "ND_green_nir" : self.ND_green_nir, \
							 "ND_green_swir1" : self.ND_green_swir1, \
							 "ND_green_swir2" : self.ND_green_swir2, \
							 "ND_red_swir1" : self.ND_red_swir1, \
							 "ND_red_swir2" : self.ND_red_swir2, \
							 "ND_nir_red" : self.ND_nir_red, \
							 "ND_nir_swir1" : self.ND_nir_swir1, \
							 "ND_nir_swir2" : self.ND_nir_swir2, \
							 "ND_swir1_swir2" : self.ND_swir1_swir2, \
							 "R_swir1_nir" : self.R_swir1_nir, \
							 "R_red_swir1" : self.R_red_swir1, \
							 "EVI" : self.EVI, \
							 "SAVI" : self.SAVI, \
							 "IBI" : self.IBI}


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

	def ND_blue_green(self,img):
		img = img.addBands(img.normalizedDifference(['blue','green']).rename(['ND_blue_green']));
		return img
	
	def ND_blue_red(self,img):
		img = img.addBands(img.normalizedDifference(['blue','red']).rename(['ND_blue_red']));
		return img
	
	def ND_blue_nir(self,img):
		img = img.addBands(img.normalizedDifference(['blue','nir']).rename(['ND_blue_nir']));
		return img
	
	def ND_blue_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['blue','swir1']).rename(['ND_blue_swir1']));
		return img
	
	def ND_blue_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['blue','swir2']).rename(['ND_blue_swir2']));
		return img

	def ND_green_red(self,img):
		img = img.addBands(img.normalizedDifference(['green','red']).rename(['ND_green_red']));
		return img
	
	def ND_green_nir(self,img):
		img = img.addBands(img.normalizedDifference(['green','nir']).rename(['ND_green_nir']));  # NDWBI
		return img
	
	def ND_green_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['green','swir1']).rename(['ND_green_swir1']));  # NDSI, MNDWI
		return img
	
	def ND_green_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['green','swir2']).rename(['ND_green_swir2']));
		return img
		
	def ND_red_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['red','swir1']).rename(['ND_red_swir1']));
		return img
			
	def ND_red_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['red','swir2']).rename(['ND_red_swir2']));
		return img

	def ND_nir_red(self,img):
		img = img.addBands(img.normalizedDifference(['nir','red']).rename(['ND_nir_red']));  # NDVI
		return img
	
	def ND_nir_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['nir','swir1']).rename(['ND_nir_swir1']));  # NDWI, LSWI, -NDBI
		return img
	
	def ND_nir_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['nir','swir2']).rename(['ND_nir_swir2']));  # NBR, MNDVI
		return img

	def ND_swir1_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['swir1','swir2']).rename(['ND_swir1_swir2']));
		return img
  
	def R_swir1_nir(self,img):
		# Add ratios
		img = img.addBands(img.select('swir1').divide(img.select('nir')).rename(['R_swir1_nir']));  # ratio 5/4
		return img
			
	def R_red_swir1(self,img):
		img = img.addBands(img.select('red').divide(img.select('swir1')).rename(['R_red_swir1']));  # ratio 3/5
		return img

	def EVI(self,img):
		#Add Enhanced Vegetation Index (EVI)
		evi = img.expression(
			'2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
			  'NIR': img.select('nir'),
			  'RED': img.select('red'),
			  'BLUE': img.select('blue')
		  }).float();
	
		img = img.addBands(evi.rename(['EVI']));

		return img
	  
	def SAVI(self,img):
		# Add Soil Adjust Vegetation Index (SAVI)
		# using L = 0.5;
		savi = img.expression(
			'(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
			  'NIR': img.select('nir'),
			  'RED': img.select('red')
		  }).float();
		img = img.addBands(savi.rename(['SAVI']));

		return img
	  
	def IBI(self,img):
		# Add Index-Based Built-Up Index (IBI)
		ibi_a = img.expression(
			'2*SWIR1/(SWIR1 + NIR)', {
			  'SWIR1': img.select('swir1'),
			  'NIR': img.select('nir')
			}).rename(['IBI_A']);
	

		ibi_b = img.expression(
			'(NIR/(NIR + RED)) + (GREEN/(GREEN + SWIR1))', {
			  'NIR': img.select('nir'),
			  'RED': img.select('red'),
			  'GREEN': img.select('green'),
			  'SWIR1': img.select('swir1')
			}).rename(['IBI_B']);
		
		ibi_a = ibi_a.addBands(ibi_b);
		ibi = ibi_a.normalizedDifference(['IBI_A','IBI_B']);
		img = img.addBands(ibi.rename(['IBI']));
		
		return img
