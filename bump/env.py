# main environment package for bump


class environment(object):
	
	def __init__(self):
		
		"""Initialize the environment."""   
		# Initialize the Earth Engine object, using the authentication credentials.            
		ee.Initialize()
		
		# set dates
        self.startYear = 2015
        self.endYear = 2015
        self.startJulian = 10
        self.endJulian = 50
        
        # define study area
        self.studyArea = [[103.876,18.552],[105.806,18.552],[105.806,19.999],[103.876,19.999],[103.876,18.552]]
	
