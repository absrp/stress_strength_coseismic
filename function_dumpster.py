import matplotlib.colors as colors
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import pyproj
import numpy as np

class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def dot(vA, vB):
	"""
	Takes the dot product of two vectors
	e.g. dot(line1,line2)
	"""
	return vA[0]*vB[0]+vA[1]*vB[1]

# from stack exchange
def measure_angle(lineA, lineB):
	"""
	Measure the angle between two line segments
	e.g. measure_angle(line1,line2)
	"""
	# Get nicer vector form
	vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
	vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
	# Get dot prod
	dot_prod = dot(vA, vB)
	# Get magnitudes
	magA = dot(vA, vA)**0.5
	magB = dot(vB, vB)**0.5
	# Get cosine value
	cos_ = dot_prod/magA/magB
	#    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
	dot_prod = dot(vA, vB)
	# Get magnitudes
	magA = dot(vA, vA)**0.5
	magB = dot(vB, vB)**0.5
	# Get cosine value
	cos_ = dot_prod/magA/magB
	# Get angle in radians and then convert to degrees
	angle = math.acos(dot_prod/magB/magA)
	# Basically doing angle <- angle mod 360
	ang_deg = math.degrees(angle)%360
	if ang_deg-180>=0:
		# As in if statement
		return 360 - ang_deg
	else: 			
		return ang_deg

def extract_coordinates_shapefile(shapefile_path, target_utm_zone, step=1):
	"""
	Extracts the lat,lon coordinates from a shapefile and converts them to UTM coordinates	
	"""
	gdf = gpd.read_file(shapefile_path)
	gdf = gdf.set_crs('EPSG:4326')
	utm_crs = f"+proj=utm +zone={target_utm_zone} +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
	gdf_utm = gdf.to_crs(utm_crs)
	x_coords = []
	y_coords = []

	counter = 0
	for line in gdf_utm['geometry']:
		for point in line.coords:
			x, y = point
			if counter % step == 0:
				# Only save every nth point as defined by the step parameter
				x_coords.append(x)
				y_coords.append(y)
			counter += 1

	return x_coords, y_coords

def extract_coordinates_stress(utm_crs, ystress, xstress):
	"""
	Extracts the lat,lon coordinates of the stress field and converts them to UTM coordinates matching the event coordinate system	
	"""
	# Perform the coordinate transformation
	source_crs = 'EPSG:4326'
	transformer = pyproj.Transformer.from_crs(source_crs, utm_crs, always_xy=True)
	lat_lon_points = [(lat, lon) for lat, lon in zip(ystress, xstress)]
	xstress_utm = []
	ystress_utm = []
	for lat, lon in lat_lon_points:
		utm_x, utm_y = transformer.transform(lon, lat)
		xstress_utm.append(utm_x)
		ystress_utm.append(utm_y)
	return xstress_utm, ystress_utm

def measure_angle_between_segments(x_coords, y_coords):
	"""
	Measure the angle between two line segments using measure_angle() function
	"""
	angles = []
	for i in range(1, len(x_coords) - 1):
		p0 = [x_coords[i - 1], y_coords[i - 1]]
		p1 = [x_coords[i], y_coords[i]]
		p2 = [x_coords[i + 1], y_coords[i + 1]]
		angles.append(measure_angle([p0, p1], [p1, p2]))
	return angles

def plot_segments_color_angle(vmin,vmax,x_coords,y_coords,angles,eventname):
    # Create a scatter plot of the extracted coordinates with colors based on angles in degrees
    plt.scatter(x_coords[1:-1], y_coords[1:-1], c=angles, cmap='magma', marker='o',  edgecolors='none')
    # Set the desired range for the colorbar
    plt.title(eventname)
    plt.axis('equal')
    plt.colorbar(label='Angle difference (degrees)')
    plt.show()

def crop_stress_to_event(x_coords, y_coords, xstress_utm, ystress_utm, SHmax):
	"""
	Crops the stress field to match the extent of the surface rupture of interest
	"""
	# Create a scatter plot of the extracted coordinates with colors based on angles in degrees
	x_range = [min(x_coords) - 10000, max(x_coords) + 10000]
	y_range = [min(y_coords) - 10000, max(y_coords) + 10000]

	# Initialize lists to store the cropped data
	cropped_xstress_utm = []
	# END: ed8c6549bwf9
	cropped_ystress_utm = []
	cropped_SHmax = []

	# Crop xstress_utm, ystress_utm, and SHmax to the specified range while keeping the same size
	for x, y, shmax in zip(xstress_utm, ystress_utm, SHmax):
		if x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]:
			cropped_xstress_utm.append(x)
			cropped_ystress_utm.append(y)
			cropped_SHmax.append(shmax)
		
	return cropped_xstress_utm, cropped_ystress_utm, cropped_SHmax

def unit_vector_stress(cropped_SHmax, mag=2000):
	"""
	Transforms the SHmax angles into unit vectors for visualization
	"""
	unit_vectors = []
    # Convert each angle to its corresponding unit vector
	for SHmax_degrees in cropped_SHmax:
		SHmax_radians = np.radians(SHmax_degrees-90)
		# Calculate the unit vector in 2D space
		x = np.cos(SHmax_radians)
		y = np.sin(SHmax_radians)
		z = 0  # Assuming 2D
		# Normalize the unit vector
		magnitude = np.sqrt(x**2 + y**2 + z**2)
		unit_vector = np.array([x / magnitude, y / magnitude, z / magnitude])
		unit_vectors.append(unit_vector)
		vectors_with_desired_magnitude = [vector * mag for vector in unit_vectors]
		u = [vector[0] for vector in vectors_with_desired_magnitude]  # x-components of unit vectors
		v = [vector[1] for vector in vectors_with_desired_magnitude]  # y-components of unit vectors

		return u, v

def find_event_hypocenter(slip_event, utm_crs):
	lat_hypo = slip_event['hypocenter_latitude_degrees']
	lon_hypo = slip_event['hypocenter_longitude_degrees']
	# select first
	lat_hypo = lat_hypo.iloc[0]
	lon_hypo = lon_hypo.iloc[0]
	source_crs = 'EPSG:4326'
	transformer = pyproj.Transformer.from_crs(source_crs, utm_crs, always_xy=True)
	utm_hypox, utm_hypoy = transformer.transform(lon_hypo, lat_hypo)
	return utm_hypox, utm_hypoy
