import matplotlib.colors as colors
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import pyproj
import numpy as np
from mpmath import sec
from mpmath import *
from itertools import chain
from scipy import stats

# from stack exchange
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
	"""
	Find event hypocenter and convert to UTM coordinates of event coordinate system
	"""	
	lat_hypo = slip_event['hypocenter_latitude_degrees']
	lon_hypo = slip_event['hypocenter_longitude_degrees']
	# select first
	lat_hypo = lat_hypo.iloc[0]
	lon_hypo = lon_hypo.iloc[0]
	source_crs = 'EPSG:4326'
	transformer = pyproj.Transformer.from_crs(source_crs, utm_crs, always_xy=True)
	utm_hypox, utm_hypoy = transformer.transform(lon_hypo, lat_hypo)
	return utm_hypox, utm_hypoy

def measure_fault_strike(x_coords,y_coords):
	"""
	Measure the strike of the fault segment
	"""
	angles = []
	for i in range(1, len(x_coords)):
		angle_radians = math.atan2(y_coords[i] - y_coords[i - 1], x_coords[i] - x_coords[i - 1])
		angle_degrees = math.degrees(angle_radians) + 90
		# Ensure the angle is between -90 and 90 degrees
		if angle_degrees < -90:
			angle_degrees += 180
		elif angle_degrees > 90:
			angle_degrees -= 180
		angles.append(angle_degrees)
	return angles


def measure_fracture_angle(x_coords, y_coords, angles, cropped_xstress_utm, cropped_ystress_utm, cropped_SHmax):
	"""
	Measure the angle between shmax and the fault segment strike
	"""
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	# x_coords = np.array(x_coords)
	# y_coords = np.array(y_coords)

	cropped_xstress_utm = np.array(cropped_xstress_utm)
	cropped_ystress_utm = np.array(cropped_ystress_utm)
	diff_angle = []
	SHmax_segment = []

	for x, y, angle in zip(x_coords, y_coords, angles):
		# Calculate the Euclidean distance to all points in cropped_xstress_utm and cropped_ystress_utm
		distances = np.sqrt((x - cropped_xstress_utm)**2 + (y - cropped_ystress_utm)**2)
		# Find the index of the minimum distance
		nearest_index = np.argmin(distances)
		SHmax_min_loc = cropped_SHmax[nearest_index]
		SHmax_segment.append(SHmax_min_loc)
		diff_angle.append(abs(SHmax_min_loc - angle))
	return diff_angle, SHmax_segment

def plot_segments_color_fracture_angle(vmin, vmax, mid_val, x_coords, y_coords, diff_angle, eventname):
	"""
	Plot segments color-coded by fracture angle (shmax-fault segment strike)
	"""
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	# Create a scatter plot of the extracted coordinates with colors based on angles in degrees
	plt.scatter(x_coords, y_coords, c=diff_angle, cmap='RdBu_r', marker='o', alpha=0.7, edgecolors='k', clim=(vmin, vmax), norm=MidpointNormalize(midpoint=mid_val, vmin=vmin, vmax=vmax))
	# Set the desired range for the colorbar
	plt.title(eventname)
	plt.axis('equal')
	plt.colorbar(label=r'Fracture angle $\phi_s$ (degrees)')

def plot_angle_distribution(subplot_ID,angles,nbins,bincol,xlabel):
	"""
	Plot histograms of measured angles
	"""
	plt.subplot(subplot_ID)  # 2 rows, 1 column, first subplot
	plt.hist(angles, bins=nbins, color=bincol, edgecolor='none',alpha=0.5)
	mean_angle = np.mean(angles)
	std_angle = np.std(angles) 
	plt.axvline(mean_angle, color=bincol, linestyle='solid', linewidth=1)
	plt.axvline(mean_angle+std_angle, color=bincol, linestyle='dashed', linewidth=1)
	plt.axvline(mean_angle-std_angle, color=bincol, linestyle='dashed', linewidth=1)
	plt.xlabel(xlabel)
	plt.ylabel('Frequency')

def prepare_slip_data(slip_event,utm_crs):
	"""
	Extracts slip data from dataframe and converts to UTM coordinates of event coordinate system
	"""
	# now load slip data 
	slip_lat = slip_event['latitude_degrees']
	slip_lon = slip_event['longitude_degrees']
	slip_preferred = abs(slip_event['recommended_net_preferred_for_analysis_meters'])
	rake_event = slip_event['fps_style']
	indices = (slip_preferred > 0) & (slip_preferred < 20) # removing zeros and errors - val 999
	slip_preferred = slip_preferred[indices]
	slip_lat = slip_lat[indices]
	slip_lon = slip_lon[indices]
	rake = rake_event[indices]
	slipx, slipy = extract_coordinates_stress(utm_crs, slip_lat, slip_lon)
	return slipx, slipy, slip_preferred, rake


def match_slip_to_segment(x_coords, y_coords, xslip_utm, yslip_utm, slip_preferred, rake, num_closest=10):
	"""
	Find closest slip values (number of values considered given by variable num_closts) to each 
	fault segment vertex, then selects the maximum from those values 
	"""
	max_slip_preferred = []
	rake_slip = []
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	for x, y in zip(x_coords, y_coords):
		# Calculate the Euclidean distances to all points in xslip_utm and yslip_utm
		distances = np.sqrt((x - xslip_utm)**2 + (y - yslip_utm)**2)          
		# Find the indices of the num_closest smallest distances
		closest_indices = np.argsort(distances)[:num_closest]   
		# Pick the largest slip_preferred value within those indices
		max_slip_preferredi = np.max(slip_preferred.iloc[closest_indices])
		rake_preferredi = stats.mode(rake[closest_indices])
		rake_preferredi = rake_preferredi[0]
		max_slip_preferred.append(max_slip_preferredi)
		rake_slip.append(rake_preferredi)
	return max_slip_preferred, rake_slip

def convert_rake_to_slip_vector(rake,strike):    
    slip_vector = np.array(strike)+np.array(rake)
    return slip_vector

def measure_slip_vector_angle(slip_vector,shmax):
    slip_vector_angle = slip_vector - shmax
    return slip_vector_angle

def match_slip_to_epicenter(x_coords, y_coords, epix, epiy, slip):
	"""
	Find closest slip value to event epicenter
	"""
	slip = slip.tolist()
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	distances = []
	for x, y in zip(x_coords, y_coords):
		distances.append(np.sqrt((x - epix)**2 + (y - epiy)**2))
	closest_index = np.argsort(distances)[0]  
	distance = distances[closest_index] # could export if interested
	slip_epi = slip[closest_index]
	return slip_epi

def match_frac_angle_to_epicenter(x_coords, y_coords, epix, epiy, frac_angle):
	"""
	Find closest fracture angle to event epicenter
	"""
	frac_angle =  list(chain.from_iterable(frac_angle))
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	distances = []
	for x, y in zip(x_coords, y_coords):
		distances.append(np.sqrt((x - epix)**2 + (y - epiy)**2))
	closest_index = np.argsort(distances)[0]  
	distance = distances[closest_index] # could export if interested
	frac_angle_epi = frac_angle[closest_index]
	return frac_angle_epi

def plot_segments_slip_color(x_coords,y_coords,max_slip_preferred,eventname,vmin,vmax):
	"""
	Plot segments color-coded by slip value (log10)
	"""
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]
	# Create a scatter plot of the extracted coordinates with colors based on angles in degrees
	log_slip_preferred = np.log10(max_slip_preferred)
	plt.scatter(x_coords, y_coords, c=log_slip_preferred, cmap='magma', marker='o', alpha=1, edgecolors='k', clim=(-2, 0.8))
	# Set the desired range for the colorbar
	plt.title(eventname)
	plt.axis('equal')
	plt.colorbar(label=r'log Slip (m)')	

def calculate_strain_drop(x_coords,y_coords,slip):
	"""
	Measure the strain drop for every segment as the diference in slip across the segment divided by the segment length
	"""
	strain_drop = []
	x_coords = x_coords[1:]
	y_coords = y_coords[1:]

	for i in range(1, len(x_coords)):
		segment_length = np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
		slip_diff = abs(slip[i] - slip[i-1])
		strain_drop.append(slip_diff / segment_length)
	return strain_drop

def calculate_shear_stress(radius,frac_angle):
	"""
	Estimate the shear stress for every segment 
	"""
	frac_angle = 2*(90-frac_angle)
	frac_angle = np.deg2rad(frac_angle)
	tau = radius * np.sin(frac_angle)
	return abs(tau)

def calculate_normal_stress(radius,frac_angle,shear_stress):
	"""
	Estimate the normal stress for every segment 
	"""
	frac_angle = 2*(90-frac_angle)
	frac_angle = np.deg2rad((frac_angle))
	center_minus_sigma = shear_stress/np.tan(frac_angle) #radius * np.cos(frac_angle)
	return center_minus_sigma


# def plot_instability(normal):
# 	minx = np.min(normal)
# 	maxx = 	np.max(normal)
# 	xvec = np.linspace(minx,maxx,100)
# 	plt.plot(xvec,xvec,'k--')
	
# # estimate delta sigma, or the maximum shear stress, which is the radius of the Morh circle for strike-slip faults
# def estimate_max_shear_stress(rot, shmax, delta_tau):
# 	shmax = np.deg2rad(shmax)
# 	rot = np.deg2rad(rot)
# 	max_shear_stress = (1/2) * delta_tau * np.cos(shmax) * (-np.cos(2*shmax)**2 * (1/np.tan(rot)) * sec(shmax)**2 + 4*np.sin(shmax) + np.tan(rot)) # from Hardebeck & Hauksson (2001) and Milliner et al. (2022) 
# 	max_shear_stress = float(max_shear_stress)
# 	return abs(max_shear_stress)

# # initial shear stress on failure plane
# def estimate_tau_initial(max_shear_stress, phi_o):
# 	phi_o = float(phi_o)
# 	phi_o = np.deg2rad(2*phi_o)
# 	tau_initial = max_shear_stress * np.sin(phi_o)
# 	return abs(tau_initial)

# def estimate_static_friction(phi_o):
# 	phi_o = np.deg2rad(phi_o)
# 	static_friction = np.tan(2*((np.pi/4)-phi_o))
# 	return abs(static_friction)

# def estimate_effective_normal_stress_initial(tau_initial,static_friction):
# 	effective_normal_stress_initial = tau_initial / static_friction
# 	return abs(effective_normal_stress_initial)

# def estimate_mean_horizontal_stress(effective_normal_stress,max_shear_stress,frac_angle):
# 	frac_angle = np.deg2rad(frac_angle*2)
# 	P = effective_normal_stress + max_shear_stress*np.cos(frac_angle)
# 	return abs(P)

# def estimate_mean_horizontal_stressb(max_shear_stress,static_friction):
# 	fric_angle = 1/np.tan(static_friction)
# 	P = max_shear_stress/fric_angle
# 	return abs(P)

# def estimate_principal_stresses(P,max_shear_stress):
# 	sigma_one = P + max_shear_stress
# 	sigma_three = P - max_shear_stress
# 	return sigma_one, sigma_three

# def estimate_dynamic_friction(stress_drop,effective_normal_stress,static_friction):
# 	dynamic_friction = static_friction - (stress_drop/effective_normal_stress)
# 	return dynamic_friction

