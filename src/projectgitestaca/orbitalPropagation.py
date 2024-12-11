import math as m
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt

# Physical Constants
MU = 3.986004418e14 						# [m**3/s**2] Gravitational parameter of the Earth
J2 = 1082.62668e-6 							# [-] Perturbation constant from gravitational potential
RE = 6378.137e3 							# [m] Equatorial radius of the Earth
DAY_SECONDS = 86400 						# [s] Seconds in a solar day

def Rot(a: float, u: int)->NDArray:	
	"""
	Rotation matrices method

	Args:
		a: rotation angle [rad], float
		u: rotation axis (0 = x, 1 = y, 2 = z), int

	Returns:
		Rot: 3x3 matrix, NDArray
	"""
	c = np.cos(a)
	s = np.sin(a)
	if u == 0:						# Rotation about the x-axis
		Rot = np.array([[1, 0, 0], 
						[0, c, -s], 
						[0, s, c]])
	elif u == 1:					# Rotation about the y-axis
		Rot = np.array([[c, 0, s], 
						[0, 1, 0], 
						[-s, 0, c]])
	elif u == 2:					# Rotation about the z-axis
		Rot = np.array([[c, -s, 0], 
						[s, c, 0], 
						[0, 0, 1]])
	else:							# Wrong axis number: error
		Rot = np.array([[1, 0, 0],
						[0, 1, 0], 
						[0, 0, 1]])
		print("Error: wrong axis number provided for rotation matrix")
	return Rot

# Acceleration function
def accelerations(pos: NDArray, vel: NDArray) -> tuple[NDArray, NDArray]:
    """
    Calculates the instantaneous acceleration and jerk vectors

	Args:
    	pos: ECI position vector, NDArray 3x1 
    	vel: ECI velocity vector, NDArray 3x1 

    Returns: 
		ECI acceleration vector 3x1 (m/s2), NDArray
		ECI jerk vector 3x1 first time deritative of acceleration (m/s3), NDArray
    """
    pos_norm = np.linalg.norm(pos)
    
    if pos_norm == 0:
        return np.zeros(3), np.zeros(3)  # Special case : position at the origin of the reference frame

    rddot2B = - MU * pos / pos_norm**3

    rddotJ2 = [
        MU / pos_norm**2 * (RE / pos_norm)**2 * J2 * (15 / 2 * (pos[2] / pos_norm)**2 - 3 / 2) * pos[0] / pos_norm,
        MU / pos_norm**2 * (RE / pos_norm)**2 * J2 * (15 / 2 * (pos[2] / pos_norm)**2 - 3 / 2) * pos[1] / pos_norm,
        MU / pos_norm**2 * (RE / pos_norm)**2 * J2 * (15 / 2 * (pos[2] / pos_norm)**2 - 9 / 2) * pos[2] / pos_norm,
    ]

    acc = rddot2B + rddotJ2
    jerk = - MU * (vel / pos_norm**3 - 3 * pos * np.dot(pos, vel) / pos_norm**5)
    return acc, jerk

def radius(a: float, ecc: float, nu: float)->float:
	"""
	Instantaneous radius of the orbit, given for the current Keplerian elements a, ecc and nu

	Args:
		a: semi major axis [m], float
		ecc: eccentricity (distance between focal point / half major axis), float
		nu: angle between the direction of the periapsis and the current position of an object [rad], float
	
	Returns: 
		Position radius for a non-circular orbit, float
	"""
	return a * (1 - ecc**2)/(1 + ecc * m.cos(nu))

def kepler2state(a: float, ecc: float, i: float, Omega: float, omega: float, nu: float)->tuple[NDArray,NDArray]:
	"""
	Converts the Keplerian elements into state vectors

	Args:
		a: semi major axis [m], float
		ecc: eccentricity (=distance between focal point / half major axis), float
		i: inclination of the orbit [rad], float
		Omega: longitude of the right ascension of the ascending node [rad], float
		omega: argument of perigee [rad], float
		nu: angle between the direction of the periapsis and the current position of an object [rad], float

	Returns: 
		r: ECI position vector 3x1, NDArray
		rdot: ECI velocity vector 3x1, NDArray
	"""
	rc = radius(a, ecc, nu)
	E = 2*m.atan2(m.tan(nu/2), m.sqrt((1 + ecc)/(1 - ecc)))
	o = rc*np.array([m.cos(nu), m.sin(nu), 0])
	odot = m.sqrt(MU*a)/rc*np.array([-m.sin(E), m.sqrt(1 - ecc**2)*m.cos(E), 0])

	r = np.matmul(Rot(Omega, 2), np.matmul(Rot(i, 0), np.matmul(Rot(omega, 2), o)))
	rdot = np.matmul(Rot(Omega, 2), np.matmul(Rot(i, 0), np.matmul(Rot(omega, 2), odot)))
	return r, rdot

# Initial conditions - Keplerian elements
TLE = ['1 55044U 23001AM 24318.44334776 .00176663 00000-0 20148-2 0 9996', 
	'2 55044 97.4024 23.3159 0005184 341.4271 18.6796 15.616308481 3798']

# Simulation parameters
nb_days = 1						# [day] Number of days of propagation
time = DAY_SECONDS*nb_days  	# [s] Total time
dt = 5    						# [s] Time step of the simulation

if __name__ == '__main__': # pragma: no cover

	tle2 = TLE[1].split(' ')
	T = DAY_SECONDS/float(tle2[7])		# [s] Orbital period
	a = (MU*(T/(2*m.pi))**2)**(1/3)		# [m] Semi-major axis, using Kepler's third law
	e = float(str(0.) + tle2[4])		# [-] Eccentricity
	i = float(tle2[2]) * m.pi/180		# [rad] Inclination
	Omega = float(tle2[3]) * m.pi/180	# [rad] Longitude of the Right Ascension of the Ascending Node (RAAN)
	omega = float(tle2[5]) * m.pi/180	# [rad] Argument of perigee
	nu = float(tle2[6]) * m.pi/180		# [rad] True anomaly

	# Initial conditions - State vectors
	r, rdot = kepler2state(a, e, i, Omega, omega, nu)	# [m, m/s]
	rddot, jerk = accelerations(r, rdot)				# [m/s2, m/s3]

	# Initialize arrays to store position and velocity
	r_array = []
	rdot_array = []

	# Main simulation loop
	for i in range(int(time / dt)):
		r_new = r + rdot * dt + rddot * dt**2 / 2 + jerk * dt**3 / 6	# [m] Forward Euler, order 3 integration
		rdot_new = rdot + rddot * dt + jerk * dt**2 / 2					# [m/s] Forward Euler, order 2 integration
		rddot, jerk = accelerations(r_new, rdot_new)					# [m/s2, m/s3] New acceleration and jerk vectors
		r_array.append(r_new)											# [m] Position array
		rdot_array.append(rdot_new)										# [m/s] Velocity array
		r = r_new														# [m] Resetting position vector
		rdot = rdot_new													# [m/s] Resetting velocity vector

	# Convert lists to numpy arrays
	r_array = np.array(r_array)
	rdot_array = np.array(rdot_array)

	# Plotting the orbit in 3D
	fig2 = plt.figure()
	plt.title("Orbit in ECI")

	# Creating a wireframe for the Earth
	phi = np.linspace(0, 2 * m.pi, 36)	# Angular meshing
	theta = np.linspace(0, m.pi, 18)	# Angular meshing
	A, B = np.meshgrid(theta, phi)
	X = RE * np.sin(A) * np.cos(B)  # X coordinates
	Y = RE * np.sin(A) * np.sin(B)  # Y coordinates
	Z = RE * np.cos(A)              # Z coordinates

	# Plotting the Earth and the trajectory
	eci = fig2.add_subplot(111, projection = '3d')
	# Sphere
	eci.plot_wireframe(X, Y, Z, color = 'c', zorder = 1, alpha = 0.25)  # type: ignore
	# Trajectory
	eci.plot3D(r_array[:, 0], r_array[:, 1], r_array[:, 2], color = 'm') # type: ignore
	
	# Plotting the ECI coordinates in 2D
	plt.figure()
	# X_ECI
	plt.plot(r_array[:, 0], color = 'r') # type: ignore
	# Y_ECI 
	plt.plot(r_array[:, 1], color = 'g') # type: ignore
	# Z_ECI
	plt.plot(r_array[:, 2], color = 'b') # type: ignore
	plt.grid()
	plt.legend(['$X_{ECI}$', '$Y_{ECI}$', '$Z_{ECI}$'])
	plt.title('Coordinates in ECI of the satellite versus time')
	plt.xlabel(f"Time, one step = {dt} sec")
	plt.ylabel('ECI coordinates in m')
	plt.show()
