import math as m
import numpy as np
import matplotlib.pyplot as plt

# Physical Constants
MU = 3.986004418e14 						# [m**3/s**2] Gravitational parameter of the Earth
J2 = 1082.62668e-6 							# Perturbation constant from gravitational potential
RE = 6378.137e3 							# [m] Equatorial radius of the Earth
DAY_SECONDS = 86400 						# [s] Seconds in a solar day

def Rot(a, u):
	"""
	Rotation matrices method

	:param a: rotation angle [rad]
	:param u: rotation axis (0 = x, 1 = y, 2 = z)

	:returns: matrix 3x3
	"""
	c = np.cos(a)
	s = np.sin(a)
	if u == 0:
		Rot = np.array([[1, 0, 0], 
					[0, c, -s], 
					[0, s, c]])
	elif u == 1:
		Rot = np.array([[c, 0, s], 
					[0, 1, 0], 
					[-s, 0, c]])
	elif u == 2:
		Rot = np.array([[c, -s, 0], 
					[s, c, 0], 
					[0, 0, 1]])
	else:
		Rot = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		print("Error: wrong axis number provided for rotation matrix")

	return Rot


# Acceleration function
def accelerations(pos, vel):
	"""
	param: position & velocity
	return: acceleration (m/s²) & jerk: time deritative of acceleration (m/s3)
	"""
	pos_norm = np.linalg.norm(pos)
	rddot2B = - MU * pos / pos_norm**3

	rddotJ2 = [MU/pos_norm**2*(RE/pos_norm)**2*J2*(15/2*(pos[2]/pos_norm)**2 - 3/2)*pos[0]/pos_norm, 
			MU/pos_norm**2*(RE/pos_norm)**2*J2*(15/2*(pos[2]/pos_norm)**2 - 3/2)*pos[1]/pos_norm, 
			MU/pos_norm**2*(RE/pos_norm)**2*J2*(15/2*(pos[2]/pos_norm)**2 - 9/2)*pos[2]/pos_norm]
	
	acc = rddot2B + rddotJ2
	
	jerk = - MU * (vel/pos_norm**3 - 3*pos*np.dot(pos, vel)/pos_norm**5)
	return acc, jerk

def radius(a, ecc, nu):
	"""
	param : a,ecc,nu 

	a: semi major axis [m], defined as a float
	ecc: eccenntricity (=distance between focal point / half major axis), defined as a float
	nu: angle between the direction of the periapsis and the current position of an object [rad], defined as a float
	
	return: bare position radius for a non-circular orbit 
	"""
	return a * (1 - ecc**2)/(1 + ecc * m.cos(nu))

def kepler2state(a, e, i, Omega, omega, nu):
	"""
	param: a,e,i,Omega,omega,nu
	a: semi major axis [m], defined as a float
	ecc: eccenntricity (=distance between focal point / half major axis), defined as a float
	nu: angle between the direction of the periapsis and the current position of an object [rad], defined as a float
	i: Incidence, defined as a float [rad]
	omega: argument from perigee, is the position of perigee relative to ascending node, defined as a float [rad]

	return: r:position vector & rdot: velocity vector
	"""
	rc = radius(a, e, nu)
	E = 2*m.atan2(m.tan(nu/2), m.sqrt((1 + e)/(1 - e)))
	o = rc*np.array([m.cos(nu), m.sin(nu), 0])
	odot = m.sqrt(MU*a)/rc*np.array([-m.sin(E), m.sqrt(1 - e**2)*m.cos(E), 0])

	r = np.matmul(Rot(Omega, 2), np.matmul(Rot(i, 0), np.matmul(Rot(omega, 2), o)))
	rdot = np.matmul(Rot(Omega, 2), np.matmul(Rot(i, 0), np.matmul(Rot(omega, 2), odot)))
	return r, rdot

# Initial conditions
def tle2kepler(TLE):
	"""
	param : TLE 
	TLE : standardized representation of orbital parameters 

	return : orbital parameters : a,e,i,Omega,omega
	"""
	tle2 = TLE[1].split(' ')
	T = DAY_SECONDS/float(tle2[7])
	a = (MU*(T/(2*m.pi))**2)**(1/3)

	e = float(str(0.) + tle2[4])
	i = float(tle2[2]) * m.pi/180
	Omega = float(tle2[3]) * m.pi/180
	omega = float(tle2[5]) * m.pi/180
	nu = float(tle2[6]) * m.pi/180

	return a, e, i, Omega, omega, nu

TLE = ['1 55044U 23001AM 24318.44334776 .00176663 00000-0 20148-2 0 9996', 
	   '2 55044 97.4024 23.3159 0005184 341.4271 18.6796 15.616308481 3798']

a, e, i, Omega, omega, nu = tle2kepler(TLE)
r, rdot = kepler2state(a, e, i, Omega, omega, nu)
rddot, jerk = accelerations(r, rdot)

# Simulation parameters
time = DAY_SECONDS*1  # total time
dt = 10    # time step

#Here is a test to see whether the problem is fixed

# Initialize arrays to store position and velocity
r_array = []
rdot_array = []

# Main simulation loop
for i in range(int(time / dt)):
    r_new = r + rdot * dt + rddot * dt**2 / 2 + jerk * dt**3 / 6
    rdot_new = rdot + rddot * dt + jerk * dt**2 / 2
    rddot, jerk = accelerations(r_new, rdot_new)
    r_array.append(r_new)
    rdot_array.append(rdot_new)
    r = r_new
    rdot = rdot_new

# Convert lists to numpy arrays
r_array = np.array(r_array)
rdot_array = np.array(rdot_array)

# Plotting the orbit in 3D
fig2 = plt.figure()
plt.title("Orbit in ECI")

# Creating a wireframe for the Earth
phi = np.linspace(0, 2 * m.pi, 36)
theta = np.linspace(0, m.pi, 18)
A, B = np.meshgrid(theta, phi)
X = RE * np.sin(A) * np.cos(B)  # X coordinates
Y = RE * np.sin(A) * np.sin(B)  # Y coordinates
Z = RE * np.cos(A)              # Z coordinates

# Plotting the Earth and the trajectory
eci = fig2.add_subplot(111, projection = '3d')
eci.plot_wireframe(X, Y, Z, color = 'c', zorder = 1, alpha = 0.25)  # Sphere
eci.plot3D(r_array[:, 0], r_array[:, 1], r_array[:, 2], color = 'm')  # Trajectory

# Plotting the ECI coordinates in 2D
plt.figure()
plt.plot(r_array[:, 0], color = 'r')	# X_ECI
plt.plot(r_array[:, 1], color = 'g')	# Y_ECI
plt.plot(r_array[:, 2], color = 'b')	# Z_ECI
plt.grid()
plt.legend(['$X_{ECI}$', '$Y_{ECI}$', '$Z_{ECI}$'])
plt.title('Coordinates in ECI of the satellite versus time')
plt.xlabel(f"Time, one step = {dt} sec")
plt.ylabel('ECI coordinates in m')
plt.show()
