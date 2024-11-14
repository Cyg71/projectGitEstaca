import math as m
import numpy as np
import matplotlib.pyplot as plt

# Physical Constants
mu = 3.986004418e14 						# [m**3/s**2] gravitational parameter of the Earth
J2 = 1082.62668e-6 							# Perturbation constant from gravitational potential
Re = 6378.137e3 							# [m] equatorial radius of the Earth
Rp = 6356.752e3 							# [m] polar radius of the Earth
day_seconds = 86400 						# [s] seconds in a solar day
sidday_seconds = 86164 						# [s] seconds in a sidereal day
g0 = 9.80665								# [m/s2] Intensity of gravity at mean sea level

def Rot(a):
	c = np.cos(a)
	s = np.sin(a)
	Rotx = np.array([[1, 0, 0], 
				  [0, c, -s], 
        		  [0, s, c]])

	Roty = np.array([[c, 0, s], 
		          [0, 1, 0], 
        		  [-s, 0, c]])

	Rotz = np.array([[c, -s, 0], 
		          [s, c, 0], 
	      		  [0, 0, 1]])

	return Rotx, Roty, Rotz

def Mult4(a, b, c, d):
	return np.matmul(a, np.matmul(b, np.matmul(c, d)))

# Acceleration function
def accelerations(pos, vel):
	pos_norm = np.linalg.norm(pos)
	rddot2B = - mu * pos / pos_norm**3
	#rddotJ2 =  TO DO acceleration du J2
	xdd = - mu/pos_norm**2*(Re/pos_norm)**2*(-J2)*(15/2*(pos[2]/pos_norm)**2-3/2)*pos[0]/pos_norm
	ydd = - mu/pos_norm**2*(Re/pos_norm)**2*(-J2)*(15/2*(pos[2]/pos_norm)**2-3/2)*pos[1]/pos_norm
	zdd = - mu/pos_norm**2*(Re/pos_norm)**2*(-J2)*(15/2*(pos[2]/pos_norm)**2-9/2)*pos[2]/pos_norm
	rddotJ2 = [xdd, ydd, zdd]
	acc = rddot2B + rddotJ2
	
	jerk = - mu * (vel/pos_norm**3 - 3*pos*np.dot(pos, vel)/pos_norm**5)
	return acc, jerk

def radius(a, ecc, nu):
	return a * (1 - ecc**2)/(1 + ecc * m.cos(nu))

def kepler2state(a, e, i, Omega, omega, nu):
	rc = radius(a, e, nu)
	E = 2*m.atan2(m.tan(nu/2), m.sqrt((1+e)/(1-e)))
	o = rc*np.array([m.cos(nu), m.sin(nu), 0])
	odot = m.sqrt(mu*a)/rc*np.array([-m.sin(E), m.sqrt(1-e**2)*m.cos(E), 0])
	r = Mult4(Rot(Omega)[2], Rot(i)[0], Rot(omega)[2], o)
	rdot = Mult4(Rot(Omega)[2], Rot(i)[0], Rot(omega)[2], odot)
	return r, rdot

# Initial conditions
def tle2kepler(TLE):
    tle2 = TLE[1].split(' ')
    T = day_seconds/float(tle2[7])
    a = (mu*(T/(2*m.pi))**2)**(1/3)
	
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
time = day_seconds*1  # total time
dt = 10    # time step

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

# Convert lists to NumPy arrays for easy slicing and plotting
r_array = np.array(r_array)
rdot_array = np.array(rdot_array)

# Plotting the orbit in 3D
fig2 = plt.figure()
plt.title("Orbit in ECI")

# Creating a wireframe for the Earth
phi = np.linspace(0, 2 * m.pi, 360)
theta = np.linspace(0, m.pi, 180)
A, B = np.meshgrid(theta, phi)
X = Re * np.sin(A) * np.cos(B)  # X coordinates
Y = Re * np.sin(A) * np.sin(B)  # Y coordinates
Z = Re * np.cos(A)              # Z coordinates

# Plotting the Earth and the trajectory
eci = fig2.add_subplot(111, projection='3d')
eci.plot_wireframe(X, Y, Z, color='c', zorder=1, alpha=0.25)  # Sphere
eci.plot3D(r_array[:, 0], r_array[:, 1], r_array[:, 2], color='m')  # Trajectory
plt.show()
