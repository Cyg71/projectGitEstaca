from orbitalPropagation import *

def test_Rot():
    """Rotational matrices test."""
    assert np.linalg.det(Rot(0, 0)) == 1
    assert np.linalg.det(Rot(0, 1)) == 1
    assert np.linalg.det(Rot(0, 2)) == 1
    id_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert np.array_equal(Rot(0, 3), id_matrix)

def test_radius():
    """Radius of ellipse test."""
    a = 6900e3
    assert radius(a, 0, 0) == a

def test_accelerations():
    """Acceleration and Jerk function test."""
    null_nd = np.array([0, 0, 0])
    acc_null, jerk_null = accelerations(np.array([0, 0, 0]), np.array([0, 0, 0]))
    acc = accelerations(np.array([RE, 0, 0]), np.array([0, 0, 0]))[0]
    acc[0] = round(acc[0], 6)
    assert np.array_equal(acc_null, null_nd)
    assert np.array_equal(jerk_null, null_nd)
    assert np.allclose(acc, np.array([round(-MU/RE**2, 6), 0, 0]), rtol = 1)

def test_kepler2state():
    """Test the kepler2state function for accuracy."""

    # Define test inputs
    a = 7000e3          # Semi-major axis [m]
    e = 0               # Eccentricity
    i = 90*m.pi/180     # Inclination [rad]
    Omega = 0           # Longitude of ascending node [rad]
    omega = 0           # Argument of perigee [rad]
    nu = 0              # True anomaly [rad]

    # Expected outputs (manually calculated or from a trusted tool)
    expected_r = np.array([a, 0, 0])  # Example position vector [m]
    expected_rdot = np.array([0, 0, m.sqrt(MU/a)])  # Example velocity vector [m/s]

    # Call the function
    r, rdot = kepler2state(a, e, i, Omega, omega, nu)

    # Assert results
    assert np.allclose(r, expected_r, atol = 1e-3), f"Position vector mismatch: {r} != {expected_r}"
    assert np.allclose(rdot, expected_rdot, atol = 1e-3), f"Velocity vector mismatch: {rdot} != {expected_rdot}"
