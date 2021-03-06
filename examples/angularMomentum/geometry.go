package dynamics

import "math"

// norm returns the norm of a given vector which is supposed to be 3x1.
func norm(v []float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// dot performs the inner product.
func dot(a, b []float64) float64 {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

// Spherical2Cartesian returns the provided spherical coordinates vector in Cartesian.
func Spherical2Cartesian(a []float64) (b []float64) {
	b = make([]float64, 3)
	sθ, cθ := math.Sincos(a[1])
	sφ, cφ := math.Sincos(a[2])
	b[0] = a[0] * sθ * cφ
	b[1] = a[0] * sθ * sφ
	b[2] = a[0] * cθ
	return
}

// Cartesian2Spherical returns the provided Cartesian coordinates vector in spherical.
func Cartesian2Spherical(a []float64) (b []float64) {
	b = make([]float64, 3)
	if norm(a) == 0 {
		return []float64{0, 0, 0}
	}
	b[0] = norm(a)
	b[1] = math.Acos(a[2] / b[0])
	b[2] = math.Atan2(a[1], a[0])
	return
}

// Deg2rad converts degrees to radians.
func Deg2rad(a float64) float64 {
	return a / 360.0 * 2 * math.Pi
}

// Rad2deg converts radians to degrees.
func Rad2deg(a float64) float64 {
	return a / (2 * math.Pi) * 360.0
}
