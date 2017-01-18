package ode

// RK4 defines an RK4 integrator using math.Big Floats.
type RK4 struct {
	X0        float64    // The initial x0.
	StepSize  float64    // The step size.
	Integator Integrable // What is to be integrated.
}

// NewRK4 returns a new RK4 integrator instance.
func NewRK4(x0 float64, stepSize float64, inte Integrable) (r *RK4) {
	if stepSize <= 0 {
		panic("config StepSize must be positive")
	}
	if inte == nil {
		panic("config Integator may not be nil")
	}
	r = &RK4{X0: x0, StepSize: stepSize, Integator: inte}
	return
}

// Solve solves the configured RK4.
// Returns the number of iterations performed and the last X_i, or an error.
func (r *RK4) Solve() (uint64, float64, error) {
	iterNum := uint64(0)
	xi := r.X0
	for !r.Integator.Stop(iterNum) {
		halfStep := r.StepSize * 0.5
		state := r.Integator.GetState()
		newState := make([]float64, len(state))
		z := make([]float64, len(state)) // a temporary variable

		// Step 1
		f1 := r.Integator.Func(xi, state)

		// Step 2
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + halfStep*f1[i]
		}
		f2 := r.Integator.Func(xi+halfStep, z)

		// Step 3
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + halfStep*f2[i]
		}
		f3 := r.Integator.Func(xi+halfStep, z)

		// Step 4
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + r.StepSize*f3[i]
		}
		f4 := r.Integator.Func(xi+r.StepSize, z)

		for i := 0; i < len(state); i++ {
			newState[i] = state[i] + r.StepSize*(f1[i]+2*f2[i]+2*f3[i]+f4[i])/6
		}
		xi += r.StepSize
		r.Integator.SetState(iterNum, newState)
		iterNum++ // Don't forget to increment the number of iterations.
	}

	return iterNum, xi, nil
}
