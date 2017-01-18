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
	const (
		half     = 1 / 2.0
		oneSixth = 1 / 6.0
		oneThird = 1 / 3.0
	)

	iterNum := uint64(0)
	xi := r.X0
	for !r.Integator.Stop(iterNum) {
		//halfStep := xi * half // WHAT?!
		state := r.Integator.GetState()
		newState := make([]float64, len(state))
		z := make([]float64, len(state))
		//k1, k2, k3, k4 are used as buffers AND result variables.
		/*k1 := make([]float64, len(state))
		k2 := make([]float64, len(state))
		k3 := make([]float64, len(state))
		k4 := make([]float64, len(state))
		tState := make([]float64, len(state))

		// Compute the k's.
		k1 := r.Integator.Func(xi, state)
		for i, y := range k1 {
			k1[i] = y * r.StepSize
			tState[i] = state[i] + k1[i]*half
		}
		k2 := r.Integator.Func(xi+halfStep, tState)
		for i, y := range k2 {
			k2[i] = y * r.StepSize
			tState[i] = state[i] + k2[i]*half
		}
		k3 := r.Integator.Func(xi+halfStep, tState)
		for i, y := range k3 {
			k3[i] = y * r.StepSize
			tState[i] = state[i] + k3[i]
		}
		k4 := r.Integator.Func(xi+halfStep, tState)
		for i, y := range k4 {
			k4[i] = y * r.StepSize
			newState[i] = state[i] + oneSixth*(k1[i]+k4[i]) + oneThird*(k2[i]+k3[i])
		}
		r.Integator.SetState(iterNum, newState)*/

		// Step 1
		//dydx(x, y, f1)
		f1 := r.Integator.Func(xi, state)

		// Step 2
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + r.StepSize*f1[i]/2
		}
		//dydx(x+h/2, z, f2)
		f2 := r.Integator.Func(xi+r.StepSize/2, z)

		// Step 3
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + r.StepSize*f2[i]/2
		}
		//dydx(x+r.StepSize/2, z, f3)
		f3 := r.Integator.Func(xi+r.StepSize/2, z)

		// Step 4
		for i := 0; i < len(state); i++ {
			z[i] = state[i] + r.StepSize*f3[i]
		}
		//dydx(x+r.StepSize, z, f4)
		f4 := r.Integator.Func(xi+r.StepSize, z)

		//ynew := ys[k*nd:]
		for i := 0; i < len(state); i++ {
			newState[i] = state[i] + r.StepSize*(f1[i]+2*f2[i]+2*f3[i]+f4[i])/6
		}
		r.Integator.SetState(iterNum, newState)
		//x += h
		//y = ynew

		xi += r.StepSize
		iterNum++ // Don't forget to increment the number of iterations.
	}

	return iterNum, xi, nil
}
