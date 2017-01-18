package ode

import (
	"fmt"
	"math"
	"testing"

	"github.com/ChristopherRabotin/ode/examples/angularMomentum"
	"github.com/gonum/floats"
)

func TestPanics(t *testing.T) {
	assertPanic(t, "negative step", func() {
		NewRK4(1, -1, nil)
	})

	assertPanic(t, "nil integrator", func() {
		NewRK4(1, 1, nil)
	})
}

type Balbasi1D struct {
	state  []float64 // Note that we don't have a state history here.
	prevIt uint
}

func NewBalbasi1D() (b *Balbasi1D) {
	b = &Balbasi1D{}
	b.state = []float64{1200.0}
	b.prevIt = 0
	return
}

func (b *Balbasi1D) GetState() []float64 {
	return b.state
}

func (b *Balbasi1D) SetState(i uint64, s []float64) {
	b.state = s
	if i != 0 && b.prevIt+1 != uint(i) {
		panic(fmt.Errorf("expected i=%d, got i=%d", b.prevIt+1, i))
	}
	b.prevIt = uint(i)
}

func (b *Balbasi1D) Stop(i uint64) bool {
	return i*30 >= 480
}

func (b *Balbasi1D) Func(t float64, s []float64) []float64 {
	val := []float64{(-2.2067 * 1e-12) * (math.Pow(s[0], 4) - 81*1e8)}
	return val
}

func TestRK4In1D(t *testing.T) {
	inte := NewBalbasi1D()
	if _, _, err := NewRK4(1, 30, inte).Solve(); err != nil {
		t.Fatalf("err: %+v\n", err)
	}
	if diff := math.Abs(inte.GetState()[0] - 647.5720536920); diff >= 1e-10 {
		t.Fatalf("expected final state %4.10f is different by %4.10f", inte.GetState()[0], diff)
	}
}

// AttitudeTest tests that the energy is conserved for a given body.
type AttitudeTest struct {
	*dynamics.Attitude
}

func (a *AttitudeTest) Stop(i uint64) bool {
	return float64(i)*1e-6 >= 1e-1
}

func NewAttitudeTest() (a *AttitudeTest) {
	a = &AttitudeTest{dynamics.NewAttitude([3]float64{0.3, -0.4, 0.5}, [3]float64{0.1, 0.4, -0.2},
		[]float64{10.0, 0, 0, 0, 5.0, 0, 0, 0, 2.0})}
	return
}

func TestRK4Attitude(t *testing.T) {
	inte := NewAttitudeTest()
	initMom := inte.Momentum()
	for _, step := range []float64{1e-4, 1e-6, 1e-8, 1e-10, 1e-12} {
		if _, _, err := NewRK4(0, step, inte).Solve(); err != nil {
			t.Fatalf("err: %+v\n", err)
		}
		if diff := math.Abs(initMom - inte.Momentum()); diff > 1e-8 {
			t.Fatalf("angular momentum changed by %4.12f", diff)
		}
	}
}

type VSimple struct {
	state []float64
}

func (v *VSimple) GetState() []float64 {
	return v.state
}

func (v *VSimple) SetState(i uint64, s []float64) {
	v.state = s
}

func (v *VSimple) Stop(i uint64) bool {
	return i >= 189
}

func (v *VSimple) Func(x float64, y []float64) []float64 {
	return []float64{y[1], -y[0]}
}

func TestRK4Simple(t *testing.T) {
	inte := new(VSimple)
	inte.state = []float64{0, 1}
	iterNum, xi, err := NewRK4(0, 0.2, inte).Solve()
	if err != nil {
		t.Fatalf("err: %+v\n", err)
	}
	if xi != 37.8 {
		t.Fatalf("xi=%f != 37.8", xi)
	}
	if iterNum != 189 {
		t.Fatalf("iterNum=%d != 189", iterNum)
	}
	exp := []float64{+1.0021441571397413e-01, +9.9488186473553231e-01}
	state := inte.GetState()
	if exp[0] != state[0] || exp[1] != state[1] {
		t.Fatalf("\nstate=%+v\n  exp=%+v", state, exp)
	}
}

type KraichnanOrszag struct {
	steps uint64
	state []float64
}

func (v *KraichnanOrszag) GetState() []float64 {
	return v.state
}

func (v *KraichnanOrszag) SetState(i uint64, s []float64) {
	v.state = s
}

func (v *KraichnanOrszag) Stop(i uint64) bool {
	return i >= v.steps
}

func (v *KraichnanOrszag) Func(x float64, y []float64) []float64 {
	return []float64{y[0] * y[2], -y[1] * y[2], -y[0]*y[0] + y[1]*y[1]}
}

func TestRK4KO(t *testing.T) {
	tolerance := 1e-10
	for _, steps := range []uint64{30, 3000} {
		inte := &KraichnanOrszag{steps, []float64{1, .4, .2}}
		iterNum, xi, err := NewRK4(0, 0.01, inte).Solve()
		if err != nil {
			t.Fatalf("err: %+v\n", err)
		}
		if iterNum != steps {
			t.Fatalf("iterNum=%d != %d", iterNum, steps)
		}
		var exp []float64
		if steps == 30 {
			if !floats.EqualWithinAbs(xi, 0.3, tolerance) {
				t.Fatalf("xi=%.16f != 0.3", xi)
			}
			exp = []float64{1.0209861554390987e+00, 3.9177808423412647e-01, -6.4009398764259762e-02}
		} else {
			if !floats.EqualWithinAbs(xi, 30, tolerance) {
				t.Fatalf("xi=%.16f != 30", xi)
			}
			exp = []float64{4.8745696934931565e-01, 8.2058525186654796e-01, -5.3761096276439480e-01}
		}
		state := inte.GetState()
		if !floats.EqualApprox(exp, state, tolerance) {
			t.Fatalf("\nstate=%+v\n  exp=%+v", state, exp)
		}
	}
}
