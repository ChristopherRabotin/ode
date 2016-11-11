package ode

import "testing"

func assertPanic(t *testing.T, msg string, f func()) {
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("code did not panic with %s", msg)
		}
	}()
	f()
}
