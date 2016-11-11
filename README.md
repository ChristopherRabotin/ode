# ode
An [ordinary differential equation](https://en.wikipedia.org/wiki/Ordinary_differential_equation) solving library in golang.

# Features
+ Multi-dimensional state vector (i.e. extended states)
+ Channel based stopping condition, which takes full advantage of [goroutines](https://en.wikipedia.org/wiki/Go_(programming_language#Concurrency)
+ Trivial usage (cf. the examples directory)
+ Easily expendable to other ODE solving methods

# Numerical methods
- [x] [RK4](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method)
- [ ] [Dormant-Prince](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method)

# Concurrency example
*The following comes from a closed source project.*

In `main.go`:
```golang
astro, wg := dynamics.NewAstro(sc, orbit, start, end, name)
go astro.Propagate()
```
