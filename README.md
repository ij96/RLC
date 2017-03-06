# RLC
Plotting RLC series circuit response using [Runge-Kutta methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) in MATLAB.

## Supported Runge-Kutta methods:
- Euler
- Midpoint
- Heun
- Ralston
- RK3
- RK4
- 3/8

## Usage

###`RLC(Vin_str,method)`

Parameters
- `Vin_str`: string type, Vin(t) as a function of `t`.
- `method`: string type, name of the RK method to use.

Examples:
- `RLC('sin(2*pi*100*t)','heun')` (input is a 100Hz sine wave, Heun's method)
- `RLC('Square(t,100)','rk4')` (input is a 100Hz square wave, RK4)

## Update
Added RL.m, a modified version of RLC.m for RL series circuit.