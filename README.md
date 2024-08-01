# Ion Acoustic Soliton Solver

This repository contains a Python script for solving the nonlinear normalized Poisson equation for electric potential in ion acoustic solitons using the Sagdeev method. The solver implements a set of equations to model the behavior of ion acoustic waves in plasma.

## Features

- Solves the nonlinear normalized Poisson equation for electric potential
- Implements the Sagdeev method for soliton analysis
- Calculates ion bulk velocity and normalized ion density
- Compares the solution with the Korteweg-de Vries (KdV) equation approximation
- Generates plots for visualization of results

## Dependencies

- NumPy
- SciPy
- Matplotlib

## Usage

1. Ensure all dependencies are installed.
2. Run the script:

```
python sagdeev_nonlinear_solve.py
```

3. The script will solve the equations and generate a plot showing:
   - The electric potential (Phi)
   - The KdV approximation of Phi
   - The ion bulk velocity
   - The normalized ion density

## Key Functions

- `sagdeev_ode_solver(phi0, v_d)`: Solves the ODE for the electric potential
- `solve_for_vd_func(phi, initial_vd_guess)`: Calculates the soliton Mach number
- `solve_for_u(phi, v_d)`: Computes the ion bulk velocity
- `solve_for_n(u, v_d)`: Calculates the normalized ion density

## Equations

The solver is based on the following set of equations:

1. Relation between drift velocity, ion velocity, and ion density:

v_d - u = v_d / n

2. Electric potential in terms of velocities:

φ = uv_d - u^2 / 2

3. Nonlinear Poisson equation:

∂^2φ / ∂ζ^2 = e^φ - n

4. Integrated form (in the absence of source term):

(1/2) * (∂φ/∂ζ)^2 - e^φ - v_d * √(v_d^2 - 2φ) + C = 0

with C = (1 + v_d^2)

These equations are transformed and solved to obtain the electric potential profile of the ion acoustic soliton.

## Notes

- The source term S(ζ) is set to 0 in this implementation.
- The initial conditions and parameters can be adjusted in the main script section.

## Contributing

Contributions to improve the solver or extend its capabilities are welcome. Please submit a pull request or open an issue to discuss proposed changes.


