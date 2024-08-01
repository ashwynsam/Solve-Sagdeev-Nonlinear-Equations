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

   ![equation](https://latex.codecogs.com/gif.latex?v_d%20-%20u%20%3D%20%5Cfrac%7Bv_d%7D%7Bn%7D)

2. Electric potential in terms of velocities:

   ![equation](https://latex.codecogs.com/gif.latex?%5Cphi%20%3D%20uv_d%20-%20%5Cfrac%7Bu%5E2%7D%7B2%7D)

3. Nonlinear Poisson equation:

   ![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%5E2%5Cphi%7D%7B%5Cpartial%5Czeta%5E2%7D%20%3D%20e%5E%5Cphi%20-%20n)

4. Integrated form (in the absence of source term):

   ![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%5Cfrac%7B%5Cpartial%5Cphi%7D%7B%5Cpartial%5Czeta%7D%5Cright%29%5E2%20-%20e%5E%5Cphi%20-%20v_d%5Csqrt%7Bv_d%5E2%20-%202%5Cphi%7D%20&plus;%20C%20%3D%200)

These equations are transformed and solved to obtain the electric potential profile of the ion acoustic soliton.

## Notes

- The source term S(ζ) is set to 0 in this implementation.
- The initial conditions and parameters can be adjusted in the main script section.

## Contributing

Contributions to improve the solver or extend its capabilities are welcome. Please submit a pull request or open an issue to discuss proposed changes.


