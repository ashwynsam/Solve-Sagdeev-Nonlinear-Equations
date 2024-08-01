import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy import optimize

def sagdeev_ode_solver(phi0, v_d): #this function solves the nonlinear normalized Poisson equation for electric potential. This follows the Sagdeev method
    C = 1 + v_d**2

    def odefun(zeta, phi):
        return -np.sqrt(2 * (np.exp(phi) + v_d * np.sqrt(v_d**2 - 2*phi) - C))

    zeta_span = [0, 20]
    sol = solve_ivp(odefun, zeta_span, [phi0], dense_output=True, 
                    method='RK45', rtol=1e-6, atol=1e-7, max_step=0.1)

    zeta = sol.t
    phi = sol.y[0]

    zeta_mirrored = np.concatenate((-zeta[1:][::-1], zeta))
    phi_mirrored = np.concatenate((phi[1:][::-1], phi))

    return zeta_mirrored, phi_mirrored

def solve_for_vd_func(phi, initial_vd_guess):
    # Check if initial_guess^2 < 2*phi
    if initial_vd_guess**2 < 2*phi:
        print(f"Error: initial_guess^2 ({initial_vd_guess**2:.4f}) is less than 2*phi ({2*phi:.4f})")
        print("This violates the condition vd^2 >= 2*phi. Stopping the program.")
        sys.exit(1)  # Exit the program with an error code
    
    # Define the function
    def f(vd):
        return np.exp(phi) + vd * np.sqrt(vd**2 - 2*phi) - (1 + vd**2)

    # Use scipy.optimize.root_scalar to find the root 
    result = optimize.root_scalar(f, x0=initial_vd_guess, method='secant')
    vd = result.root

    # Display the result
    print(f'For phi = {phi:.4f}, v_d = {vd:.4f}')

    # Verify the solution
    lhs = np.exp(phi) + vd * np.sqrt(vd**2 - 2*phi)
    rhs = 1 + vd**2
    print(f'Verification: LHS = {lhs:.6f}, RHS = {rhs:.6f}')

    # Round vd to 4 decimal places and return. This rounding is to ensure a small perturbation 
    return round(vd, 4)

def solve_for_u(phi, v_d): #this function is the analytical solution to u if phi is given in the soliton frame. Note we are using the -sqrt() solution 
    phi_array = np.array(phi)
    
    # Check for valid input
    if np.any(v_d**2 < 2*phi_array):
        raise ValueError("Invalid input: v_d^2 must be greater than or equal to 2*phi for all phi values")
    
    u = v_d - np.sqrt(v_d**2 - 2*phi_array)
    return u

def solve_for_n(u, v_d): #this gives the normalized ion density with the ion bulk velocity provided in the soliton frame
    u = np.array(u) #doing this just in case u is manually put in as non numpy array
    n = v_d / (v_d - u)
    return n

# Main script
phi0 = 0.52
vd_guess = 1.2
v_d = solve_for_vd_func(phi0, vd_guess)
zeta, phi = sagdeev_ode_solver(phi0, v_d)
u = solve_for_u(phi, v_d)
n = solve_for_n(u, v_d)
phi_kdv = phi0 * np.power(np.cosh(np.sqrt(phi0/6) * zeta), -2)


# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(zeta, phi, 'b-', linewidth=2, label='Phi')
plt.plot(zeta, phi_kdv, ':', label='KdV Phi')
plt.plot(zeta, u, label= 'Velocity')
plt.plot(zeta, n - 1, label = 'Ion Density')
plt.xlabel(r'$\zeta$')
plt.ylabel(r'$\phi$')
plt.title('Solution of the ODE')
plt.grid(True)
plt.legend()
plt.show()

