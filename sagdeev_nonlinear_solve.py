import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy import optimize

def sagdeev_ode_solver(phi0, v_d): #this function solves the nonlinear normalized Poisson equation for electric potential. This follows the Sagdeev method
    C = 1 + v_d**2
    if (v_d**2 - 2*phi0) < 0:
        raise ValueError("Error: Taking the squqare root of a negative number. Check if v_d was calculated outside of its range")

    def odefun(zeta, phi):
        arg = 2 * (np.exp(phi) + v_d * np.sqrt(v_d**2 - 2*phi) - C)
        if arg <= 0:
            arg = 1e-5 # This epsilon is to ensure a small perturbation for where the derivative equals 0 or if argument is negative due to rounding errors

        return -np.sqrt(arg)

    zeta_span = [0, 20]
    sol = solve_ivp(odefun, zeta_span, [phi0], dense_output=True, 
                    method='RK45', rtol=1e-6, atol=1e-7, max_step=0.1)

    zeta = sol.t
    phi = sol.y[0]

    zeta_mirrored = np.concatenate((-zeta[1:][::-1], zeta))
    phi_mirrored = np.concatenate((phi[1:][::-1], phi))

    return zeta_mirrored, phi_mirrored

def solve_for_vd_func(phi):
    v_d = np.sqrt((np.exp(phi) - 1)**2 / (2 * (np.exp(phi) - 1 - phi)))
    return v_d 

def solve_for_u(phi, v_d): #this function is the analytical solution to u if phi is given in the soliton frame. Note we are using the -sqrt() solution 
    phi_array = np.array(phi)
    
    # Check for valid input
    if np.any(v_d**2 < 2*phi_array):
        raise ValueError("Error: Taking the squqare root of a negative number. Check if v_d was calculated outside of its range")
    
    u = v_d - np.sqrt(v_d**2 - 2*phi_array)
    return u

def solve_for_n(u, v_d): #this gives the normalized ion density with the ion bulk velocity provided in the soliton frame
    u = np.array(u) #doing this just in case u is manually put in as non numpy array
    n = v_d / (v_d - u)
    return n

# Main script
phi0 = 0.9
v_d = solve_for_vd_func(phi0)
print(v_d)
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
plt.title('Solution of the full set of nonlinear fluid equations')
plt.grid(True)
plt.legend()
plt.show()

