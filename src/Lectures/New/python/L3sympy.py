#%%
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.matrices import Matrix 

# --- Gauss Quadrature Points and Weights ---
def get_gauss_points(n):
    """
    Returns Gauss points (xi) and weights (w) for n-point quadrature
    on the standard interval [-1, 1].
    """
    if n == 1:
        xi = [0.0]
        w = [2.0]
    elif n == 2:
        xi = [-1/np.sqrt(3), 1/np.sqrt(3)]
        w = [1.0, 1.0]
    elif n == 3:
        xi = [-np.sqrt(3/5), 0.0, np.sqrt(3/5)]
        w = [5/9, 8/9, 5/9]
    else:
        raise ValueError("Only 1, 2, or 3 Gauss points are implemented")
    return xi, w

def solve_1d_bar_fem(L, N, young_moduli, area, left_disp, right_traction, num_gauss_points=2,verifiaction=False):
    """
    Solves a 1D bar problem with the finite element method using shape functions
    on [-1,1] and Gauss quadrature for numerical integration.

    Parameters:
    -----------
    L : float
        Total length of the bar
    N : int
        Number of elements
    young_moduli : list
        List of Young's modulus values to be assigned to elements. If len(young_moduli) < N,
        the values will be cycled through the elements.
    area : float
        Cross-sectional area (assumed constant)
    left_disp : float
        Prescribed displacement at the left end
    right_traction : float
        Prescribed traction at the right end (Force/Area)
    num_gauss_points : int
        Number of Gauss points to use for integration (default=2)
    verifiaction : bool
        If True, verifies the stiffness matrix against analytical values (default=False)

    Returns:
    --------
    x_nodes : numpy array
        Nodal positions
    u : numpy array
        Nodal displacements
    element_strains : numpy array
        Strains evaluated at element center (xi=0)
    element_stresses : numpy array
        Stresses evaluated at element center (xi=0)
    element_centers : numpy array
        Positions of the center of each element
    reaction_force : float
        Reaction force at the left end
    element_E : numpy array
        Young's modulus assigned to each element
    """

    # Element length (constant for this discretization)
    h = L / N
    if h <= 0:
        raise ValueError("Length L must be positive and N must be >= 1.")

    # Node coordinates
    x_nodes = np.linspace(0, L, N+1)

    # Element centers for visualization
    element_centers = [(x_nodes[i] + x_nodes[i+1])/2 for i in range(N)]

    # Assign Young's modulus to each element
    element_E = np.zeros(N)
    if len(young_moduli) < N:
        for i in range(N):
            element_E[i] = young_moduli[i % len(young_moduli)]
    else:
        element_E[:N] = young_moduli[:N]

    # --- Symbolic Shape Functions & Derivatives (Master element [-1, 1]) ---
    xi = sp.Symbol('xi') # Local coordinate in [-1, 1]
    # Linear shape functions N(xi) 
    N1 = (1 - xi) / 2
    N2 = (1 + xi) / 2

    # Derivatives w.r.t. local coordinate xi 
    
    dN1_dxi = sp.diff(N1, xi) # -1/2
    dN2_dxi = sp.diff(N2, xi) # +1/2
    # Matrix of derivatives [dN1/dxi, dN2/dxi]
    dN_dxi = Matrix([[dN1_dxi, dN2_dxi]])

    # --- Get Gauss Points (standard interval [-1, 1]) ---
    gauss_xi, gauss_w = get_gauss_points(num_gauss_points)
    

    # --- Initialize Global Matrices ---
    
    K_global = Matrix.zeros(N+1)

    # --- Assemble Global Stiffness Matrix using Gauss Quadrature ---
    for e in range(N):
        # Get element properties
        E_e = element_E[e]

        # Element stiffness matrix initialization
        k_e = Matrix.zeros(2)

        # Jacobian of transformation: dx/dxi = h / 2 (Constant for linear element)
        J = h / 2.0
        if J <= 0: continue # Skip zero-length elements if they somehow occur

        # Loop over Gauss points
        for i in range(num_gauss_points):
            xi_i = gauss_xi[i]
            w_i = gauss_w[i]

            # Evaluate dN/dxi at Gauss point (constants for linear element)
            dN_dxi_i = dN_dxi.subs(xi, xi_i)

            # Calculate dN/dx = (dN/dxi) / (dx/dxi) = (dN/dxi) / J
            # Convert to float for B matrix calculation to avoid potential type issues later
            try:
                 dN_dx_i = (1.0 / J) * dN_dxi_i.evalf()
            except AttributeError: # in case dN_dxi_i is already float 
                 dN_dx_i = (1.0 / J) * dN_dxi_i


            # Form B matrix [dN1/dx, dN2/dx] at Gauss point
            B_i = dN_dx_i 
            if not isinstance(B_i, Matrix): B_i = Matrix(B_i)
            # Calculate contribution to stiffness matrix integral 

            k_contribution = B_i.transpose() * E_e * B_i * area * J * w_i
            k_e += k_contribution

        # --- Verification  ---
        if verifiaction:
            k00_numeric = float(k_e[0,0])
            expected_k_val = (E_e * area) / h
            # Use relative tolerance 
            if expected_k_val != 0:
                relative_diff = abs(k00_numeric - expected_k_val) / abs(expected_k_val)
                assert relative_diff < 1e-9, \
                    f"Element {e} stiffness error. Expected: {expected_k_val:.6g}, Calculated: {k00_numeric:.6g}"
            else: # 
                assert abs(k00_numeric) < 1e-10, \
                    f"Element {e} stiffness error. Expected: 0, Calculated: {k00_numeric:.6g}"


        # Get global indices
        idx1 = e
        idx2 = e + 1

        # Assemble into global stiffness matrix
        indices = [idx1, idx2]
        for r in range(2):
            for c in range(2):
                K_global[indices[r], indices[c]] += k_e[r, c]

    # --- Apply Loads and Boundary Conditions ---
    F_applied_external = Matrix.zeros(N+1, 1) 
    F_applied_external[N, 0] += right_traction * area

    # Store original stiffness matrix for reaction force calculation
    K_original = K_global.copy()

    # Apply fixed displacement at the left end 
    
    K_global[0, :] = Matrix.zeros(1, N+1) # Zero out row 0
    K_global[0, 0] = 1.0                  # Set diagonal to 1
    F_modified_for_solve = F_applied_external.copy() # Start with external forces
    F_modified_for_solve[0, 0] = left_disp         # Set BC value

    # --- Solve System ---
    try:
        # Use SymPy's LUsolve
        u_sympy = K_global.LUsolve(F_modified_for_solve)
        # Convert final displacements to NumPy array
        u = np.array(u_sympy).astype(float).flatten()

    except Exception as e:
        print(f"Error solving the system: {e}")
        return None, None, None, None, None, None, None # Match return tuple

    # --- Calculate Reaction Force (Fixed end)---
    
    u_col_sympy = Matrix(u) 
    internal_forces = K_original * u_col_sympy
    reaction_force = float(internal_forces[0,0] - F_applied_external[0,0])

    # --- Post-processing: Calculate Strains and Stresses ---
    element_strains = []
    element_stresses = []
    # Evaluate at element center (xi=0 for [-1, 1] interval)
    
    dN_dxi_center = Matrix([[sp.Rational(-1, 2), sp.Rational(1, 2)]])

    for e in range(N):
    
        J_e = h / 2.0
        if J_e <= 0:
            element_strains.append(np.nan)
            element_stresses.append(np.nan)
            continue

        dN_dx_center = (1.0 / J_e) * dN_dxi_center 

        # B matrix is [dN1/dx, dN2/dx]
        B_center = dN_dx_center

        # Get nodal displacements for this element
        u_e = Matrix([[u[e]], [u[e+1]]])

        # Calculate strain: epsilon = B * u_e (1x2 * 2x1 = 1x1)
        strain_e = (B_center * u_e)[0,0] 

        # Calculate stress: sigma = E * epsilon
        stress_e = element_E[e] * strain_e

        element_strains.append(float(strain_e))
        element_stresses.append(float(stress_e))

    return (x_nodes, u, np.array(element_strains), np.array(element_stresses),
            np.array(element_centers), reaction_force, element_E)


def plot_fem_results(x_nodes, u, element_strains, element_stresses, element_centers, 
                    young_moduli, reaction_force=None, scale_factor=10.0):
    """
    Plot the results of the FEM analysis for the 1D bar problem.
    
    Parameters:
    -----------
    x_nodes : numpy array
        Nodal positions
    u : numpy array
        Nodal displacements
    element_strains : numpy array
        Strains at the center of each element
    element_stresses : numpy array
        Stresses at the center of each element
    element_centers : numpy array
        Positions of the center of each element
    young_moduli : list or numpy array
        Young's modulus values used for each element
    reaction_force : float
        Reaction force at the left end (optional)
    scale_factor : float
        Factor to scale displacements for visibility (default=1.0)
    """
    
    if u is None:
        print("Cannot plot results as the solution failed.")
        return
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle('1D Bar Finite Element Analysis', fontsize=16)
    
    if reaction_force is not None:
        fig.text(0.5, 0.96, f'Reaction Force at Left End = {reaction_force:.4g} N', 
                 ha='center', va='center', fontsize=12, 
                 bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
    
    # 1. Displacement plot 
    ax1.plot(x_nodes, u, 'bo-', label='Nodal Displacement')
    ax1.set_ylabel('Displacement (u)')
    ax1.set_title('Displacement along the bar')
    ax1.grid(True)
    
    # For the deformed mesh
    ax1_twin = ax1.twinx()
    
    ax1_twin.plot(x_nodes, np.zeros_like(x_nodes), 'k--', label='Original Mesh')
    
    ax1_twin.plot(x_nodes + u * scale_factor, np.zeros_like(x_nodes), 'r-', 
                 label=f'Deformed Mesh (u×{scale_factor})')
    
    ax1_twin.set_ylabel('Mesh Visualization')
    ax1_twin.set_yticks([])
    

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
    
    # 2. Strain plot
    ax2.plot(element_centers, element_strains, 'ro-', label='Element Strain')
    ax2.set_ylabel('Strain (ε)')
    ax2.set_title('Strain distribution (element centers)')
    ax2.grid(True)
    ax2.legend()
    
    # 3. Stress plot 

    element_length = x_nodes[1] - x_nodes[0]
    ax3.bar(element_centers, element_stresses, width=element_length*0.8, 
            align='center', alpha=0.7, color='g', label='Element Stress')
    
    ax3.set_ylabel('Stress (σ)')
    ax3.set_xlabel('Position (x)')
    ax3.set_title('Stress distribution (element centers)')
    ax3.grid(True)
    ax3.legend(loc='upper left')
    
    # 4. Add Young's modulus distribution 
    ax3_twin = ax3.twinx()
    

    N = len(element_centers)
    L = x_nodes[-1]  
    
    x_boundaries = np.zeros(N+1)
    for i in range(N+1):
        x_boundaries[i] = i * element_length
        
    
    ax3_twin.step(x_boundaries, np.append(young_moduli, young_moduli[-1]), 'b-', where='post', 
                 label="Young's Modulus", alpha=0.5)
    ax3_twin.set_ylabel("Young's Modulus (E)")
    ax3_twin.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax3_twin.legend(loc='upper right')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def demonstrate_1d_bar_fem():
    """
    Demonstrate the solution of a 1D bar problem with FEM.
    
    It's intended as a quick example.
    """
    # Define problem parameters
    L = 2.0  # Length of the bar
    N = 20    # Number of elements
    
    # Alternating Young's modulus
    young_moduli =[1e7]*5+ 5*[1e6]+10*[1e7]
    
    area = 0.01          # Cross-sectional area (m²)
    left_disp = 0.0      # Fixed left end (m)
    right_traction = 1000 # Traction at right end (Pa)
    
    # Number of Gauss points to use
    num_gauss_points = 2  # 2-point Gauss quadrature
    
    
    # Solve the problem
    x_nodes, u, element_strains, element_stresses, element_centers, reaction_force, element_E = solve_1d_bar_fem(
        L, N, young_moduli, area, left_disp, right_traction, num_gauss_points
    )
    
    # Plot the results
    if u is not None:
        print("FEM Solution Completed Successfully!")
        print(f"Number of Elements: {N}")
        print(f"Number of Nodes: {N+1}")
        print(f"Number of Gauss Points: {num_gauss_points}")
        print("\nNodal Coordinates:")
        print(x_nodes)
        print("\nNodal Displacements:")
        print(u)
        print("\nElement Strains:")
        print(element_strains)
        print("\nElement Stresses:")
        print(element_stresses)
        print("\nElement Young's Moduli:")
        print(element_E)
        print(f"\nReaction Force at Left End: {reaction_force:.4g} N")
        
        # Compare with analytical solution for uniform bar
        if len(set(young_moduli)) == 1:  # Check if all moduli are the same
            E = young_moduli[0]
            analytical_strain = right_traction / E
            analytical_disp_end = right_traction * L / E
            print("\nAnalytical Solution (uniform bar):")
            print(f"Strain: {analytical_strain}")
            print(f"Displacement at x=L: {analytical_disp_end}")
            print(f"Numerical Error: {abs(u[-1] - analytical_disp_end)/analytical_disp_end:.6%}")
        
        plot_fem_results(x_nodes, u, element_strains, element_stresses, 
                        element_centers, element_E, reaction_force, scale_factor=100)


if __name__ == "__main__":
    demonstrate_1d_bar_fem()
    
#%%