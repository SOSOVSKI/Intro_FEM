# %%
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.special import legendre
#%%
import sympy as sp
from sympy import symbols, Function, Derivative, Integral, diff, integrate, Matrix, Symbol, Idx, IndexedBase, Sum, Eq, Q, Wild, MatrixSymbol
from sympy.printing import pprint, latex
from IPython.display import Math, display, HTML, Markdown
import inspect 

# Use HTML for styled print
def print_styled(text, size="1.1em",htm=False):
    if htm:
        display(HTML(f"<div style='margin: 5px 0; font-size: {size};'>{text}</div>"))
    else:
        # display(Markdown(f"<div style='margin: 5px 0; font-size: {size};'>{text}</div>"))
        display(Markdown(text))

def Lprint(expr):
    """
    Displays SymPy expressions or raw LaTeX strings using MathJax.
    """
    if not isinstance(expr, str) or hasattr(expr, '_repr_latex_'):
        latex_str = sp.latex(expr, mode='plain')
    else:
        latex_str = expr
    if not latex_str.startswith('$') and not latex_str.startswith('\\('):
         latex_str = f"${latex_str}$"
    display(Math(latex_str))

# ---- FACTORY FUNCTION FOR CUSTOM INDEXED FUNCTIONS ----
def create_latex_indexed_function(symbolic_name: str, latex_symbol: str):
    """
    Factory to create a SymPy Function subclass that renders as symbol_{index}(var).

    """
    class LatexIndexedFunc(Function):
        
        _latex_symbol_str = latex_symbol
        is_commutative = False # Functions are generally not commutative

        def _latex(self, printer):

            if len(self.args) == 2:

                func_name = self._latex_symbol_str
                index_latex = printer._print(self.args[0])
                var_latex = printer._print(self.args[1])
                return rf"{func_name}_{{{index_latex}}}({var_latex})"
            else:
                cls_name = self.__class__.__name__
                args_latex = ", ".join(printer._print(arg) for arg in self.args)
                return rf"{cls_name}({args_latex})"

    LatexIndexedFunc.__name__ = symbolic_name
    return LatexIndexedFunc

def create_fem_diagram(
    width=10,                  # Domain width
    height=1.0,                # Domain height
    num_elements=5,            # Number of elements
    element_color='#a0a0a0',   # Base color if not using gradient
    domain_color='none',       # Face color of the domain rectangle ('none' for transparent)
    use_gradient=True,         # Use gradient coloring for elements
    gradient_cmap='Blues',     # Colormap for gradient
    show_element_labels=True,  # Control visibility of E_i labels and arrows
    save_path=None,            # Path to save figure, None for display only
    dpi=300                    # Resolution for saved figure
):
    """
    Creates an improved Finite Element Method (FEM) 1D domain diagram
    with element labels and arrows.

    Args:
        width (float): The total width of the domain (Omega).
        height (float): The height of the rectangular representation.
        num_elements (int): The number of finite elements to display.
        element_color (str): The color for elements if use_gradient is False.
        domain_color (str): The background color of the main domain rectangle.
                            Use 'none' for no fill.
        use_gradient (bool): Whether to apply a color gradient across elements.
        gradient_cmap (str): The name of the matplotlib colormap for the gradient.
        show_element_labels (bool): If True, display E_i labels and arrows.
        save_path (str, optional): File path to save the figure. If None, displays the plot.
        dpi (int): Dots per inch for the saved figure.
    """
    if num_elements <= 0:
        print("Number of elements must be positive.")
        return

    # --- Figure Setup ---
    fig_width = 8
    # Adjust fig_height calculation slightly to give more vertical space
    fig_height = height * (fig_width / width) * 3.5
    plt.figure(figsize=(fig_width, fig_height), dpi=150)
    ax = plt.gca()

    # --- Domain Geometry ---
    x_min, x_max = 0, width
    y_min, y_max = 0, height
    element_width = width / num_elements
    node_positions = np.linspace(x_min, x_max, num_elements + 1)

    # --- Draw Domain Outline ---
    domain_rect = patches.Rectangle(
        (x_min, y_min),
        width,
        height,
        facecolor=domain_color,
        edgecolor='black',
        linewidth=2,
        zorder=1 # Ensure it's behind elements if filled
    )
    ax.add_patch(domain_rect)

    # --- Draw Elements, Labels, and Arrows ---
    if use_gradient:
        cmap = cm.get_cmap(gradient_cmap)
        colors = cmap(np.linspace(0.3, 0.9, num_elements))
    else:
        colors = []
        color_1 = '#81a2be'  # Blue for elements that are odd or multiples of 3
        color_2 = '#cc6666'  # Red for elements that are even or multiples of 5
        
        for i in range(num_elements):
            element_num = i + 1
            if element_num % 2 == 1 or element_num % 3 == 0:
                colors.append(color_1)
            elif element_num % 2 == 0 or element_num % 5 == 0:
                colors.append(color_2)

    # Define arrow properties
    arrow_props = dict(arrowstyle='->', linewidth=1.5, color='black')

    for i in range(num_elements):
        x_left = node_positions[i]
        element = patches.Rectangle(
            (x_left, y_min),
            element_width,
            height,
            facecolor=colors[i],
            edgecolor='#303030', # Darker edge for elements
            linewidth=1,
            zorder=2 # Elements on top of domain rect
        )
        ax.add_patch(element)

        if show_element_labels:
            # --- Element Labels (Inside, near top) ---
            label_x = x_left + element_width / 2
            label_y = y_max - 0.2 * height # Position inside near top
            # Determine which label to show based on the rules:
            # E_1 if i+1 is odd or a multiple of 3
            # E_2 if i+1 is even or a multiple of 5
            element_num = i+1
            if element_num % 2 == 1 or element_num % 3 == 0:
                label = 'E_1'
            elif element_num % 2 == 0 or element_num % 5 == 0:
                label = 'E_2'
            else:
                label = f'E_{{{element_num}}}'  # Fallback (shouldn't be reached with these rules)
                
            plt.text(label_x, label_y, rf'${label}$',
                     fontsize=11, ha='center', va='center', zorder=4) # zorder 4 for label

            # # --- Arrow from Label to Element Center ---
            # arrow_start_y = label_y - 0.15 * height # Start slightly below label text
            # arrow_end_y = y_min + 0.3 * height    # Point towards lower part of element
            # plt.annotate('', xy=(label_x, arrow_end_y), xytext=(label_x, arrow_start_y),
            #              arrowprops=arrow_props, zorder=3) # zorder 3 for arrow

    # --- Draw Nodes ---
    node_y = y_min + height / 2 # Center vertically
    ax.plot(node_positions, np.full(num_elements + 1, node_y), 'ko', # Black circles
            markersize=4, zorder=3) # Nodes on top

    # --- Boundary Markers (Diagonal Lines) ---
    num_markers = 4
    marker_length = 0.15 * height
    # Left boundary (Gamma_u)
    for y_pos in np.linspace(y_min + marker_length, y_max - marker_length, num_markers):
        ax.plot([x_min - marker_length, x_min], [y_pos + marker_length, y_pos],
                'k-', linewidth=1.5, zorder=3)
    # --- Boundary Markers (traction arrow) ---
    # Right boundary (Gamma_t)
    # Draw right-facing arrow from middle of right edge
    arrow_length = 0.3 * height
    arrow_y = y_min + height/2  # Middle of right edge
    ax.arrow(x_max, arrow_y, 
             arrow_length, 0,  # x direction only, no y component
             head_width=0.2*height, 
             head_length=0.1*height,
             fc='black', ec='black',
             linewidth=1.5, zorder=3)

    # --- Text Labels ---
    # Y positions adjusted slightly if needed
    omega_y_offset = y_max + 1.0 * height
    bc_text_y_offset = y_max + 0.5 * height
    gamma_label_y_offset = y_min - 0.8 * height # Keep Gamma labels below

    # Domain Label (Omega)
    plt.text(width / 2, omega_y_offset, r'$\Omega$',
             fontsize=16, ha='center', va='center', zorder=3)

    # Boundary Condition Labels (Above)
    gamma_x_offset = marker_length + 0.1 * width
    plt.text(x_min -gamma_x_offset, bc_text_y_offset, r'$u=u^*$ specified',
             fontsize=11, ha='left', va='center', zorder=3)
    plt.text(x_max +gamma_x_offset, bc_text_y_offset, r'$t=t^*$ specified',
             fontsize=11, ha='right', va='center', zorder=3)

    # Gamma Labels (Below Boundaries)
    gamma_x_offset = marker_length + 0.1 * width
    plt.text(x_min - gamma_x_offset, gamma_label_y_offset, r'$\Gamma_u$',
             fontsize=14, ha='center', va='center', zorder=3)
    plt.text(x_max + gamma_x_offset, gamma_label_y_offset, r'$\Gamma_t$',
             fontsize=14, ha='center', va='center', zorder=3)
    # --- Arrow from Gamma_u label to left boundary ---
    plt.annotate('', 
                 xy=(x_min, y_min + 0.5 * height), # arrow end: left boundary, middle height
                 xytext=(x_min - gamma_x_offset, gamma_label_y_offset), # start from Gamma_u label
                 arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'),
                 zorder=3)
    # arrow from Gamma_t label to right boundary
    plt.annotate('', 
                 xy=(x_max, y_min + 0.5 * height), # arrow end: right boundary, middle height
                 xytext=(x_max + gamma_x_offset, gamma_label_y_offset), # start from Gamma_t label
                 arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'),
                 zorder=3)
    # --- Domain Width Arrow ---
    arrow_y = y_max + 0.75 * height # Positioned between Omega and BC text
    arrow_props_wide = dict(arrowstyle='<->', linewidth=1.5, color='black')
    plt.annotate('', xy=(x_min, arrow_y), xytext=(x_max, arrow_y),
                 arrowprops=arrow_props_wide, zorder=3)

    # --- Final Plot Adjustments ---
    x_pad = gamma_x_offset + 0.3 # Adjust padding based on Gamma label position
    y_pad_top = omega_y_offset - y_max + 0.2 * height # Space needed above y_max
    y_pad_bottom = abs(gamma_label_y_offset) + 0.2 * height # Space needed below y_min

    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad_bottom, y_max + y_pad_top) # Use calculated padding

    ax.set_aspect('equal', adjustable='box') # Maintain aspect ratio
    plt.axis('off') # Hide axes
    plt.tight_layout(pad=0.5) # Adjust layout to prevent overlap

    # --- Save or Show ---
    if save_path:
        fig_facecolor = 'white' if domain_color == 'none' else 'auto'
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor=fig_facecolor, transparent=False)
        print(f"Figure saved to {save_path}")
    else:
        plt.show()

#%%


class FEMShapeFunctions:
    def __init__(self):
        """Initialize the FEM shape functions visualization tool."""
        pass
    
    def lagrange_basis(self, xi, degree, node_index):
        """
        Compute Lagrange basis function of specified degree at point xi.
        
        Parameters:
        -----------
        xi : float or numpy.ndarray
            Point(s) at which to evaluate the basis function (-1 <= xi <= 1)
        degree : int
            Degree of the polynomial (0 = constant, 1 = linear, 2 = quadratic, etc.)
        node_index : int
            Index of the node associated with this basis function (0 to degree)
            
        Returns:
        --------
        float or numpy.ndarray
            Value of the basis function at xi
        """
        # Create the node locations (equidistant in the reference element)
        nodes = np.linspace(-1, 1, degree + 1)
        
        # Initialize the basis function value to 1
        phi = np.ones_like(xi, dtype=float)
        
        # Lagrange polynomial construction
        for j in range(degree + 1):
            if j != node_index:
                phi *= (xi - nodes[j]) / (nodes[node_index] - nodes[j])
                
        return phi
    
    def hierarchical_basis(self, xi, degree, node_index):
        """
        Compute hierarchical basis function of specified degree at point xi.
        
        Parameters:
        -----------
        xi : float or numpy.ndarray
            Point(s) at which to evaluate the basis function (-1 <= xi <= 1)
        degree : int
            Maximum degree of the polynomial
        node_index : int
            Index of the basis function (0 to degree)
            
        Returns:
        --------
        float or numpy.ndarray
            Value of the basis function at xi
        """
        if node_index == 0:
            # Left node basis function
            return (1 - xi) / 2
        elif node_index == 1 and degree >= 1:
            # Right node basis function
            return (1 + xi) / 2
        else:
            # Interior (bubble) functions for higher-order elements
            n = node_index - 1  # Index of the interior mode
            # Calculate the Legendre polynomial of degree n
            P_n = legendre(n)
            # Integrate to get the integrated Legendre polynomial
            return np.sqrt(2 * n + 1) * (1 - xi**2) * P_n(xi)
    
    def plot_multielement_shape_functions(self, num_elements=3, element_type='lagrange', degrees=[1, 2, 3], figsize=(15, 15)):
        """
        Plot shape functions across multiple elements.
        
        Parameters:
        -----------
        num_elements : int
            Number of elements to include
        element_type : str
            Type of basis functions ('lagrange' or 'hierarchical')
        degrees : list
            List of polynomial degrees to plot
        figsize : tuple
            Figure size
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        fig = plt.figure(figsize=figsize)
        
        # Create a grid layout
        n_rows = len(degrees)
        gs = GridSpec(n_rows, 1)
        
        # Function to use based on element type
        if element_type.lower() == 'lagrange':
            basis_function = self.lagrange_basis
            title_prefix = "Lagrange"
        elif element_type.lower() == 'hierarchical':
            basis_function = self.hierarchical_basis
            title_prefix = "Hierarchical"
        else:
            raise ValueError(f"Unknown element type: {element_type}")
        
        for i, degree in enumerate(degrees):
            ax = fig.add_subplot(gs[i, 0])
            
            # Create nodes for the global mesh
            global_nodes = np.linspace(0, 1, num_elements*degree + 1)
            
            # Create a fine evaluation grid
            x = np.linspace(0, 1, 500)
            
            # Dictionary to store all shape function values for later sum calculation
            all_phi_values = {}
            
            # Plot each global shape function
            for node_idx in range(len(global_nodes)):
                # Initialize shape function values
                phi_values = np.zeros_like(x)
                
                # Determine which elements this node contributes to
                for e in range(num_elements):
                    element_start = e
                    element_end = e + 1
                    
                    # Local nodes in this element
                    local_nodes = np.linspace(element_start/num_elements, 
                                             element_end/num_elements, 
                                             degree + 1)
                    
                    # If this global node is in this element
                    if node_idx in range(e*degree, (e+1)*degree + 1):
                        # Map from global to local node index
                        local_index = node_idx - e*degree
                        
                        # Points in this element
                        element_mask = (x >= element_start/num_elements) & (x <= element_end/num_elements)
                        
                        if any(element_mask):
                            # Map x to the reference element [-1, 1]
                            xi = 2 * (x[element_mask] - element_start/num_elements) / (1/num_elements) - 1
                            
                            # Compute shape function values
                            phi_values[element_mask] = basis_function(xi, degree, local_index)
                
                # Plot the shape function
                ax.plot(x, phi_values, label=f'$\\phi_{{{node_idx}}}$')
                
                # Store for sum calculation
                all_phi_values[node_idx] = phi_values
            
            # Calculate sum of all shape functions
            sum_phi = np.zeros_like(x)
            for node_idx in range(len(global_nodes)):
                sum_phi += all_phi_values[node_idx]
            
            # Plot the sum
            ax.plot(x, sum_phi, 'k--', lw=3, label='$\\sum \\phi_i = 1$')
            
            # Plot element boundaries
            for e in range(num_elements+1):
                ax.axvline(x=e/num_elements, color='r', linestyle='--', alpha=0.5)
            
            # Plot node positions
            ax.scatter(global_nodes, np.zeros_like(global_nodes), c='red', s=50, zorder=10)
            
            # Set axis limits and title
            ax.set_ylim(-0.1, 1.3)
            ax.set_title(f"{title_prefix} Shape Functions - {num_elements} Elements (Degree {degree})")
            
            # Add legend and axis labels
            if len(global_nodes) > 8:
                # For many global nodes, just show a simplified legend
                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                         ncol=min(8, len(global_nodes)+1),
                         labels=['Shape Functions', 'Sum = 1'])
            else:
                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                         ncol=min(8, len(global_nodes)+1))
            
            ax.set_xlabel('x (Physical Coordinate)')
            ax.set_ylabel('$\\phi(x)$')
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
        plt.tight_layout()
        return fig
    
    def plot_multielement_shape_function_derivatives(self, num_elements=3, element_type='lagrange', degrees=[1, 2], figsize=(15, 10)):
        """
        Plot derivatives of shape functions across multiple elements.
        
        Parameters:
        -----------
        num_elements : int
            Number of elements to include
        element_type : str
            Type of basis functions ('lagrange' or 'hierarchical')
        degrees : list
            List of polynomial degrees to plot
        figsize : tuple
            Figure size
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        fig = plt.figure(figsize=figsize)
        
        # Create a grid layout
        n_rows = len(degrees)
        gs = GridSpec(n_rows, 1)
        
        # Function to use based on element type
        if element_type.lower() == 'lagrange':
            basis_function = self.lagrange_basis
            title_prefix = "Lagrange"
        elif element_type.lower() == 'hierarchical':
            basis_function = self.hierarchical_basis
            title_prefix = "Hierarchical"
        else:
            raise ValueError(f"Unknown element type: {element_type}")
        
        for i, degree in enumerate(degrees):
            ax = fig.add_subplot(gs[i, 0])
            
            # Create nodes for the global mesh
            global_nodes = np.linspace(0, 1, num_elements*degree + 1)
            
            # Create a fine evaluation grid
            x = np.linspace(0, 1, 500)
            h = 1e-6  # Step size for numerical differentiation
            
            # Plot each global shape function derivative
            for node_idx in range(len(global_nodes)):
                # Initialize derivative values
                dphi_dx_values = np.zeros_like(x)
                
                # Determine which elements this node contributes to
                for e in range(num_elements):
                    element_start = e
                    element_end = e + 1
                    
                    # If this global node is in this element
                    if node_idx in range(e*degree, (e+1)*degree + 1):
                        # Map from global to local node index
                        local_index = node_idx - e*degree
                        
                        # Points in this element
                        element_mask = (x >= element_start/num_elements) & (x <= element_end/num_elements)
                        
                        if any(element_mask):
                            # Map x to the reference element [-1, 1]
                            xi = 2 * (x[element_mask] - element_start/num_elements) / (1/num_elements) - 1
                            
                            # Calculate the Jacobian
                            J = (1/num_elements) / 2  # dx/dxi
                            
                            # Compute shape function derivatives using central finite difference
                            phi_plus = basis_function(xi + h, degree, local_index)
                            phi_minus = basis_function(xi - h, degree, local_index)
                            dphi_dxi = (phi_plus - phi_minus) / (2 * h)
                            
                            # Convert to physical derivative
                            dphi_dx_values[element_mask] = dphi_dxi / (2 * J)
                
                # Plot the shape function derivative
                ax.plot(x, dphi_dx_values, label=f'$\\frac{{d\\phi_{{{node_idx}}}}}{{dx}}$')
            
            # Plot element boundaries
            for e in range(num_elements+1):
                ax.axvline(x=e/num_elements, color='r', linestyle='--', alpha=0.5)
            
            # Plot horizontal line at y=0
            ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            
            # Set title
            ax.set_title(f"{title_prefix} Shape Function Derivatives - {num_elements} Elements (Degree {degree})")
            
            # Add legend and axis labels
            if len(global_nodes) > 8:
                # For many global nodes, just show a simplified legend
                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                         ncol=min(8, 4),
                         labels=['Shape Function Derivatives'])
            else:
                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                         ncol=min(8, len(global_nodes)))
            
            ax.set_xlabel('x (Physical Coordinate)')
            ax.set_ylabel('$\\frac{d\\phi}{dx}$')
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
        plt.tight_layout()
        return fig
    
    def plot_multielement_partition_of_unity_per_element(self, num_elements=3, degree=1, figsize=(15, 8)):
        """
        Plot the sum of shape functions within each element to illustrate local partition of unity.
        
        Parameters:
        -----------
        num_elements : int
            Number of elements to include
        degree : int
            Polynomial degree of the shape functions
        figsize : tuple
            Figure size
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
        
        # Create nodes for the global mesh
        global_nodes = np.linspace(0, 1, num_elements*degree + 1)
        
        # Create a fine evaluation grid
        x = np.linspace(0, 1, 500)
        
        # Dictionary to store all shape function values
        all_phi_values = {}
        
        # Plot each global shape function
        for node_idx in range(len(global_nodes)):
            # Initialize shape function values
            phi_values = np.zeros_like(x)
            
            # Determine which elements this node contributes to
            for e in range(num_elements):
                element_start = e
                element_end = e + 1
                
                # If this global node is in this element
                if node_idx in range(e*degree, (e+1)*degree + 1):
                    # Map from global to local node index
                    local_index = node_idx - e*degree
                    
                    # Points in this element
                    element_mask = (x >= element_start/num_elements) & (x <= element_end/num_elements)
                    
                    if any(element_mask):
                        # Map x to the reference element [-1, 1]
                        xi = 2 * (x[element_mask] - element_start/num_elements) / (1/num_elements) - 1
                        
                        # Compute shape function values
                        phi_values[element_mask] = self.lagrange_basis(xi, degree, local_index)
            
            # Plot the shape function
            ax1.plot(x, phi_values, label=f'$\\phi_{{{node_idx}}}$')
            
            # Store for sum calculation
            all_phi_values[node_idx] = phi_values
        
        # Calculate sum of all shape functions
        sum_phi = np.zeros_like(x)
        for node_idx in range(len(global_nodes)):
            sum_phi += all_phi_values[node_idx]
        
        # Plot the sum on the first subplot
        ax1.plot(x, sum_phi, 'k--', lw=3, label='$\\sum \\phi_i = 1$')
        
        # Plot element boundaries
        for e in range(num_elements+1):
            ax1.axvline(x=e/num_elements, color='r', linestyle='--', alpha=0.5)
            ax2.axvline(x=e/num_elements, color='r', linestyle='--', alpha=0.5)
        
        # Add labels for the first subplot
        ax1.set_title(f"Global Shape Functions - {num_elements} Elements (Degree {degree})")
        ax1.set_xlabel('x (Physical Coordinate)')
        ax1.set_ylabel('$\\phi(x)$')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), 
                   ncol=min(8, len(global_nodes)+1))
        
        # Now create a plot showing the contribution of shape functions within each element
        colors = plt.cm.tab20(np.linspace(0, 1, num_elements))
        
        for e in range(num_elements):
            element_start = e / num_elements
            element_end = (e + 1) / num_elements
            
            # Points in this element
            element_mask = (x >= element_start) & (x <= element_end)
            x_element = x[element_mask]
            
            # Calculate sum of shape functions active in this element
            element_sum = np.zeros_like(x_element)
            
            # Plot the local shape functions for this element
            for local_idx in range(degree + 1):
                global_idx = e * degree + local_idx
                
                # Get the values of this shape function in the element
                phi_element = all_phi_values[global_idx][element_mask]
                
                # Add to the element sum
                element_sum += phi_element
                
                # Only plot a subset of local shape functions to avoid clutter
                if e == 0 or e == 1 or e == num_elements - 1:
                    ax2.plot(x_element, phi_element, 
                             label=f'$\\phi_{{{global_idx}}}$ (Element {e+1})',
                             color=colors[e], linestyle=['solid', 'dashed', 'dotted', 'dashdot'][local_idx % 4],
                             alpha=0.7)
            
            # Plot the element sum
            ax2.plot(x_element, element_sum, 'k--', lw=2)
            
            # Annotate the element
            mid_x = (element_start + element_end) / 2
            ax2.annotate(f"Element {e+1}", 
                        xy=(mid_x, 0.5), 
                        xytext=(mid_x, 0.5),
                        ha='center', va='center',
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
        
        # Set title and labels for the second subplot
        ax2.set_title("Partition of Unity Within Each Element")
        ax2.set_xlabel('x (Physical Coordinate)')
        ax2.set_ylabel('$\\phi(x)$')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(-0.1, 1.3)
        
        # Create a custom legend for the second subplot
        handles, labels = ax2.get_legend_handles_labels()
        if handles:
            ax2.legend(handles[:6], labels[:6], loc='upper center', 
                       bbox_to_anchor=(0.5, -0.15), ncol=3)
        
        plt.tight_layout()
        return fig


#%%
x, L = symbols('x L', real=True, positive=True)
i = Symbol('i', integer=True, positive=True)
j = Symbol('j', integer=True, positive=True)
N = Symbol('N', integer=True, positive=True)
y = Function('y')(x)
v = Function('v')(x)
A, B, C = symbols('A B C', real=True)
F = Function('F')(x) 
a = IndexedBase('a')
b = IndexedBase('b')
sigma = Function('sigma')(x)
u = Function('u')(x)
f = Symbol('f', real=True)
r = Symbol('r', real=True)
epsilon = Function('epsilon')(x)
E = Symbol('E', real=True, positive=True)
u_star = symbols('u^*')
t_star = symbols('t^*')
Gamma_u = symbols('\\Gamma_u')
Gamma_t = symbols('\\Gamma_t')
t_star = symbols('t^*')