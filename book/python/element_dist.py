#%%
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def get_default_c3d20_nodes():
    """
    Generates the coordinates for a C3D20 element
    (-1 to 1 in each direction).
    Node ordering follows  VTK and Abaqus conventions (e.g., Abaqus).
    Nodes 1-8 are corners.
    Nodes 9-20 are midside nodes.
    """
    return( np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], 
                          [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
                          [0, -1, -1], [1, 0, -1], [0, 1, -1], [-1, 0, -1], 
                          [0, -1, 1], [1, 0, 1], [0, 1, 1], [-1, 0, 1],
                          [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]],dtype=np.float64))
def get_SF(x,y,z):
    xm = 1.0 - x
    ym = 1.0 - y
    zm = 1.0 - z
    xp = 1.0 + x
    yp = 1.0 + y
    zp = 1.0 + z
    omx2 = 1.0 - x*x
    omy2 = 1.0 - y*y
    omz2= 1.0 - z*z
    return(np.asarray([xm * ym * zm * (-2.0 - x - y - z)/8.0,
                                xp * ym * zm * (-2.0 + x - y - z)/8.0,
                                xp * yp * zm * (-2.0 + x + y - z)/8.0,
                                xm * yp * zm * (-2.0 - x + y - z)/8.0,
                                xm * ym * zp * (-2.0 - x - y + z)/8.0,
                                xp * ym * zp * (-2.0 + x - y + z)/8.0,
                                xp * yp * zp * (-2.0 + x + y + z)/8.0,
                                xm * yp * zp * (-2.0 - x + y + z)/8.0,
                                omx2 * ym * zm/4.0,
                                xp * omy2 * zm/4.0,
                                omx2 * yp * zm/4.0,
                                xm * omy2 * zm/4.0,
                                omx2 * ym * zp/4.0,
                                xp * omy2 * zp/4.0,
                                omx2 * yp * zp/4.0,
                                xm * omy2 * zp/4.0,
                                xm * ym * omz2/4.0,
                                xp * ym * omz2/4.0,
                                xp * yp * omz2/4.0,
                                xm * yp * omz2/4.0],dtype=np.float64))
def get_dN(x,y,z):
        xm = 1.0 - x
        ym = 1.0 - y
        zm = 1.0 - z
        xp = 1.0 + x
        yp = 1.0 + y
        zp = 1.0 + z
        omx2 = 1.0 - x*x
        omy2 = 1.0 - y*y
        omz2= 1.0 - z*z
        return np.asarray([[-xm * ym * zm / 8.0 - ym * (-2.0 - x - y - z) * zm / 8.0,
                         -xm * ym * zm / 8.0 - xm * (-2.0 - x - y - z) * zm / 8.0,
                          -xm * ym * (-2.0 - x - y - z) / 8.0 - xm * ym * zm / 8.0],
                          
                          [xp * ym * zm / 8.0 + ym * (-2.0 + x - y - z) * zm / 8.0,
                            -xp * ym * zm / 8.0 - xp * (-2.0 + x - y - z) * zm / 8.0,
                           -xp * ym * (-2.0 + x - y - z) / 8.0 - xp * ym * zm / 8.0 ],

                           [xp * yp * zm / 8.0 + yp * (-2.0 + x + y - z) * zm / 8.0,
                            xp * yp * zm / 8.0 + xp * (-2.0 + x + y - z) * zm / 8.0,
                            -xp * yp * (-2.0 + x + y - z) / 8.0 - xp * yp * zm / 8.0],

                           [-xm * yp * zm / 8.0 - yp * (-2.0 - x + y - z) * zm / 8.0,
                            xm * yp * zm / 8.0 + xm * (-2.0 - x + y - z) * zm / 8.0,
                            -xm * yp * (-2.0 - x + y - z) / 8.0 - xm * yp * zm / 8.0],

                           [-xm * ym * zp / 8.0 - ym * (-2.0 - x - y + z) * zp / 8.0,
                            -xm * ym * zp / 8.0 - xm * (-2.0 - x - y + z) * zp / 8.0,
                            xm * ym * (-2.0 - x - y + z) / 8.0 + xm * ym * zp / 8.0],

                           [xp * ym * zp / 8.0 + ym * (-2.0 + x - y + z) * zp / 8.0,
                            -xp * ym * zp / 8.0 - xp * (-2.0 + x - y + z) * zp / 8.0,
                            xp * ym * (-2.0 + x - y + z) / 8.0 + xp * ym * zp / 8.0],

                           [xp * yp * zp / 8.0 + yp * (-2.0 + x + y + z) * zp / 8.0,
                            xp * yp * zp / 8.0 + xp * (-2.0 + x + y + z) * zp / 8.0,
                            xp * yp * (-2.0 + x + y + z) / 8.0 + xp * yp * zp / 8.0],

                           [-xm * yp * zp / 8.0 - yp * (-2.0 - x + y + z) * zp / 8.0,
                            xm * yp * zp / 8.0 + xm * (-2.0 - x + y + z) * zp / 8.0,
                            xm * yp * (-2.0 - x + y + z) / 8.0 + xm * yp * zp / 8.0],

                           [-x * ym * zm / 2.0,
                            -omx2 * zm / 4.0,
                            -omx2 * ym / 4.0],

                           [omy2 * zm / 4.0,
                            -xp * y * zm / 2.0,
                            -xp * omy2 / 4.0],

                            [-x * yp * zm / 2.0,
                             omx2 * zm / 4.0,
                             -omx2* yp / 4.0],
    
                            [-omy2 * zm / 4.0,
                             -xm * y * zm / 2.0,
                             -xm * omy2 / 4.0],
    
                            [-x * ym * zp / 2.0,
                             -omx2 * zp / 4.0,
                             omx2* ym / 4.0],
    
                            [omy2 * zp / 4.0,
                             -xp * y * zp / 2.0,
                             xp * omy2 / 4.0],
    
                            [-x * yp * zp / 2.0,
                             omx2 * zp / 4.0,
                             omx2* yp / 4.0],
    
                            [-omy2 * zp / 4.0,
                             -xm * y * zp / 2.0,
                             xm * omy2 / 4.0],
    
                            [-ym * omz2 / 4.0,
                             -xm * omz2 / 4.0,
                             -xm * ym * z / 2.0],
    
                            [ym * omz2 / 4.0,
                             -xp * omz2 / 4.0,
                             -xp * ym * z / 2.0],
    
                            [yp * omz2 / 4.0,
                             xp * omz2 / 4.0,
                             -xp * yp * z / 2.0],
                            
                            [-yp * omz2 / 4.0,
                             xm * omz2 / 4.0,
                             -xm * yp * z / 2.0]], dtype=np.float64)
def fill_GPs():
        xyzw=np.empty([8,4])
        c=1./np.sqrt(3.)
        xyzw[:,0]=[-c,c,-c,c,-c,c,-c,c]
        xyzw[:,1]=[-c,-c,c,c,-c,-c,c,c]
        xyzw[:,2]=[-c,-c,-c,-c,c,c,c,c]
        xyzw[:,3]=[1.,1.,1.,1.,1.,1.,1.,1.]
        return xyzw.astype(np.float64)
    
def fillGP_full():
    xyzw=np.empty([27,4])
    c1 = 5/9 ; c2 = 8/9 ; alpha = np.sqrt(3/5)
    wg1 = c1**3; wg2 = c1**2*c2 ; wg3 = c1*c2**2; wg4 = c2**3
    
    xyzw[0:27:3,0] = -alpha
    xyzw[1:27:3,0] = 0.0
    xyzw[2:27:3,0] = alpha

    xyzw[0:3,1]   = -alpha
    xyzw[3:6,1]   =  0.0
    xyzw[6:9,1]   = +alpha
    xyzw[9:12,1] = -alpha
    xyzw[12:15,1] =  0.0
    xyzw[15:18,1] = +alpha
    xyzw[18:21,1] = -alpha
    xyzw[21:24,1] =  0.0
    xyzw[24:27,1] = +alpha

    xyzw[0:9,2]   = -alpha
    xyzw[9:18,2] =  0.0
    xyzw[18:27,2] = +alpha
    
    xyzw[:,3]=[wg1,wg2,wg1,wg2,wg3,wg2,wg1,wg2,
                wg1,wg2,wg3,wg2,wg3,wg4,wg3,wg2,
                wg3,wg2,wg1,wg2,wg1,wg2,wg3,wg2,wg1,wg2,wg1]
        
    return xyzw.astype(np.float64)

def calculate_jacobian_at_gp(gp_natural_coords, element_node_coords_global):
    """
    Calculates the Jacobian matrix J and its determinant det(J) at a Gauss point.

    """
    zeta1, zeta2, zeta3 = gp_natural_coords
    dN_local = get_dN(zeta1, zeta2, zeta3)  # Shape: (20 nodes, 3 natural_dims)

    J = np.zeros((3,3), dtype=np.float64)
    for i in range(3): # global_dim x, y, z
        for j in range(3): # natural_dim zeta1, zeta2, zeta3
            J[i, j] = np.sum(dN_local[:, j] * element_node_coords_global[:, i])

    det_J = np.linalg.det(J)
    return J, det_J, dN_local


def calculate_strain_at_gp(gp_natural_coords, element_node_coords_global, element_nodal_displacements):
    """
    Calculates the small strain vector (6 components) at a Gauss point.
    Strain vector: [eps_x, eps_y, eps_z, gam_xy, gam_yz, gam_zx] (Voigt Notation)
    """
    J, det_J, dN_local = calculate_jacobian_at_gp(gp_natural_coords, element_node_coords_global)

    if det_J <= 1e-12: # 
        raise ValueError(f"Jacobian determinant is non-positive or too small ({det_J:.2e}) at GP {gp_natural_coords}. Element may be distorted.")

    inv_J = np.linalg.inv(J)
    inv_J_T = inv_J.T  

    B_matrix = np.zeros((6, 60), dtype=np.float64) # 6 strain components, 20 nodes * 3 DOF/node

    for k in range(20): 
        dNk_local_row_vec = dN_local[k, :] # [dNk/d_zeta1, dNk/d_zeta2, dNk/d_zeta3]

        dNk_global_xyz = inv_J_T @ dNk_local_row_vec # [dNx/dx, dNx/dy, dNx/dz] for node k
        
        dNx_dx = dNk_global_xyz[0]
        dNx_dy = dNk_global_xyz[1]
        dNx_dz = dNk_global_xyz[2]

        col_start = k * 3
        
        B_matrix[0, col_start]     = dNx_dx
        B_matrix[1, col_start + 1] = dNx_dy
        B_matrix[2, col_start + 2] = dNx_dz

        B_matrix[3, col_start]     = dNx_dy
        B_matrix[3, col_start + 1] = dNx_dx

        B_matrix[4, col_start + 1] = dNx_dz
        B_matrix[4, col_start + 2] = dNx_dy

        B_matrix[5, col_start]     = dNx_dz
        B_matrix[5, col_start + 2] = dNx_dx

    # Ensure nodal displacements is a column vector (60,1)
    element_nodal_displacements_col = element_nodal_displacements.reshape(60, 1)
    
    strain_vector = B_matrix @ element_nodal_displacements_col
    return strain_vector.flatten() 


def plot_c3d20r_element(node_coords, GPs=None,show_labels=True, title="C3D20R Element"):
    """
    Plots a 20-node hexahedral element (C3D20R) with node labels.
                            Nodes 1-8: Corners
                            Nodes 9-20: Midside nodes
    """
    nodes = np.array(node_coords)
    if nodes.shape != (20, 3):
        raise ValueError("node_coords must be a 20x3 array or equivalent list of lists.")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot nodes
    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c='blue', marker='o', s=50, label="Nodes")

    # Add node labels
    if show_labels:
        for i in range(20):
            ax.text(nodes[i, 0] * 1.05, nodes[i, 1] * 1.05, nodes[i, 2] * 1.05,
                    str(i + 1), color='red', fontsize=9)

    
    edges = [
        (0, 8, 1),  # Edge 1-9-2
        (1, 9, 2),  # Edge 2-10-3
        (2, 10, 3), # Edge 3-11-4
        (3, 11, 0), # Edge 4-12-1

        (4, 12, 5), # Edge 5-13-6
        (5, 13, 6), # Edge 6-14-7
        (6, 14, 7), # Edge 7-15-8
        (7, 15, 4), # Edge 8-16-5

        (0, 16, 4), # Edge 1-17-5 (bottom face to top face)
        (1, 17, 5), # Edge 2-18-6
        (2, 18, 6), # Edge 3-19-7
        (3, 19, 7)  # Edge 4-20-8
    ]

    for n1_idx, mid_idx, n2_idx in edges:

        ax.plot([nodes[n1_idx, 0], nodes[mid_idx, 0]],
                [nodes[n1_idx, 1], nodes[mid_idx, 1]],
                [nodes[n1_idx, 2], nodes[mid_idx, 2]], 'k-', lw=1.5)

        ax.plot([nodes[mid_idx, 0], nodes[n2_idx, 0]],
                [nodes[mid_idx, 1], nodes[n2_idx, 1]],
                [nodes[mid_idx, 2], nodes[n2_idx, 2]], 'k-', lw=1.5)

    if GPs is not None:
        # Plot Gauss points
        ax.scatter(GPs[:, 0], GPs[:, 1], GPs[:, 2], c='green', marker='x', s=50, label="Gauss Points")
        for i in range(GPs.shape[0]):
            ax.text(GPs[i, 0] * 1.05, GPs[i, 1] * 1.05, GPs[i, 2] * 1.05,
                    f"GP {i + 1}", color='orange', fontsize=9)
    ax.set_title(title)

    x_min, x_max = nodes[:,0].min(), nodes[:,0].max()
    y_min, y_max = nodes[:,1].min(), nodes[:,1].max()
    z_min, z_max = nodes[:,2].min(), nodes[:,2].max()
    
    max_range = np.array([x_max-x_min, y_max-y_min, z_max-z_min]).max() / 2.0

    mid_x = (x_max+x_min) * 0.5
    mid_y = (y_max+y_min) * 0.5
    mid_z = (z_max+z_min) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    ax.legend()
    plt.show()

def from_voigt(s):
    return np.asarray([[s[0],s[3]/2, s[5]/2],
                       [s[3]/2, s[1], s[4]/2],
                       [s[5]/2, s[4]/2, s[2]]], dtype=np.float64)

#%%
default_nodes = get_default_c3d20_nodes()
plot_c3d20r_element(default_nodes, title="Default Canonical C3D20R Element")
GPs= fill_GPs()
plot_c3d20r_element(default_nodes, GPs,title="Default Canonical C3D20R Element")
GPs_full = fillGP_full()
plot_c3d20r_element(default_nodes, GPs_full, title="Default Canonical C3D20R Element with Full Gauss Points")
# %%
nodal_displacements = np.zeros(60, dtype=np.float64)
nodal_displacements[1*3 + 0] = 0.1 # ux for node 2
nodal_displacements[2*3 + 0] = 0.2 # ux for node 3
nodal_displacements[6*3 + 1] = 0.005 # uy for node 7
strain_gp1 = calculate_strain_at_gp(GPs[0,:3], default_nodes, nodal_displacements)

# %%
print(from_voigt(strain_gp1))
# %%
plot_c3d20r_element(default_nodes+nodal_displacements.reshape((20,3)), title="Default Canonical C3D20R Element")
# %%
