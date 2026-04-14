# FEniCSx Homework Assignment

## Question 1: Heat Equation with Time Stepping Methods in a 2D Plate

### Problem Setup
Consider a 2D trapezoidal plate $\Omega$ defined by the vertices:

A=(0,0), B=(2,0), C=(1.5,2), D=(0.5,2) (dimesnions in m).


The transient heat conduction problem is:

$$\rho c \frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = 0 \text{ in } \Omega \times (0,T]$$
$$u(x,y,0) = u_0(x,y) \text{ in } \Omega_0 $$
$$u(x,y,t) = 100 \text{ on } \partial\Omega \times (0,T]$$

where:

- $\Omega_0$ is the subdomain defined by $x=[0.5,1.5]$, $y=[0.5,1.5]$ on which we apply the initial condition

- $\partial\Omega$ is the boundary of the domain $\Omega$

- $\rho = 2700 \text{ kg/m}^3$ is the density

- $c = 500 \text{ J/(kg·K)}$ is the specific heat capacity

- $k = 45 \text{ W/(m·K)}$ is the thermal conductivity

- $T$ is the total simulation time (e.g., $T = 10 \text{ s}$)

- $u_0(x,y) = 100+\exp(-a(x-1)^2-b^(y-1)^2) \text{ K}$ (initial temperature distribution) with $a = 10$ and $b = 5$

### Tasks

1. **Derive the weak formulation of the heat equation.**
   

2. **Derive the time-discretized weak form using the time discretization shown in class for:**

   - $\theta=0$ 

   - $\theta=1$ 
   
   - $\theta=0.5$ 

3. **Implement the FEniCSx solution using the provided template for each time-stepping method. Use P1 elements for spatial discretization.**

4. **For each time-stepping method:**
   - Solve with the following different time steps: $\Delta t = 2.5\text{ s}$, $1.25\text{ s}$, $0.75\text{ s}$ and $0.25\text{ s}$

   - Plot the temperature distribution at $t = 2.5\text{ s}$, $t = 5\text{ s}$, and $t = T$ along the line $x = 1$

   - Compare the solutions obtained for different time stepping schemes and time steps. 
   

5. **Analyze and discuss:**
   
   - The stability of each method based on your numerical results
   
   - The computational efficiency and accuracy tradeoff between methods
   

## Question 2: Linear Elasticity with Distributed Load and Different Element Types

### Problem Setup

Consider a rectangular domain $\Omega = [0,L] \times [0,H]$ representing a simply supported beam with length $L = 10$ and height $H = 1$ (dimensions are in meters). The beam is supported at the bottom corners $(x = 0, y = 0)$ and $(x = L, y = 0)$. 

A non-uniform distributed load is applied on the top edge given by:

$$p(x) = p_0 \cdot \sin(\pi x/L)$$

where $p_0 = 2 \text{ MPa}$ is the maximum load intensity.

The governing equations for linear elasticity are:

$$-\nabla \cdot \sigma(\mathbf{u}) = \mathbf{f} \text{ in } \Omega$$
$$u_y = 0 \text{ at points } (0,0) \text{ and } (L,0)$$
$$u_x = 0 \text{ at point } (0,0) $$
$$\sigma(\mathbf{u}) \cdot \mathbf{n} = (0, -p(x)) \text{ on the top edge } (y = H)$$
$$\sigma(\mathbf{u}) \cdot \mathbf{n} = \mathbf{0} \text{ on other edges}$$

where:

- $\sigma(\mathbf{u}) = \lambda(\nabla \cdot \mathbf{u})\mathbf{I} + 2\mu\epsilon(\mathbf{u})$ is the stress tensor

- $\epsilon(\mathbf{u}) = (\nabla \mathbf{u} + (\nabla \mathbf{u})^T)/2$ is the strain tensor

- $\mathbf{f} = (0, -\rho g)$ with $\rho = 1500 \text{ kg/m}^3$ 

- $\lambda = \frac{E \cdot \nu}{(1+\nu)(1-2\nu)}$ and $\mu = \frac{E}{2(1+\nu)}$ are the Lamé parameters

- $E = 2.1 \times 10^5 \text{ MPa}$ is Young's modulus

- $\nu = 0.3$ is Poisson's ratio

### Tasks

1. **Derive the weak formulation for the linear elasticity problem with the non-uniform distributed load.**

   - Start with the strong form and multiply by a test function

   - Apply integration by parts

   - Incorporate the boundary conditions

   - Obtain the final weak form

2. **Derive an expression for the strain energy and the work done by the external forces.**

   - Strain energy: $U = \frac{1}{2} \int_{\Omega} \sigma(\mathbf{u}):\epsilon(\mathbf{u}) \, d\Omega$

   - Work done by external forces: $W = \int_{\Omega} \mathbf{f} \cdot \mathbf{u} \, d\Omega + \int_{\partial\Omega} \mathbf{g} \cdot \mathbf{u} \, ds$

3. **Implement the FEniCSx solution using the provided template with:**

   - P1 elements (linear)

   - P2 elements (quadratic)

4. **For each element type:**

   - Solve the problem with four different mesh refinement levels (start with 20×4 elements and refine by doubling elements in each direction)

   - Calculate the maximum deflection along the top edge

   - Calculate the total strain energy

   - Visualize the displacement magnitude and von Mises stress

   - Compute the reaction forces at the support points

5. **Analyze and discuss:**

   - The order of convergence in displacement and energy for each element type

   - How well different element orders capture the stress concentrations near the supports

   - Compare your numerical results with the analytical solution. 

   
