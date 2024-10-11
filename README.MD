# Explanation of the Belousov-Zhabotinsky Reaction Simulation

## 1. Reaction-Diffusion Equations

The Belousov-Zhabotinsky reaction is modeled using a system of partial differential equations known as reaction-diffusion equations. We use the Oregonator model, which describes the concentrations of two key species, U and V.

### 1.1 Governing Equations

The evolution of species U and V is described by:

$$
\begin{align}
\frac{\partial U}{\partial t} &= D_U \nabla^2 U + \frac{1}{\epsilon}(U - U^2 - fV\frac{U-q}{U+q}) \\
\frac{\partial V}{\partial t} &= D_V \nabla^2 V + U - V
\end{align}
$$

Where:
- $U$ and $V$ are the concentrations of the two species
- $t$ is time
- $D_U$ and $D_V$ are diffusion coefficients
- $\nabla^2$ is the Laplacian operator
- $\epsilon$, $f$, and $q$ are model parameters

### 1.2 Discretized Form

For numerical simulation, we discretize these equations in space and time:

$$
\begin{align}
\frac{U_{i,j}^{n+1} - U_{i,j}^n}{\Delta t} &= D_U \frac{U_{i+1,j}^n + U_{i-1,j}^n + U_{i,j+1}^n + U_{i,j-1}^n - 4U_{i,j}^n}{(\Delta x)^2} + \frac{1}{\epsilon}(U_{i,j}^n - (U_{i,j}^n)^2 - fV_{i,j}^n\frac{U_{i,j}^n-q}{U_{i,j}^n+q}) \\
\frac{V_{i,j}^{n+1} - V_{i,j}^n}{\Delta t} &= D_V \frac{V_{i+1,j}^n + V_{i-1,j}^n + V_{i,j+1}^n + V_{i,j-1}^n - 4V_{i,j}^n}{(\Delta x)^2} + U_{i,j}^n - V_{i,j}^n
\end{align}
$$

Where:
- $i,j$ are spatial indices
- $n$ is the time index
- $\Delta t$ is the time step
- $\Delta x$ is the spatial step

## 2. Numerical Integration Method

We use the fourth-order Runge-Kutta (RK4) method for time integration, which provides higher accuracy compared to simpler methods like Euler integration.

### 2.1 RK4 Algorithm

For a general ODE $\frac{dy}{dt} = f(t, y)$, the RK4 method is:

$$
\begin{align}
k_1 &= f(t_n, y_n) \\
k_2 &= f(t_n + \frac{\Delta t}{2}, y_n + \frac{\Delta t}{2}k_1) \\
k_3 &= f(t_n + \frac{\Delta t}{2}, y_n + \frac{\Delta t}{2}k_2) \\
k_4 &= f(t_n + \Delta t, y_n + \Delta t k_3) \\
y_{n+1} &= y_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{align}
$$

### 2.2 Application to BZ Reaction

For our system, we apply RK4 to both U and V concentrations:

$$
\begin{align}
U_{i,j}^{n+1} &= U_{i,j}^n + \frac{\Delta t}{6}(k_{1U} + 2k_{2U} + 2k_{3U} + k_{4U}) \\
V_{i,j}^{n+1} &= V_{i,j}^n + \frac{\Delta t}{6}(k_{1V} + 2k_{2V} + 2k_{3V} + k_{4V})
\end{align}
$$

Where each $k$ term is calculated using the discretized reaction-diffusion equations.

## 3. Initial Conditions

The simulation starts with small perturbations around the steady-state concentrations:

$$
\begin{align}
U_{steady} &\approx 1.0 \\
V_{steady} &\approx \frac{1}{q}
\end{align}
$$

Initial concentrations are set as:

$$
\begin{align}
U_{i,j}^0 &= U_{steady} + \delta U_{i,j} \\
V_{i,j}^0 &= V_{steady} + \delta V_{i,j}
\end{align}
$$

Where $\delta U_{i,j}$ and $\delta V_{i,j}$ are small random perturbations.

## 4. Boundary Conditions

We use periodic boundary conditions:

$$
\begin{align}
U_{0,j} &= U_{N,j}, \quad U_{N+1,j} = U_{1,j} \\
U_{i,0} &= U_{i,N}, \quad U_{i,N+1} = U_{i,1} \\
V_{0,j} &= V_{N,j}, \quad V_{N+1,j} = V_{1,j} \\
V_{i,0} &= V_{i,N}, \quad V_{i,N+1} = V_{i,1}
\end{align}
$$

Where $N$ is the grid size.

## 5. Visualization

The concentrations U and V are visualized as 2D color maps, where the color intensity represents the concentration level at each grid point $(i,j)$ at a given time step $n$.

This mathematical framework forms the basis of our BZ reaction simulation, capturing the complex spatio-temporal dynamics of this fascinating chemical system.