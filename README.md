# Poisson-solver-2D



<p float="left">
<a href = "https://github.com/zaman13/Poisson-solver-2D/tree/master/Code"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Poisson-solver-2D/blob/master/LICENSE"> <img src="https://img.shields.io/github/license/zaman13/Poisson-solver-2D" alt="alt text"></a>
<a href = "https://github.com/zaman13/Poisson-solver-2D/tree/master/Code"> <img src="https://img.shields.io/badge/version-1.5-red" alt="alt text"> </a>
</p>

<p>
Finite difference solution of 2D Poisson equation $\nabla^2u(x,y) = f(x,y)$

Detials about the work can be found in the following tutorial paper: 

Zaman, M.A. "Numerical Solution of the Poisson Equation Using Finite Difference Matrix Operators", Electronics 2022, 11, 2365. https://doi.org/10.3390/electronics11152365

<p float="right">
  <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_3.svg"  width = "440" />
  
  <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_4.svg"  width = "440" /> 
</p>  
  
  
Current version can handle Dirichlet, Neumann, and mixed (combination of Dirichlet and Neumann) boundary conditions:

$u(x=x_L,y) = u_L$ (Dirichlet left boundary value)

$u(x=x_R,y) = u_R$ (Dirichlet left boundary value)

$u(x,y=y_T) = u_T$ (Dirichlet left boundary value)

$u(x,y=y_B) = u_B$ (Dirichlet left boundary value)

$u(x_l < x < x_h, y_l < y < y_h) = u_b$  (Dirichlet interior boundary value) 
  
$\frac{du}{dx}(x=x_L,y) = u_L$  (Neumann left boundary value)

$\frac{du}{dx}(x=x_R,y) = u_R$  (Neumann right boundary value)

$\frac{du}{dy}(x,y=y_T) = u_T$  (Neumann top boundary value)

$\frac{du}{dy}(x,y=y_B) = u_B$  (Neumann bottom boundary value)

  
  
</p>

The boundary values themselves can be functions of (x,y). In addition to the boundaries being at the edge of the solution domain, boundary values imposed on interior regions (i.e. regions surrounded by points where the equation is to be solved) can be also be solved using this code. 

## Package requirements
  - NumPy 
  - SciPy (sparse matrices, sparse linear algebra) <a href = "https://docs.scipy.org/doc/scipy/reference/sparse.html"> <img src="https://img.shields.io/badge/Pkg-sparse-yellow"> </a> <a href = "https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html"> <img src="https://img.shields.io/badge/Pkg-sparse.linalg-yellow"> </a>


## Version notes

- version 1.5.2
    - Created geo_source_definition file. Geometry and source definition functions have been moved there
    - Created utilities_definition file. Moved the plotting (and color) definition functions there  
    
- version 1.5
    - Generalized the boundary definition process. A separate function file is used to find all the outer and inner boundary indices.
      All boundary quantities are in list form. This makes it easier to implement all the boundary operations in one go (rather than
      treating each boundary separately).
      
- version 1.4.4
    - Fixed aliasing problem in contour plot export
    
- version 1.3
  - It is now possible to apply Neumann and mixed boundary conditions

- version 1.2
  - It is now possible to define arbitrary Dirichlet boundary points at the interior of the solution domain
  
- version 1.1
  - Fixed a bug regarding the right-hand function
  - Figure size and font size adjusted

- version 1.0 notes
  - Sparse matrix implementation. CSR format (Compressed sparse row matrix) matrix.


## Graphically Defined Geometry
In addition to algebraically defining the solution domain and the boundary regions, it is possible to import the geometry from bitmap image (bmp) file. Differere colors in the bmp file are taken as different regions (i.e. Dirichlet boundary 1, Direchlet boundary 2, solution domain, Neumann boundary etc.). The different regions can be handled appropriately by defining the color mapping in the code. This feature can be useful when working with complex geometries that are difficult/cumbersome to define algebraically. 

## Sample Output
#### Dirichlet boundary conditions at outer walls
Solution of $\nabla^2u(x,y) = 0$ with boundary conditions $u(-6,y) = 0.5$, $u(6,y) = 1.2$, $u(x,-3) = -0.75$, $u(x,3) = -1$ is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output.png"  width = "400">

#### Dirichlet boundary conditions at outer wall and inner regions
Solution of $\nabla^2u(x,y) = 0$ with boundary conditions $u(-6,y) = 0.5$, $u(6,y) = 1.2$, $u(x,-3) = -0.75$, $u(x,3) = -1$, $u(1 < x < 1.4, -0.5 < y < 0.2)=1.5$  is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_3.svg"  width = "400">

#### Mixed boundary conditions (both Dirichlet and Neumann boundary conditions)
Solution of $\nabla^2u(x,y) = 0$ with boundary conditions $\frac{du}{dx}(-6,y) = 0$, $\frac{du}{dx}(6,y) = 0$, $u(x,-3) = 0$, $\frac{du}{dy}(x,3) = 0$, $u(1 < x < 1.4, -0.5 < y < 0.2)=1.5$ (Dirichlet boundary condition on the left wall and in the region $1 < x< 1.4$, $-0.5 < y < 0.2$. Neumann boundary conditions on the right, top and bottom walls.) is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_4.svg"  width = "400">
 
 
 #### Graphically defined geometry (both Dirichlet and Neumann boundary conditions)
Solution $\nabla^2u(x,y) = 0$ with the following boundary conditions:

Center circular region: $u = 1$

Left and right circular region: $u = -2$

Left and right rectangular region: $u = 2$

All the outer boundaries have Neumann boundary conditions: $du/dx = 0$ (left and right boundary), $du/dy = 0$ (top and bottom boundary)

The results are shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Fig_test_pot.png"  width = "400">
 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Fig_test_field.png"  width = "400">
 
 
### References
  - Zaman, M.A. "Numerical Solution of the Poisson Equation Using Finite Difference Matrix Operators", Electronics 2022, 11, 2365.     https://doi.org/10.3390/electronics11152365
  - Sparse matrices: https://docs.scipy.org/doc/scipy/reference/sparse.html
  - Sparse matrix linear algebra: https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html
