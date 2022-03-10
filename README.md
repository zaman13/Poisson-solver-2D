# Poisson-solver-2D



<p float="left">
<a href = "https://github.com/zaman13/Particle-Swarm-Optimization-Fortran-95/tree/master/Fortran%20codes"> <img src="https://img.shields.io/badge/Language-Python-blue" alt="alt text"> </a>
<a href = "https://github.com/zaman13/Poisson-solver-2D/blob/master/LICENSE"> <img src="https://img.shields.io/github/license/zaman13/Poisson-solver-2D" alt="alt text"></a>
<a href = "https://github.com/zaman13/Poisson-solver-2D/tree/master/Code"> <img src="https://img.shields.io/badge/version-1.3-red" alt="alt text"> </a>
</p>

<p>
Finite difference solution of 2D Poisson equation <img src="https://render.githubusercontent.com/render/math?math=\nabla^2u(x,y) = g(x,y)">


<p float="right">
  <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_3.svg"  width = "440" />
  
  <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_4.svg"  width = "440" /> 
</p>  
  
  
Current version can handle Dirichlet, Neumann, and mixed (combination of Dirichlet and Neumann) boundary conditions:

<img src="https://render.githubusercontent.com/render/math?math=u(x=x_L,y) = u_L">  (Dirichlet left boundary value)

<img src="https://render.githubusercontent.com/render/math?math=u(x=x_R,y) = u_R">  (Dirichlet right boundary value)

<img src="https://render.githubusercontent.com/render/math?math=u(x,y=y_T) = u_T">  (Dirichlet top boundary value)

<img src="https://render.githubusercontent.com/render/math?math=u(x,y=y_B) = u_B">  (Dirichlet bottom boundary value)
  
<img src="https://render.githubusercontent.com/render/math?math=u(x_l<x<x_h,y_l<y<y_h) = u_b">  (Dirichlet interior boundary value) 

<img src="https://render.githubusercontent.com/render/math?math=\frac{du}{dx}(x=x_L,y) = u_L">  (Neumann left boundary value)

<img src="https://render.githubusercontent.com/render/math?math=\frac{du}{dx}(x=x_R,y) = u_R">  (Neumann right boundary value)
  
<img src="https://render.githubusercontent.com/render/math?math=\frac{du}{dy}(x,y=y_T) = u_T">  (Neumann left boundary value)

<img src="https://render.githubusercontent.com/render/math?math=\frac{du}{dy}(x,y=y_B) = u_B">  (Neumann right boundary value)
  
  
</p>

The boundary values themselves can be functions of (x,y).

## Package requirements
  - NumPy 
  - SciPy (sparse matrices, sparse linear algebra) <a href = "https://docs.scipy.org/doc/scipy/reference/sparse.html"> <img src="https://img.shields.io/badge/Pkg-sparse-yellow"> </a> <a href = "https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html"> <img src="https://img.shields.io/badge/Pkg-sparse.linalg-yellow"> </a>


## Version notes
- version 1.3
  - It is now possible to apply Neumann and mixed boundary conditions

- version 1.2
  - It is now possible to define arbitrary Dirichlet boundary points at the interior of the solution domain
  
- version 1.1
  - Fixed a bug regarding the right-hand function
  - Figure size and font size adjusted

- version 1.0 notes
  - Sparse matrix implementation. CSR format (Compressed sparse row matrix) matrix.

## Sample Output
#### Dirichlet boundary conditions at outer walls
Solution of <img src="https://render.githubusercontent.com/render/math?math=\nabla^2u(x,y) = 0"> with boundary conditions <img src="https://render.githubusercontent.com/render/math?math=u(-6,y) = 0.5, u(6,y) = 1.2, u(x,-3) = -0.75, u(x,3) = -1"> is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output.png"  width = "400">

#### Dirichlet boundary conditions at outer wall and inner regions
Solution of <img src="https://render.githubusercontent.com/render/math?math=\nabla^2u(x,y) = 0"> with boundary conditions <img src="https://render.githubusercontent.com/render/math?math=u(-6,y) = 0.5, u(6,y) = 1.2, u(x,-3) = -0.75, u(x,3) = -1, u(1<x<1.4,-0.5<y<0.2)=1.5">  is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_3.svg"  width = "400">

#### Mixed boundary conditions (both Dirichlet and Neumann boundary conditions)
Solution of <img src="https://render.githubusercontent.com/render/math?math=\nabla^2u(x,y) = 0"> with boundary conditions <img src="https://render.githubusercontent.com/render/math?math=\frac{du}{dx}(-6,y) = 0, \frac{du}{dx}(6,y) = 0, u(x,-3) = 0, \frac{du}{dy}(x,3) = 0, u(1<x<1.4,-0.5<y<0.2)=1.5"> (Dirichlet boundary condition on the left wall and in the region 1<x<1.4, -0.5<y<0.2. Neumann boundary conditions on the right, top and bottom walls.) is shown below:

 <img src="https://github.com/zaman13/Poisson-solver-2D/blob/master/Laplace_figure_output_4.svg"  width = "400">
 
 
### References
  - Sparse matrices: https://docs.scipy.org/doc/scipy/reference/sparse.html
  - Sparse matrix linear algebra: https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html
