%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This MATLAB MATLAB code has been developed as part of my master's
thesis at the Norwegian University of Science and Technology during
the spring term of 2016.

	     Øystein Strengehagen Klemetsdal (oystein.klemetsdal@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This repository contains an implementation of the virtual element
method for solving differential equations in two and three dimensions
using the open source MATLAB reservoir simulation toolbox (MRST, see
http://www.sintef.no/MRST). It is implemented for the propose of
testing the method, and comparing it to industry standard reservoir
simulation methods and grids.

VEM2D contains the following main functions:

---------------------------------------------------------------------------
- VEM2D_addBC: Adds boundary condition to existing or empty boundary
  condition struct.
	
- VEM2D_bc: Incorporates boundary conditions in global stiffness
  matrix and global load vector.

- VEM2D_glob: Assembles global stiffness matrix and load vector from
  local ones.

- VEM2D_loc: Computes local stiffness matrices and load vectors.

- VEM2D: Solves a given differential equation using the above
  functions.
---------------------------------------------------------------------------

In addition, it contains some auxiliary functions. See each
individual file for details and syntax.

Examples of usage of the software can be found in VEM2D/examples.m.

VEM3D contains the following main functions:

---------------------------------------------------------------------------
- VEM3D_addBC: Adds boundary condition to existing or empty boundary
  condition struct.
	
- VEM3D_bc: Incorporates boundary conditions in global stiffness
  matrix and global load vector.

- VEM3D_faceProjectors: Calculates the VEM projection operators for
  each face in the grid.

- VEM3D_glob: Assembles global stiffness matrix and load vector from
  local ones.

- VEM3D_loc: Computes local stiffness matrices and load vectors.

- VEM3D: Solves a given differential equation using the above
  functions.
---------------------------------------------------------------------------

In addition, it contains some auxiliary functions. See each
individual file for details and syntax.

Examples of usage of the software can be found in VEM3D/examples.m.

For further details, see:

Øystein Strengehagen Klemetsdal,
The virtual element method as a common framework for finite element
and finite difference methods - Numerical and theoretical analysis,
June 2016,
Masters thesis,
Department of Mathematics,
Norwegian University of Science and Technology.
