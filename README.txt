%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This MATLAB matlab code has been developed as part of my Master's
thesis at the Norwegian University of Science and Technology during
the spring te rm 2016.

	     Øystein Strengehagen Klemetsdal (oystein.klemetsdal@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Voronoi-2D contains several functions for generating unstructured Voronoi
grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus 
of this software is the construction of Pebi grids conforming to faults, 
fractures or other geological structures. 

Voronoi-2D supports two different types of structures.
For the first structure type one wants cell edges to follow the 
structure. For the second type, one wants cell centroids to follow 
the structure. In Voronoi-2D is the first type of structure is exclusively 
called for faults, and the second type for wells. 


Voronoi-2D contains the following functions
Functions:
    compositePebiGrid
    createFaultGridPoints
    createWellGridPoints
    pebiGrid
    removeConflicPoints2
    splitAtInt

- compositePebiGrid is an interface function that can be called to create 
a valid MRST grid structure. It creates a Pebi grid conforming to faults 
and wells. The functions creates a semi-structured grid, by inserting 
voronoi seeds around wells and fractures.

- createFaultGridPoints creates points equiv-distant on each side of the
given faults. 

- createWellGridPoints places points along given well-lines.

- pebiGrid is an interface function that can be called to create a valid 
MRST grid structure. It creates a fully unstructured pebi grid that 
conforms to faults and wells. It uses the software DistMesh to create the 
background grid. DistMesh is a software for creating unstructured delaunay 
triangulations: Per-Olof Persson and Gilbert Strang, "A Simple Mesh 
Generator in MATLAB," SIAM Review Vol. 46 (2) 2004.

- removeConflicPoints2 removes any points from one set that are too close 
to points from another set. This function can be used to remove small or
constricted cells.

- splitAtInt is a function that splits a set of paths at each intersection.
It can be used to split all faults and wells at their intersections.
