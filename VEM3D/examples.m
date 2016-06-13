%   Poisson problem examples for VEM3D.
%
%
%   Example 1: Solution to 
%
%             -\Delta u = 4 \pi^2 x cos(2\pi y) x \in \Omega,
%                    gD = x \cos(2 \pi y)       x \in \partial \Omega.
%   
%
%   Example 2: Simple 3D problem
%
%              -\Delta u = -y.*e^z(x^2+2),
%
%              with Dirichlet and Neumann boundary conditions.
%              Visualization of contribution to L^2 error.
%
%   Example 3: Point source problem.
%-----------------------------------------------------------------ØSK-2016-

clc; clear; close all;

%   Change to your local version of voronoi3D. voroni3D can be downloaded
%   from git: https://github.com/92runarlb/pebiGridding.git
%   Alternatively, use the pregenerated grid in voronoiGrid.mat.

addpath('~/Documents/master/pebiGridding/voronoi3D/')

ex = 3;

switch ex
    
    case 1
        
        %%  Generate grid
        
        %   Uncomment cartGrid to use regular Cartesian grid. Uncomment
        %   load(.. to use pregenerated voronoi grid.
        
        G = voronoiCube(300,[1,1,1]);
%         G = cartGrid([7,7,7], [1,1,1]);
%         load('voronoiGrid.mat');
        G = computeVEM3DGeometry(G);
        
        %%  Define bounday functions and source term
        
        f  = @(X) 4*pi^2*X(:,1).*cos(2*pi*X(:,2)); % source term.
        gD = @(X) X(:,1).*cos(2*pi*X(:,2));        % Dirichlet bc function.
            
        %%  Set boundary conditions
    
        boundaryEdges = find(any(G.faces.neighbors == 0,2));
        bc = VEM3D_addBC([], boundaryEdges, 'pressure', gD);
        
        %%  Solve problem
        
        k   = 2;    %   Method order.
        sol = VEM3D(G, f, bc, k, 'cellProjectors', true);
        
        %%  Plot results
        
        %   Visualize solution using the first moments over the faces. Ball
        %   with radius r, center c carved out  to see the solution
        %   away from Dirichlet boundary.
        
        Kc = G.cells.centroids;
        cells = 1:G.cells.num;
        r = .7; c = [1,1,0];
        cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
        faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
        faces = G.cells.faces(faceNum);
                
        plotFaces(G, faces, sol.faceMoments(faces))
        colorbar;
        axis equal; axis([0,1,0,1,0,1]);
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
    case 2

        %%  Generate grid
        
        G = voronoiCube(200,[1,1,1]);
        G = computeVEM3DGeometry(G);
        
        %%  Define bounday functions and source term
        
        f  = @(X) -X(:,2).*exp( X(:,3) ).*( X(:,1).^2 + 2 );
        gD = @(X) X(:,1).^2.*X(:,2).*exp( X(:,3) );
        gN = @(X) X(:,1).^2.*X(:,2).*exp( X(:,3) );
            
        %%  Set boundary conditions
    
        boundaryEdges = find(any(G.faces.neighbors == 0,2));
        tol   = 1e-6;
        isNeu = abs(G.faces.centroids(boundaryEdges,3)-1) < tol;
        bc    = VEM3D_addBC([], boundaryEdges(~isNeu), 'pressure', gD);
        bc    = VEM3D_addBC(bc, boundaryEdges(isNeu) , 'flux'    , gN);
        
        %%  Solve and calculate error
        
        %   Calculate solutions. cellProjectors = true in order to stor
        %   projection operators \Pi^\nabla for each cell to grid G. used
        %   for approximating L^2 error. faceAverages = true for first
        %   order in order to visualize solution.
        
        k            = 1;
        [sol1, G]    = VEM3D(G, f, bc, k, 'cellProjectors', true, ...
                                          'faceAverages'  , true);
        l2Err1       = l2Error3D(G, sol1, gD ,k);
        [l2Err1, I1] = sort(l2Err1, 'descend');
        
        k            = 2;
        [sol2, G]    = VEM3D(G, f, bc, k, 'cellProjectors', true);
        l2Err2       = l2Error3D(G, sol2, gD ,k);
        [l2Err2, I2] = sort(l2Err2, 'descend');
        
        %%  Plot results
        
        %   Plots results and the 10 cells with largest contribution to L^2
        %   error.
        
        Kc    = G.cells.centroids;
        cells = 1:G.cells.num;
        r = .7; c = [1,1,0];
        cells   = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
        faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
        faces   = G.cells.faces(faceNum);
                
        subplot(2,2,1)
        plotFaces(G, faces, sol1.faceMoments(faces));
        axis equal; axis([0,1,0,1,0,1]);
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
        subplot(2,2,3)
        plotCellData(G,  l2Err1(1:10), I1(1:10));
        axis equal; axis([0,1,0,1,0,1]);
        colorbar;
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
        subplot(2,2,2)
        plotFaces(G, faces, sol2.faceMoments(faces));
        axis equal; axis([0,1,0,1,0,1]);
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
        subplot(2,2,4)
        plotCellData(G,  l2Err2(1:10), I2(1:10));
        axis equal; axis([0,1,0,1,0,1]);
        colorbar;
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
      
    case 3
        
        %%  Generate grid

        G = voronoiCube(200,[1,1,1]);
        G = computeVEM3DGeometry(G);

        %%  Define bounday functions and source term
        
        %   Find cell with centroid closest to C
        C = [0.5, 0.5, 0.5];
        diff = sum(bsxfun(@minus, G.cells.centroids, C).^2, 2);
        srcCell = find(diff == min(diff));
        C = G.cells.centroids(srcCell(1),:);
        
        %   Set source cell with rate Q = 1
        Q = 1;
        src = addSource([], srcCell(1), Q);
        gD = @(X) 1./(4*pi*sqrt(sum(bsxfun(@minus, X, C).^2,2)));
        
        %%  Set boundary conditions

        boundaryEdges = find(any(G.faces.neighbors == 0,2));
        bc = VEM3D_addBC([], boundaryEdges, 'pressure', gD);

        %%  Solve problem

        k   = 2;    %   Method order.
        sol = VEM3D(G, 0, bc, k, 'src', src, 'faceAverages', true);

        %%  Plot results

        Kc = G.cells.centroids;
        cells = 1:G.cells.num;
        r = .9; c = [1,1,0];
        cells = cells(sum(bsxfun(@minus, Kc, c).^2,2) > r^2);
        faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
        faces = G.cells.faces(faceNum);

        plotFaces(G, faces, sol.faceMoments(faces))
        colorbar;
        axis equal; axis([0,1,0,1,0,1]);
        view([150, 30])
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
                
end