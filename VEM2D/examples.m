%   Poisson problem examples for VEM2D.
%
%
%   Example 1: Solution to 
%
%              -\Delta u = 25 \pi^2 x \cos(5 \pi y) x \in \Omega
%                     gD = x \cos(5 \pi y)          x \in \partial \Omega
%   
%               on exotic grid.
%
%   Example 2: Simple 2D problem
%
%              -\Delta u = (4\pi^2-4)*e^{2x}\cos(2\pi y)
%
%              with Dirichlet and Neumann boundary conditions.
%              Visualization of contribution to L^2 error for each cell.
%
%   Example 3: Point source in the middle of the domain, with Neumann and
%              Dirichelt boundary conditions.
%-----------------------------------------------------------------ØSK-2016-

clc; clear; close all;

ex = 3; %   Example number

switch ex
    
    case 1
    
        %%  Load grid
        
        %   Functions sortEdges and computeVEM2DGeometry must be called for
        %   the grid proir to using VEM2D.
        
        load('elephant.mat');
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        %%  Set source term and boundary functions
        
        f  = @(X) 25*pi^2*X(:,1).*cos(5*pi*X(:,2));
        gD = @(X) X(:,1).*cos(5*pi*X(:,2));
        
        %%  Set boundary conditions
        
        bcEdges = find(any(G.faces.neighbors == 0, 2));
        bc      = VEM2D_addBC([], G, bcEdges, 'pressure', gD);
        
        %%  Solve problem
        
        k   = 2;
        sol = VEM2D(G,f,bc,k);
        
        %%  Plot solutions
        
        plotVEM2D(G,sol,2);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        axis equal
        view([0,90])
        colorbar;
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
    
    case 2
        
        %%  Generate grid
        
        %   Unit square generates a grid of n x n polygons. Uncomment
        %   cartGrid to use a regular Cartesian grid instead.
        
        G = unitSquare([10,10],[1,1]);
%         G = cartGrid([10,10], [1,1]);
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        %%  Set source term and boundary functions
        
        %   gD is also the exact solution, and we can use this to compare
        %   with the arrpoximated solutions.
        
        a  = 2;
        f  = @(X) (4*pi^2-a^2)*exp(-a*X(:,1)).*cos(2*pi*X(:,2));
        gD = @(X) exp(-a*X(:,1)).*cos(2*pi*X(:,2));
        gN = @(X) a*gD(X);
        
        %%  Set Boundary conditions
        
        %   Dirichlet conditions are given as 'pressure', Neumann are given
        %   as 'Neumann'.
        
        bE = find(any(G.faces.neighbors == 0,2));   % All boundary edges.
        tol   = 1e-10;
        isNeu = abs(G.faces.centroids(bE,1)) < tol; % Neumann edges.
        bc    = VEM2D_addBC([], G, bE(~isNeu), 'pressure', gD);
        bc    = VEM2D_addBC(bc, G, bE(isNeu) , 'flux'    , gN);
        
        %%  Solve problem and calculate L^2 error
        
        %   Set projectors = true to store cell projection operators
        %   \Pi^\nabla in grid struct G. Used for approximating L^2 error.
        %   l2Error2D approximates the square of the L^2 error for each
        %   cell.

        k         = 1; % Method order.
        [sol1, G] = VEM2D(G, f, bc, k, 'projectors', true);
        l2Err1    = l2Error2D(G, sol1, gD, 1);
        
        k         = 2; % Method order.
        [sol2, G] = VEM2D(G, f, bc, k, 'projectors', true);
        l2Err2    = l2Error2D(G, sol2, gD, 2);
        
        %%  Plot results
        
        subplot(2,2,1);
        plotVEM2D(G,sol1,1);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
        
        subplot(2,2,2);
        plotVEM2D(G,sol2,2);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');
        
        subplot(2,2,3);
        plotCellData(G, l2Err1)
        colorbar;
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$');
        
        subplot(2,2,4);
        plotCellData(G, l2Err2)
        colorbar;
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$');
        
    case 3
        
        %%  Load grid
        
        load singularityGrid.mat;      
        G = sortEdges(G);
        G = computeVEM2DGeometry(G);
        
        %%  Set source term and boundary functions
        
        f = 0;
        C = [0.5, 0.5];
        gD = @(X) -log(sqrt(sum(bsxfun(@minus,X,C).^2,2)))/(2*pi)    ;
        gN = @(X) -(X(:,2)-C(2))./(2*pi*sum(bsxfun(@minus,X,C).^2,2));
        
        %%  Set boundary conditions

        tol = 1e-6;
        xMax = 1; yMax = 1;

        bDir     = find(abs(G.faces.centroids(:,1))      < tol | ...
                        abs(G.faces.centroids(:,1)-xMax) < tol);
        bNeuS    = find(abs(G.faces.centroids(:,2))      < tol);
        bNeuN    = find(abs(G.faces.centroids(:,2)-yMax) < tol);

        bc = VEM2D_addBC([], G, bDir , 'pressure', gD         );
        bc = VEM2D_addBC(bc, G, bNeuS, 'flux'    , @(X) -gN(X));
        bc = VEM2D_addBC(bc, G, bNeuN, 'flux'    , gN         );
        
        %%  Set source cell
        
        Q = 1;
        srcCells = find(G.cells.tag);
        src = addSource([],srcCells(1),Q);
        
        %%  Solve problem
        
        k   = 2;
        sol = VEM2D(G, f, bc, k, 'src', src);
        
        %%  Plot solution
        
        plotVEM2D(G,sol,k);
        axis([0,1,0,1]);
        set(gcf, 'defaultTextInterpreter', 'LaTex');
        xlabel('$x$'); ylabel('$y$'); zlabel('$u_h$');   
        
end

