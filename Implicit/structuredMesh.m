% structuredMesh - Create structured quadrilateral mesh.
%
% FILE: structuredMesh.m
% DESCRIPTION:
% Generate a regular structured quadrilateral mesh over a rectangular
% domain and return element connectivity and nodal coordinates.
%
% Inputs:
%   numElementsX (variable): Number of elements in the x-direction.
%   numElementsY (variable): Number of elements in the y-direction.
%   lengthX (variable): Total domain length in the x-direction.
%   lengthY (variable): Total domain length in the y-direction.
%   xMin (variable): Minimum x-coordinate of the mesh domain.
%   yMin (variable): Minimum y-coordinate of the mesh domain.
% Outputs:
%   connectivityMatrix_mat : Element connectivity matrix (nElements x nodesPerElement)
%   xCoord_vec : Nodal x-coordinates of the mesh.
%   yCoord_vec : Nodal y-coordinates of the mesh.
function [connectivityMatrix_mat,xCoord_vec,yCoord_vec] = structuredMesh(numElementsX,numElementsY,lengthX,lengthY,xMin,yMin)
% STRUCTUREDMESH creates a structured quadrilateral mesh on a rectangle.
% It builds the element connectivity and the nodal x and y coordinates for a
% regular grid.
% Inputs: numElementsX, numElementsY, lengthX, lengthY, xMin, yMin.
% Outputs: connectivityMatrix_mat, xCoord_vec, yCoord_vec.

    %% Build the connectivity matrix
    
    numNodesX = numElementsX + 1;
    
    totalNumElements = numElementsX *numElementsY;
    connectivityMatrix_mat = zeros(totalNumElements,4);
    
    indexElementX_vec = 1:numElementsX;
    indexElementY_vec = 1:numElementsY;
    
    node1_mat = indexElementX_vec' + (indexElementY_vec - 1) * numNodesX;
    connectivityMatrix_mat(:,1) = node1_mat(:);
    connectivityMatrix_mat(:,2) = connectivityMatrix_mat(:,1) + 1;
    connectivityMatrix_mat(:,3) = connectivityMatrix_mat(:,2) + numNodesX;
    connectivityMatrix_mat(:,4) = connectivityMatrix_mat(:,3) - 1;
    
    %% Build the coordinate grid
    
    dX = lengthX /  numElementsX;
    dY = lengthY /  numElementsY;
    
    indexNodeX_vec = [indexElementX_vec indexElementX_vec(end)+1];
    indexNodeY_vec = [indexElementY_vec indexElementY_vec(end)+1];
    
    xCoord_mat  = xMin + (dX * ones(length(indexNodeY_vec),1) .* (indexNodeX_vec-1));
    xCoord_mat  = xCoord_mat';
    yCoord_mat  = yMin + (dY * ones(length(indexNodeX_vec),1) .*(indexNodeY_vec-1));
    
    xCoord_vec = xCoord_mat(:)';
    yCoord_vec = yCoord_mat(:)';

    %% Optional mesh plot
    % figure;
    % nodes = [connectivityMatrix_mat connectivityMatrix_mat(:,1)];
    % plot(xCoord_vec(nodes)',yCoord_vec(nodes)','k')
    % title('Mesh of Two Connected Rectangular Domains','Interpreter','latex')
    % xlabel('X','Interpreter','latex')
    % ylabel('Y','Interpreter','latex')
    % axis equal


end
