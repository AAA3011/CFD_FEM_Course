% newVTK - Write VTK output files.
%
% Description: Write an unstructured VTK file (ASCII) with point data for
% velocity and pressure for visualization in ParaView.
% Inputs:
%   xCoord_vec (column vector): Nodal x-coordinates of the mesh.
%   yCoord_vec (column vector): Nodal y-coordinates of the mesh.
%   un_col (column vector): Nodal u-velocity column vector (N x 1)
%   vn_col (column vector): Nodal v-velocity column vector (N x 1)
%   connectivityMatrix_mat (matrix): Element connectivity matrix (nElements x nodesPerElement)
%   pressureSolution_col (column vector): Nodal pressure column vector (N x 1)
%   filename (variable): Output file name (string).
%   time (variable): Time step (scalar)
% Outputs:
%   none
function newVTK(xCoord_vec, yCoord_vec, un_col, vn_col, connectivityMatrix_mat, pressureSolution_col, filename, time)
    % Prepare data for export
    xNodes = xCoord_vec;
    yNodes = yCoord_vec;
    zNodes = zeros(size(xNodes)); % Create a z-vector of zeros for 2D data
    
    % Combine velocity components into a vector field
    velocityVector = [un_col, vn_col, zeros(length(un_col), 1)];
    
    % Calculate velocity magnitude
    velocityMagnitude = sqrt(un_col.^2 + vn_col.^2);
    
    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Write VTK header with time information
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Flow Simulation - Time: %.6f\n', time);
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n\n');
    
    % Write points
    fprintf(fid, 'POINTS %d float\n', length(xNodes));
    for i = 1:length(xNodes)
        fprintf(fid, '%.6f %.6f %.6f\n', xNodes(i), yNodes(i), zNodes(i));
    end
    fprintf(fid, '\n');
    
    % Write cells (elements)
    numElements = size(connectivityMatrix_mat, 1);
    numNodesPerElement = size(connectivityMatrix_mat, 2);
    
    % Calculate total size: 5 integers per element (1 for count + 4 node indices)
    totalSize = numElements * (1 + numNodesPerElement);
    
    fprintf(fid, 'CELLS %d %d\n', numElements, totalSize);
    for i = 1:numElements
        % VTK uses zero-based indexing, so subtract 1 from node numbers
        nodes = connectivityMatrix_mat(i, :) - 1;
        fprintf(fid, '%d', numNodesPerElement);
        for j = 1:numNodesPerElement
            fprintf(fid, ' %d', nodes(j));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
    
    % Write cell types
    fprintf(fid, 'CELL_TYPES %d\n', numElements);
    if numNodesPerElement == 4
        cellType = 9; % VTK_QUAD
    elseif numNodesPerElement == 3
        cellType = 5; % VTK_TRIANGLE
    else
        cellType = 7; % VTK_POLYGON (generic)
    end
    fprintf(fid, '%d\n', repmat(cellType, numElements, 1));
    fprintf(fid, '\n');
    
    % Write FIELD data (global data - only 1 tuple) BEFORE point data
    fprintf(fid, 'FIELD FieldData 2\n');
    fprintf(fid, 'TIME 1 1 double\n');
    fprintf(fid, '%.6f\n', time);
    fprintf(fid, 'CYCLE 1 1 int\n');
    fprintf(fid, '%d\n', round(time/0.005));
    fprintf(fid, '\n');
    
    % Write point data
    fprintf(fid, 'POINT_DATA %d\n', length(xNodes));
    
    % Write velocity vectors
    fprintf(fid, 'VECTORS velocity float\n');
    for i = 1:length(velocityVector)
        fprintf(fid, '%.6f %.6f %.6f\n', velocityVector(i, 1), velocityVector(i, 2), velocityVector(i, 3));
    end
    fprintf(fid, '\n');
    
    % Write pressure scalar field
    fprintf(fid, 'SCALARS pressure float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i = 1:length(pressureSolution_col)
        fprintf(fid, '%.6f\n', pressureSolution_col(i));
    end
    fprintf(fid, '\n');
    
    % Write velocity magnitude scalar field
    fprintf(fid, 'SCALARS velocity_magnitude float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for i = 1:length(velocityMagnitude)
        fprintf(fid, '%.6f\n', velocityMagnitude(i));
    end
    
    % Close file
    fclose(fid);
end
