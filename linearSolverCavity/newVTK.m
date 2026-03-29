function newVTK(xCoord_vec, yCoord_vec, un_cvec, vn_cvec, connectivityMatrix_mat, pressureSolution_cvec, filename, time)
    % Prepare data for export
    xNodes = xCoord_vec;
    yNodes = yCoord_vec;
    zNodes = zeros(size(xNodes)); % Create a z-vector of zeros for 2D data

    % Combine velocity components into a vector field
    velocityVector = [un_cvec, vn_cvec, zeros(length(un_cvec), 1)];

    % Calculate velocity magnitude
    velocityMagnitude = sqrt(un_cvec.^2 + vn_cvec.^2);

    numNodesPerElement = size(connectivityMatrix_mat, 2);
    
    % Check if we need to subdivide higher-order elements (order > 2)
    if numNodesPerElement > 9
        % Subdivide higher-order elements into linear quads for ParaView compatibility
        [xNodes, yNodes, zNodes, connectivityMatrix_mat, velocityVector, pressureSolution_cvec, velocityMagnitude] = ...
            subdivideHigherOrderElements(xCoord_vec, yCoord_vec, un_cvec, vn_cvec, ...
                                        connectivityMatrix_mat, pressureSolution_cvec);
        numNodesPerElement = 4; % After subdivision
    end

    % Use ASCII output (most reliable for ParaView)
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

    % Calculate total size: (1 count + N nodes) per element
    totalSize = numElements * (1 + numNodesPerElement);

    % Reorder connectivity for VTK_BIQUADRATIC_QUAD (Quad9 only)
    vtk_conn = connectivityMatrix_mat;
    if numNodesPerElement == 9
        % VTK Quad9 ordering: corners (0-3), mid-edges (4-7), center (8)
        % Our ordering: corners (1-4), edge1 (5), edge2 (6), edge3 (7), edge4 (8), center (9)
        for e = 1:numElements
            local = connectivityMatrix_mat(e, :);
            % Reorder: [c1 c2 c3 c4 e1 e2 e3 e4 center]
            vtk_conn(e, :) = [local(1), local(2), local(3), local(4), ...
                              local(5), local(6), local(7), local(8), local(9)];
        end
    end

    fprintf(fid, 'CELLS %d %d\n', numElements, totalSize);
    for i = 1:numElements
        nodes = vtk_conn(i, :) - 1;  % Convert to zero-based indexing
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
    elseif numNodesPerElement == 9
        cellType = 28; % VTK_BIQUADRATIC_QUAD (properly ordered)
    else
        % For higher orders (16, 25, etc), use polygon with linear interpolation
        cellType = 7; % VTK_POLYGON
    end
    for i = 1:numElements
        fprintf(fid, '%d\n', cellType);
    end
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
    for i = 1:length(pressureSolution_cvec)
        fprintf(fid, '%.6f\n', pressureSolution_cvec(i));
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

function [xNodes_new, yNodes_new, zNodes_new, conn_new, velocity_new, pressure_new, velMag_new] = ...
    subdivideHigherOrderElements(xCoord, yCoord, u_vec, v_vec, connectivity, pressure)
    % Subdivide higher-order quads into linear quads for ParaView
    
    numElements = size(connectivity, 1);
    nodesPerElem = size(connectivity, 2);
    p = round(sqrt(nodesPerElem)) - 1; % polynomial order
    n = p + 1; % nodes per side
    
    % Each higher-order element subdivides into p×p linear quads
    numSubQuads = p * p;
    conn_new = zeros(numElements * numSubQuads, 4);
    
    % Map original connectivity to grid [i,j] where i,j = 1..n
    nodeCounter = size(xCoord, 1);
    xNodes_new = xCoord;
    yNodes_new = yCoord;
    zNodes_new = zeros(size(xCoord));
    velocity_new = [u_vec, v_vec, zeros(size(u_vec))];
    pressure_new = pressure;
    velMag_new = sqrt(u_vec.^2 + v_vec.^2);
    
    elemCounter = 0;
    
    for elem = 1:numElements
        local = connectivity(elem, :);
        
        % Build node grid (i=xi, j=eta, both 1..n)
        grid = zeros(n, n);
        pos = 1;
        % Corners
        grid(1,1) = local(pos); pos = pos + 1;
        grid(n,1) = local(pos); pos = pos + 1;
        grid(n,n) = local(pos); pos = pos + 1;
        grid(1,n) = local(pos); pos = pos + 1;
        
        % Edges
        for k = 2:(n-1)
            grid(k,1) = local(pos); pos = pos + 1; % bottom
        end
        for k = 2:(n-1)
            grid(n,k) = local(pos); pos = pos + 1; % right
        end
        for k = (n-1):-1:2
            grid(k,n) = local(pos); pos = pos + 1; % top
        end
        for k = (n-1):-1:2
            grid(1,k) = local(pos); pos = pos + 1; % left
        end
        
        % Interior (CCW rings)
        if n > 2
            m = n - 2;
            rings = ceil(m/2);
            for r = 1:rings
                left = r + 1;
                right = n - r;
                bottom = r + 1;
                top = n - r;
                
                if left == right && bottom == top
                    grid(left, bottom) = local(pos);
                    pos = pos + 1;
                else
                    for i = left:right
                        grid(i, bottom) = local(pos);
                        pos = pos + 1;
                    end
                    for j = bottom+1:top
                        grid(right, j) = local(pos);
                        pos = pos + 1;
                    end
                    if top ~= bottom
                        for i = (right-1):-1:left
                            grid(i, top) = local(pos);
                            pos = pos + 1;
                        end
                    end
                    if left ~= right
                        for j = (top-1):-1:(bottom+1)
                            grid(left, j) = local(pos);
                            pos = pos + 1;
                        end
                    end
                end
            end
        end
        
        % Create p×p sub-quads
        for j = 1:p
            for i = 1:p
                elemCounter = elemCounter + 1;
                % Sub-quad corners in grid coordinates
                n1 = grid(i, j);
                n2 = grid(i+1, j);
                n3 = grid(i+1, j+1);
                n4 = grid(i, j+1);
                conn_new(elemCounter, :) = [n1, n2, n3, n4];
            end
        end
    end
    
    % Extend arrays to match new node count (no new nodes, just reuse)
    zNodes_new = zeros(size(xNodes_new));
end