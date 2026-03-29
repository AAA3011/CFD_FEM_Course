function [connectivityMatrix_mat, xCoord_vec, yCoord_vec] = structuredMesh(numElementsX, numElementsY, lengthX, lengthY, xMin, yMin)

    %% 1. Node Grid Generation - FIXED
    numNodesX = numElementsX + 1;
    numNodesY = numElementsY + 1;
    
    dX = lengthX / numElementsX;
    dY = lengthY / numElementsY;
    
    % Create coordinate grid
    totalNodes = numNodesX * numNodesY;
    xCoord_vec = zeros(1, totalNodes);
    yCoord_vec = zeros(1, totalNodes);
    
    % CRITICAL: We need to number nodes row by row from BOTTOM to TOP
    % For each row (y direction), number all x positions
    nodeIdx = 1;
    for row = 1:numNodesY  % Start from bottom (yMin)
        currentY = yMin + (row - 1) * dY;
        for col = 1:numNodesX  % Left to right (xMin to xMax)
            currentX = xMin + (col - 1) * dX;
            xCoord_vec(nodeIdx) = currentX;
            yCoord_vec(nodeIdx) = currentY;
            nodeIdx = nodeIdx + 1;
        end
    end
    
    %% 2. Connectivity Matrix - Now [1 2 3 4] will be CCW
    totalNumElements = numElementsX * numElementsY;
    connectivityMatrix_mat = zeros(totalNumElements, 4);
    
    elementIdx = 1;
    for row = 1:numElementsY  % y direction (bottom to top)
        for col = 1:numElementsX  % x direction (left to right)
            % Bottom-left node of current element
            nodeBL = (row - 1) * numNodesX + col;
            
            % Assign nodes for quadrilateral element
            node1 = nodeBL;                  % Bottom-left
            node2 = nodeBL + 1;              % Bottom-right
            node3 = nodeBL + numNodesX + 1;  % Top-right
            node4 = nodeBL + numNodesX;      % Top-left
            
            % This gives [1 2 3 4] = BL→BR→TR→TL which is CCW
            connectivityMatrix_mat(elementIdx, :) = [node1, node2, node3, node4];
            elementIdx = elementIdx + 1;
        end
    end
    
    %% Debug output for verification
    fprintf('Testing with 1x1 mesh:\n');
    if numElementsX == 1 && numElementsY == 1
        fprintf('Coordinates:\n');
        for i = 1:4
            fprintf('Node %d: (%.1f, %.1f)\n', i, xCoord_vec(i), yCoord_vec(i));
        end
        fprintf('Connectivity: %d %d %d %d\n', connectivityMatrix_mat(1,:));
        
        % Verify CCW
        % Compute cross product to check orientation
        x1 = xCoord_vec(1); y1 = yCoord_vec(1);
        x2 = xCoord_vec(2); y2 = yCoord_vec(2);
        x3 = xCoord_vec(3); y3 = yCoord_vec(3);
        x4 = xCoord_vec(4); y4 = yCoord_vec(4);
        
        area = 0.5 * ((x2*y3 - x3*y2) + (x3*y4 - x4*y3) + (x4*y1 - x1*y4) + (x1*y2 - x2*y1));
        if area > 0
            fprintf('Element is CCW (positive area = %.2f)\n', area);
        else
            fprintf('Element is CW (negative area = %.2f)\n', area);
        end
    end
end