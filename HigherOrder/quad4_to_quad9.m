
function [quad9_conn, xCoord9, yCoord9] = quad4_to_quad9(quad4_conn, xCoord, yCoord)
% Corrected conversion from Quad4 to Quad9 with proper node ordering
% Input:
%   quad4_conn - Quad4 connectivity matrix [nelem x 4]
%   xCoord, yCoord - Node coordinates
% Output:
%   quad9_conn - Quad9 connectivity matrix [nelem x 9]
%   xCoord9, yCoord9 - New node coordinates

    nelem = size(quad4_conn, 1);
    nnode_original = length(xCoord);
    
    % Initialize Quad9 connectivity
    quad9_conn = zeros(nelem, 9);
    
    % Copy corner nodes (positions 1, 2, 3, 4)
    quad9_conn(:, 1:4) = quad4_conn;
    
    % Initialize new coordinate arrays
    xCoord9 = xCoord(:);
    yCoord9 = yCoord(:);
    current_node = nnode_original + 1;
    
    % Map to avoid duplicate edge nodes
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    fprintf('Converting %d Quad4 elements to Quad9...\n', nelem);
    
    for elem = 1:nelem
        corner_nodes = quad4_conn(elem, :);
        
        % Define edge pairs for Quad9 node ordering:
        % Edge 1: between node 1-2 (position 5)
        % Edge 2: between node 2-3 (position 6) 
        % Edge 3: between node 3-4 (position 7)
        % Edge 4: between node 4-1 (position 8)
        % Center: position 9
        
        edge_pairs = [1, 2;  % Edge 1 -> position 5
                      2, 3;  % Edge 2 -> position 6
                      3, 4;  % Edge 3 -> position 7
                      4, 1]; % Edge 4 -> position 8
        
        % Process each edge
        for edge = 1:4
            nodeA_idx = corner_nodes(edge_pairs(edge, 1));
            nodeB_idx = corner_nodes(edge_pairs(edge, 2));
            
            % Create unique key for this edge
            edge_key = sprintf('%d-%d', min(nodeA_idx, nodeB_idx), max(nodeA_idx, nodeB_idx));
            
            if isKey(edge_map, edge_key)
                % Use existing midpoint
                quad9_conn(elem, edge+4) = edge_map(edge_key);
            else
                % Calculate midpoint coordinates
                x_mid = (xCoord(nodeA_idx) + xCoord(nodeB_idx)) / 2;
                y_mid = (yCoord(nodeA_idx) + yCoord(nodeB_idx)) / 2;
                
                % Store new midpoint
                xCoord9(current_node) = x_mid;
                yCoord9(current_node) = y_mid;
                quad9_conn(elem, edge+4) = current_node;
                edge_map(edge_key) = current_node;
                current_node = current_node + 1;
            end
        end
        
        % Calculate element center (position 9)
        x_center = mean(xCoord(corner_nodes));
        y_center = mean(yCoord(corner_nodes));
        
        % Store center node
        xCoord9(current_node) = x_center;
        yCoord9(current_node) = y_center;
        quad9_conn(elem, 9) = current_node;
        current_node = current_node + 1;
    end
    
    fprintf('Conversion complete!\n');
    fprintf('Original nodes: %d\n', nnode_original);
    fprintf('New nodes: %d\n', length(xCoord9));
end