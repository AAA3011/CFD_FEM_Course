function [ind_i, ind_j] = generateLagrangeNodeIndices(order_1D, is_serendipity)
    % Generate indices for Lagrange or serendipity elements
    % Input:
    %   order_1D: polynomial order in 1D
    %   is_serendipity: (optional) true for serendipity elements like Q9
    %                   false for full Lagrange elements (default)
    
    if nargin < 2
        is_serendipity = false;
    end
    
    % Number of nodes in 1D = order + 1
    n_nodes_1D = order_1D + 1;
    
    % Generate 1D node positions (in the order used by N_row function)
    xi_1D = linspace(-1, 1, n_nodes_1D);
    xi_1D = xi_1D([1, n_nodes_1D, 2:n_nodes_1D-1]); % Reordering as in N_row
    
    % Initialize arrays
    nodes = [];
    
    % 1. CORNERS (always first 4 nodes)
    % Bottom-left, Bottom-right, Top-right, Top-left
    corners = [-1, -1; 1, -1; 1, 1; -1, 1];
    
    % 2. EDGES (in counter-clockwise order)
    % Bottom edge (excluding corners)
    if order_1D > 1
        for i = 2:order_1D
            nodes = [nodes; xi_1D(i+1), -1]; % +1 because xi_1D(1) is -1
        end
    end
    
    % Right edge (excluding corners)
    if order_1D > 1
        for i = 2:order_1D
            nodes = [nodes; 1, xi_1D(i+1)];
        end
    end
    
    % Top edge (excluding corners, going right to left)
    if order_1D > 1
        for i = order_1D:-1:2
            nodes = [nodes; xi_1D(i+1), 1];
        end
    end
    
    % Left edge (excluding corners, going top to bottom)
    if order_1D > 1
        for i = order_1D:-1:2
            nodes = [nodes; -1, xi_1D(i+1)];
        end
    end
    
    % 3. INTERIOR (row-major order, bottom to top, left to right)
    % For serendipity elements like Q9, we only add center node
    % For full Lagrange, add all interior nodes
    if order_1D > 1
        if is_serendipity
            % For Q9 (quadratic serendipity): only center node
            if order_1D == 2
                nodes = [nodes; 0, 0]; % Center node for Q9
            else
                % For higher-order serendipity, we need different logic
                % For now, we'll handle Q9 specifically
                warning('Higher-order serendipity elements need special handling');
            end
        else
            % For full Lagrange elements: all interior nodes
            if order_1D > 2
                for j = 2:order_1D
                    for i = 2:order_1D
                        nodes = [nodes; xi_1D(i+1), xi_1D(j+1)];
                    end
                end
            end
        end
    end
    
    % Combine corners and other nodes
    all_nodes = [corners; nodes];
    
    % Find indices for each node
    ind_i = zeros(size(all_nodes, 1), 1);
    ind_j = zeros(size(all_nodes, 1), 1);
    
    for node_idx = 1:size(all_nodes, 1)
        xi = all_nodes(node_idx, 1);
        eta = all_nodes(node_idx, 2);
        
        % Find the index in xi_1D
        [~, ind_i(node_idx)] = min(abs(xi_1D - xi));
        [~, ind_j(node_idx)] = min(abs(xi_1D - eta));
    end
end