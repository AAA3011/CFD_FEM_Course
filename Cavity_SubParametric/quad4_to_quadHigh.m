function [high_order_conn, x_new, y_new] = quad4_to_quadHigh(conn, x, y, target_order)
    % Convert quadrilateral elements to higher order with complete CCW ordering
    % including interior nodes in CCW rings/spirals
    %
    % Node ordering:
    % 1. Corner nodes: bottom-left, bottom-right, top-right, top-left
    % 2. Edge nodes: bottom (left→right), right (bottom→top), 
    %                top (right→left), left (top→bottom)
    % 3. Interior nodes: arranged in CCW rings from outer to inner
    %    - Bottom interior row (left→right)
    %    - Right interior column (bottom→top)
    %    - Top interior row (right→left)
    %    - Left interior column (top→bottom)
    %    - Then next inner ring, etc.
    
    if target_order < 1
        error('target_order must be at least 1 (Q4 element)');
    end
    
    % Determine current order
    current_nodes_per_elem = size(conn, 2);
    current_order = round(sqrt(current_nodes_per_elem)) - 1;
    
    if target_order <= current_order
        error('target_order must be greater than current element order');
    end
    
    nelem = size(conn, 1);
    nnode_original = length(x);
    
    % Number of nodes for target order quadrilateral
    n_nodes_per_edge = target_order + 1;
    n_nodes_high = n_nodes_per_edge^2;
    
    % Generate reference points in natural coordinates with complete CCW ordering
    [ref_xi, ref_eta, node_positions] = generate_quad_ccw_rings(target_order);
    
    % Initialize arrays
    high_order_conn = zeros(nelem, n_nodes_high);
    x_new = x(:)';  % Make row vectors
    y_new = y(:)';
    next_node = nnode_original + 1;
    
    % Maps for shared nodes
    edge_node_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    interior_node_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    fprintf('Converting %d quadrilateral elements from order %d to order %d...\n', ...
            nelem, current_order, target_order);
    
    for elem = 1:nelem
        % Get current element nodes
        elem_nodes = conn(elem, :);
        
        % For each new node position in CCW order
        for pos = 1:n_nodes_high
            xi = ref_xi(pos);
            eta = ref_eta(pos);
            
            % Compute physical coordinates using Lagrange interpolation
            [x_phys, y_phys] = interpolate_quad_point(elem_nodes, xi, eta, x, y, current_order);
            
            % Determine if this is a shared node
            if node_positions(pos) == 1  % Corner node
                % Use existing corner node
                if pos == 1  % Bottom-left corner
                    node_idx = elem_nodes(1);
                elseif pos == 2  % Bottom-right corner
                    node_idx = elem_nodes(2);
                elseif pos == 3  % Top-right corner
                    node_idx = elem_nodes(3);
                elseif pos == 4  % Top-left corner
                    node_idx = elem_nodes(4);
                else
                    error('Invalid corner position');
                end
                
            elseif node_positions(pos) == 2  % Edge node
                % Determine which edge and create unique key
                if abs(eta + 1) < 1e-10  % Bottom edge
                    t = (xi + 1) / 2;  % Parameter along edge [0,1]
                    edge_nodes = sort([elem_nodes(1), elem_nodes(2)]);
                    key = sprintf('b-%d-%d-%.10f', edge_nodes(1), edge_nodes(2), t);
                    
                elseif abs(xi - 1) < 1e-10  % Right edge
                    t = (eta + 1) / 2;
                    edge_nodes = sort([elem_nodes(2), elem_nodes(3)]);
                    key = sprintf('r-%d-%d-%.10f', edge_nodes(1), edge_nodes(2), t);
                    
                elseif abs(eta - 1) < 1e-10  % Top edge
                    t = (1 - xi) / 2;  % Reverse direction for CCW
                    edge_nodes = sort([elem_nodes(3), elem_nodes(4)]);
                    key = sprintf('t-%d-%d-%.10f', edge_nodes(1), edge_nodes(2), t);
                    
                elseif abs(xi + 1) < 1e-10  % Left edge
                    t = (1 - eta) / 2;  % Reverse direction for CCW
                    edge_nodes = sort([elem_nodes(4), elem_nodes(1)]);
                    key = sprintf('l-%d-%d-%.10f', edge_nodes(1), edge_nodes(2), t);
                else
                    error('Invalid edge node');
                end
                
                % Check if this edge node already exists
                if isKey(edge_node_map, key)
                    node_idx = edge_node_map(key);
                else
                    node_idx = next_node;
                    edge_node_map(key) = node_idx;
                    x_new(node_idx) = x_phys;
                    y_new(node_idx) = y_phys;
                    next_node = next_node + 1;
                end
                
            else  % Interior node (position == 3)
                % Interior nodes are unique to each element
                key = sprintf('e%d-%.10f-%.10f', elem, xi, eta);
                
                if isKey(interior_node_map, key)
                    node_idx = interior_node_map(key);
                else
                    node_idx = next_node;
                    interior_node_map(key) = node_idx;
                    x_new(node_idx) = x_phys;
                    y_new(node_idx) = y_phys;
                    next_node = next_node + 1;
                end
            end
            
            % Store in connectivity
            high_order_conn(elem, pos) = node_idx;
        end
    end
    
    % Trim coordinate arrays
    x_new = x_new(1:next_node-1);
    y_new = y_new(1:next_node-1);
    
    fprintf('Conversion complete!\n');
    fprintf('Original nodes: %d, New nodes: %d\n', nnode_original, length(x_new));
    fprintf('New element type: Q%d with %d nodes per element\n', ...
            target_order, n_nodes_high);
    
    % Display the CCW ordering pattern
    display_ccw_ordering_info(target_order, ref_xi, ref_eta, node_positions);
end

function [ref_xi, ref_eta, node_positions] = generate_quad_ccw_rings(order)
    % Generate natural coordinates for quadrilateral nodes in complete CCW ordering
    % including interior nodes arranged in CCW rings
    %
    % Ordering:
    % 1. Corners (4)
    % 2. Edge nodes (4*(order-1))
    % 3. Interior nodes in CCW rings:
    %    - First ring: nodes adjacent to edges
    %    - Second ring: next inner ring, etc.
    
    n_nodes_per_edge = order + 1;
    n_nodes = n_nodes_per_edge^2;
    
    % Generate 1D Chebyshev-Gauss-Lobatto points
    i = 0:order;
    xi_1d = -cos(pi * i / order);
    
    % Create full grid
    [XI, ETA] = meshgrid(xi_1d, xi_1d);
    
    % Initialize arrays
    ref_xi = zeros(n_nodes, 1);
    ref_eta = zeros(n_nodes, 1);
    node_positions = zeros(n_nodes, 1);  % 1=corner, 2=edge, 3=interior
    
    % Track which nodes have been assigned
    assigned = false(n_nodes_per_edge, n_nodes_per_edge);
    
    idx = 1;
    
    % 1. CORNER NODES (fixed positions)
    % Bottom-left
    ref_xi(idx) = XI(1, 1);
    ref_eta(idx) = ETA(1, 1);
    node_positions(idx) = 1;
    assigned(1, 1) = true;
    idx = idx + 1;
    
    % Bottom-right
    ref_xi(idx) = XI(1, end);
    ref_eta(idx) = ETA(1, end);
    node_positions(idx) = 1;
    assigned(1, end) = true;
    idx = idx + 1;
    
    % Top-right
    ref_xi(idx) = XI(end, end);
    ref_eta(idx) = ETA(end, end);
    node_positions(idx) = 1;
    assigned(end, end) = true;
    idx = idx + 1;
    
    % Top-left
    ref_xi(idx) = XI(end, 1);
    ref_eta(idx) = ETA(end, 1);
    node_positions(idx) = 1;
    assigned(end, 1) = true;
    idx = idx + 1;
    
    % 2. EDGE NODES (CCW order, excluding corners)
    
    % Bottom edge (left to right, excluding corners)
    for col = 2:n_nodes_per_edge-1
        ref_xi(idx) = XI(1, col);
        ref_eta(idx) = ETA(1, col);
        node_positions(idx) = 2;
        assigned(1, col) = true;
        idx = idx + 1;
    end
    
    % Right edge (bottom to top, excluding corners)
    for row = 2:n_nodes_per_edge-1
        ref_xi(idx) = XI(row, end);
        ref_eta(idx) = ETA(row, end);
        node_positions(idx) = 2;
        assigned(row, end) = true;
        idx = idx + 1;
    end
    
    % Top edge (right to left, excluding corners)
    for col = n_nodes_per_edge-1:-1:2
        ref_xi(idx) = XI(end, col);
        ref_eta(idx) = ETA(end, col);
        node_positions(idx) = 2;
        assigned(end, col) = true;
        idx = idx + 1;
    end
    
    % Left edge (top to bottom, excluding corners)
    for row = n_nodes_per_edge-1:-1:2
        ref_xi(idx) = XI(row, 1);
        ref_eta(idx) = ETA(row, 1);
        node_positions(idx) = 2;
        assigned(row, 1) = true;
        idx = idx + 1;
    end
    
    % 3. INTERIOR NODES (CCW rings from outer to inner)
    
    % Number of interior rings (for order >= 2)
    n_interior_rings = floor((order - 1) / 2);
    
    for ring = 1:n_interior_rings
        % Starting row/col for this ring
        start_idx = ring + 1;
        end_idx = n_nodes_per_edge - ring;
        
        if start_idx > end_idx
            break;  % No more interior rings
        end
        
        % Bottom interior row of this ring (left to right)
        for col = start_idx:end_idx
            if ~assigned(start_idx, col)
                ref_xi(idx) = XI(start_idx, col);
                ref_eta(idx) = ETA(start_idx, col);
                node_positions(idx) = 3;
                assigned(start_idx, col) = true;
                idx = idx + 1;
            end
        end
        
        % Right interior column of this ring (bottom to top)
        for row = start_idx+1:end_idx-1
            if ~assigned(row, end_idx)
                ref_xi(idx) = XI(row, end_idx);
                ref_eta(idx) = ETA(row, end_idx);
                node_positions(idx) = 3;
                assigned(row, end_idx) = true;
                idx = idx + 1;
            end
        end
        
        % Top interior row of this ring (right to left)
        for col = end_idx:-1:start_idx
            if ~assigned(end_idx, col)
                ref_xi(idx) = XI(end_idx, col);
                ref_eta(idx) = ETA(end_idx, col);
                node_positions(idx) = 3;
                assigned(end_idx, col) = true;
                idx = idx + 1;
            end
        end
        
        % Left interior column of this ring (top to bottom)
        for row = end_idx-1:-1:start_idx+1
            if ~assigned(row, start_idx)
                ref_xi(idx) = XI(row, start_idx);
                ref_eta(idx) = ETA(row, start_idx);
                node_positions(idx) = 3;
                assigned(row, start_idx) = true;
                idx = idx + 1;
            end
        end
    end
    
    % 4. CENTER NODE (if odd order, there's a single center)
    center_row = ceil(n_nodes_per_edge / 2);
    center_col = ceil(n_nodes_per_edge / 2);
    
    if ~assigned(center_row, center_col)
        ref_xi(idx) = XI(center_row, center_col);
        ref_eta(idx) = ETA(center_row, center_col);
        node_positions(idx) = 3;
        assigned(center_row, center_col) = true;
        idx = idx + 1;
    end
    
    % Verify all nodes are assigned
    if any(~assigned(:))
        error('Not all nodes were assigned in CCW ordering');
    end
end

function [ref_xi, ref_eta, node_positions] = generate_quad_ccw_spiral(order)
    % Alternative: Generate nodes in a true CCW spiral pattern
    % This creates a continuous spiral from the outer boundary to the center
    
    n_nodes_per_edge = order + 1;
    n_nodes = n_nodes_per_edge^2;
    
    % Generate 1D Chebyshev-Gauss-Lobatto points
    i = 0:order;
    xi_1d = -cos(pi * i / order);
    
    % Create full grid
    [XI, ETA] = meshgrid(xi_1d, xi_1d);
    
    % Initialize arrays
    ref_xi = zeros(n_nodes, 1);
    ref_eta = zeros(n_nodes, 1);
    node_positions = zeros(n_nodes, 1);
    
    % Track which nodes have been assigned
    assigned = false(n_nodes_per_edge, n_nodes_per_edge);
    
    idx = 1;
    direction = 'right';  % Start moving right
    row = 1;
    col = 1;
    
    % Define boundaries
    top = 1;
    bottom = n_nodes_per_edge;
    left = 1;
    right = n_nodes_per_edge;
    
    while idx <= n_nodes
        % Assign current node
        ref_xi(idx) = XI(row, col);
        ref_eta(idx) = ETA(row, col);
        
        % Determine node type
        if (row == 1 && col == 1) || (row == 1 && col == n_nodes_per_edge) || ...
           (row == n_nodes_per_edge && col == n_nodes_per_edge) || (row == n_nodes_per_edge && col == 1)
            node_positions(idx) = 1;  % Corner
        elseif row == 1 || col == n_nodes_per_edge || row == n_nodes_per_edge || col == 1
            node_positions(idx) = 2;  % Edge
        else
            node_positions(idx) = 3;  % Interior
        end
        
        assigned(row, col) = true;
        idx = idx + 1;
        
        % Move to next position in spiral
        switch direction
            case 'right'
                if col < right
                    col = col + 1;
                else
                    direction = 'up';
                    row = row + 1;
                    top = top + 1;
                end
            case 'up'
                if row < bottom
                    row = row + 1;
                else
                    direction = 'left';
                    col = col - 1;
                    right = right - 1;
                end
            case 'left'
                if col > left
                    col = col - 1;
                else
                    direction = 'down';
                    row = row - 1;
                    bottom = bottom - 1;
                end
            case 'down'
                if row > top
                    row = row - 1;
                else
                    direction = 'right';
                    col = col + 1;
                    left = left + 1;
                end
        end
    end
end

function [x_phys, y_phys] = interpolate_quad_point(elem_nodes, xi, eta, x, y, current_order)
    % Interpolate physical coordinates using Lagrange shape functions
    
    n_nodes = length(elem_nodes);
    
    if current_order == 1  % Q4 element
        % Standard bilinear interpolation
        N = zeros(1, 4);
        N(1) = 0.25 * (1 - xi) * (1 - eta);  % Bottom-left
        N(2) = 0.25 * (1 + xi) * (1 - eta);  % Bottom-right
        N(3) = 0.25 * (1 + xi) * (1 + eta);  % Top-right
        N(4) = 0.25 * (1 - xi) * (1 + eta);  % Top-left
        
        x_phys = N * x(elem_nodes)';
        y_phys = N * y(elem_nodes)';
        
    elseif current_order == 2 && n_nodes == 9  % Q9 element
        % Q9 shape functions
        N = zeros(1, 9);
        
        % Corner nodes
        N(1) = 0.25 * xi * eta * (xi - 1) * (eta - 1);  % Bottom-left
        N(2) = 0.25 * xi * eta * (xi + 1) * (eta - 1);  % Bottom-right
        N(3) = 0.25 * xi * eta * (xi + 1) * (eta + 1);  % Top-right
        N(4) = 0.25 * xi * eta * (xi - 1) * (eta + 1);  % Top-left
        
        % Edge nodes
        N(5) = 0.5 * eta * (eta - 1) * (1 - xi^2);      % Bottom edge
        N(6) = 0.5 * xi * (xi + 1) * (1 - eta^2);      % Right edge
        N(7) = 0.5 * eta * (eta + 1) * (1 - xi^2);      % Top edge
        N(8) = 0.5 * xi * (xi - 1) * (1 - eta^2);      % Left edge
        
        % Center node
        N(9) = (1 - xi^2) * (1 - eta^2);               % Center
        
        x_phys = N * x(elem_nodes)';
        y_phys = N * y(elem_nodes)';
        
    else
        % For general higher order elements, use tensor product Lagrange
        n_nodes_1d = current_order + 1;
        
        % Generate 1D points for current order
        i = 0:current_order;
        xi_1d = -cos(pi * i / current_order);
        
        % Compute 1D shape functions
        N_xi = lagrange_basis(xi, xi_1d);
        N_eta = lagrange_basis(eta, xi_1d);
        
        % Assume nodes are stored in tensor product order
        N = zeros(1, n_nodes);
        node_idx = 1;
        
        for j = 1:n_nodes_1d  % eta direction
            for i = 1:n_nodes_1d  % xi direction
                N(node_idx) = N_xi(i) * N_eta(j);
                node_idx = node_idx + 1;
            end
        end
        
        x_phys = N * x(elem_nodes)';
        y_phys = N * y(elem_nodes)';
    end
end

function N = lagrange_basis(x, points)
    % Compute Lagrange basis functions at point x
    n = length(points);
    N = ones(1, n);
    
    for i = 1:n
        for j = 1:n
            if i ~= j
                N(i) = N(i) * (x - points(j)) / (points(i) - points(j));
            end
        end
    end
end

function display_ccw_ordering_info(order, ref_xi, ref_eta, node_positions)
    % Display information about the CCW ordering pattern
    
    n_nodes = (order + 1)^2;
    
    fprintf('\n=== CCW Ordering Pattern (Order %d) ===\n', order);
    fprintf('Total nodes: %d\n', n_nodes);
    fprintf('Corners: 4\n');
    fprintf('Edge nodes: %d\n', 4*(order-1));
    fprintf('Interior nodes: %d\n\n', (order-1)^2);
    
    fprintf('Node ordering sequence:\n');
    fprintf('1-4: Corners (BL, BR, TR, TL)\n');
    
    if order >= 2
        edge_start = 5;
        edge_end = 4 + 4*(order-1);
        fprintf('%d-%d: Edge nodes (Bottom→Right→Top→Left)\n', edge_start, edge_end);
        
        if order >= 3
            interior_start = edge_end + 1;
            fprintf('%d-%d: Interior nodes in CCW rings\n', interior_start, n_nodes);
            
            % Show first few interior nodes
            fprintf('\nFirst few interior nodes (xi, eta):\n');
            interior_count = 0;
            for i = interior_start:min(interior_start+7, n_nodes)
                if node_positions(i) == 3
                    fprintf('  Node %d: (%.3f, %.3f)\n', i, ref_xi(i), ref_eta(i));
                    interior_count = interior_count + 1;
                end
                if interior_count >= 4
                    break;
                end
            end
        end
    end
    fprintf('\n');
end

function visualize_ccw_ordering(order)
    % Visualize the CCW ordering pattern
    
    [ref_xi, ref_eta, node_positions] = generate_quad_ccw_rings(order);
    n_nodes = length(ref_xi);
    
    figure;
    hold on;
    
    % Plot with different colors for different node types
    corners = find(node_positions == 1);
    edges = find(node_positions == 2);
    interiors = find(node_positions == 3);
    
    plot(ref_xi(corners), ref_eta(corners), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot(ref_xi(edges), ref_eta(edges), 'bo', 'MarkerSize', 8);
    plot(ref_xi(interiors), ref_eta(interiors), 'go', 'MarkerSize', 6);
    
    % Draw connecting lines to show the ordering
    for i = 1:n_nodes-1
        plot([ref_xi(i), ref_xi(i+1)], [ref_eta(i), ref_eta(i+1)], 'k-', 'LineWidth', 0.5);
    end
    
    % Add node numbers
    for i = 1:n_nodes
        text(ref_xi(i), ref_eta(i), sprintf('%d', i), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8);
    end
    
    axis equal;
    grid on;
    xlabel('\xi');
    ylabel('\eta');
    title(sprintf('Q%d Element - CCW Node Ordering', order));
    legend({'Corners', 'Edges', 'Interiors'}, 'Location', 'best');
    
    % Add arrows to show direction
    arrow_x = [-0.8, -0.6; 0.6, 0.8; 0.8, 0.6; -0.6, -0.8];
    arrow_y = [-0.8, -0.8; -0.8, -0.6; 0.8, 0.6; 0.6, 0.8];
    
    for i = 1:4
        quiver(arrow_x(i,1), arrow_y(i,1), ...
               arrow_x(i,2)-arrow_x(i,1), arrow_y(i,2)-arrow_y(i,1), ...
               0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1);
    end
    
    hold off;
end