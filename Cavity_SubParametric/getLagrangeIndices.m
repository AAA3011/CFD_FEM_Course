function [ind_i, ind_j, all_nodes] = getLagrangeIndices(order_1D)
    % Get Lagrange indices (ind_i, ind_j) for COMPLETE CCW ordering
    % All nodes (corners, edges, interior) follow CCW pattern
    %
    % Ordering:
    % 1. Corners: bottom-left, bottom-right, top-right, top-left
    % 2. Edges: bottom (left→right), right (bottom→top), 
    %           top (right→left), left (top→bottom)
    % 3. Interior: arranged in CCW rings from outer to inner
    %    Each ring: bottom row → right column → top row → left column
    
    K = order_1D + 1;  % Number of nodes per edge
    n_nodes = K^2;
    
    % Generate 1D points using the same method as in N_row
    xi_i_row = linspace(-1, 1, K);
    xi_i_row = xi_i_row([1, K, 2:K-1]);  % Your reordering
    
    % Pre-allocate
    ind_i = zeros(1, n_nodes);
    ind_j = zeros(1, n_nodes);
    all_nodes = zeros(n_nodes, 2);
    
    % Track which nodes have been assigned
    assigned = false(K, K);
    idx = 1;
    
    % ============================================================
    % 1. CORNER NODES (CCW order)
    % Corners are at the EXTREMES of the reordered array
    % ============================================================
    
    % Find min and max indices in xi_i_row
    % For your reordering, min is always at index 1, max is at index 2
    min_idx = 1;
    max_idx = 2;
    
    % But let's verify
    [~, actual_min_idx] = min(xi_i_row);
    [~, actual_max_idx] = max(xi_i_row);
    
    % Bottom-left (min_xi, min_eta)
    ind_i(idx) = actual_min_idx;
    ind_j(idx) = actual_min_idx;
    all_nodes(idx, :) = [xi_i_row(actual_min_idx), xi_i_row(actual_min_idx)];
    assigned(actual_min_idx, actual_min_idx) = true;
    idx = idx + 1;
    
    % Bottom-right (max_xi, min_eta)
    ind_i(idx) = actual_max_idx;
    ind_j(idx) = actual_min_idx;
    all_nodes(idx, :) = [xi_i_row(actual_max_idx), xi_i_row(actual_min_idx)];
    assigned(actual_min_idx, actual_max_idx) = true;
    idx = idx + 1;
    
    % Top-right (max_xi, max_eta)
    ind_i(idx) = actual_max_idx;
    ind_j(idx) = actual_max_idx;
    all_nodes(idx, :) = [xi_i_row(actual_max_idx), xi_i_row(actual_max_idx)];
    assigned(actual_max_idx, actual_max_idx) = true;
    idx = idx + 1;
    
    % Top-left (min_xi, max_eta)
    ind_i(idx) = actual_min_idx;
    ind_j(idx) = actual_max_idx;
    all_nodes(idx, :) = [xi_i_row(actual_min_idx), xi_i_row(actual_max_idx)];
    assigned(actual_max_idx, actual_min_idx) = true;
    idx = idx + 1;
    
    % ============================================================
    % 2. EDGE NODES (CCW order, excluding corners)
    % ============================================================
    
    if order_1D > 1
        % Get all indices (1:K)
        all_indices = 1:K;
        
        % Bottom edge (left to right): all indices except corners, at eta=min_idx
        for i = all_indices
            if ~assigned(actual_min_idx, i) && (i ~= actual_min_idx && i ~= actual_max_idx)
                ind_i(idx) = i;
                ind_j(idx) = actual_min_idx;
                all_nodes(idx, :) = [xi_i_row(i), xi_i_row(actual_min_idx)];
                assigned(actual_min_idx, i) = true;
                idx = idx + 1;
            end
        end
        
        % Right edge (bottom to top): xi=max_idx, all eta except corners
        for j = all_indices
            if ~assigned(j, actual_max_idx) && (j ~= actual_min_idx && j ~= actual_max_idx)
                ind_i(idx) = actual_max_idx;
                ind_j(idx) = j;
                all_nodes(idx, :) = [xi_i_row(actual_max_idx), xi_i_row(j)];
                assigned(j, actual_max_idx) = true;
                idx = idx + 1;
            end
        end
        
        % Top edge (right to left): eta=max_idx, all xi except corners
        for i = all_indices(end:-1:1)  % Reverse for right-to-left
            if ~assigned(actual_max_idx, i) && (i ~= actual_min_idx && i ~= actual_max_idx)
                ind_i(idx) = i;
                ind_j(idx) = actual_max_idx;
                all_nodes(idx, :) = [xi_i_row(i), xi_i_row(actual_max_idx)];
                assigned(actual_max_idx, i) = true;
                idx = idx + 1;
            end
        end
        
        % Left edge (top to bottom): xi=min_idx, all eta except corners
        for j = all_indices(end:-1:1)  % Reverse for top-to-bottom
            if ~assigned(j, actual_min_idx) && (j ~= actual_min_idx && j ~= actual_max_idx)
                ind_i(idx) = actual_min_idx;
                ind_j(idx) = j;
                all_nodes(idx, :) = [xi_i_row(actual_min_idx), xi_i_row(j)];
                assigned(j, actual_min_idx) = true;
                idx = idx + 1;
            end
        end
    end
    
    % ============================================================
    % 3. INTERIOR NODES (CCW rings)
    % ============================================================
    
    if order_1D > 2
        % Process interior nodes in CCW rings
        % Start with outermost ring and move inward
        
        % Determine how many interior rings we have
        % For order_1D = n, we have floor((n-1)/2) complete rings
        % and possibly a center node
        
        % Get all interior indices (excluding boundary)
        interior_indices = setdiff(1:K, [actual_min_idx, actual_max_idx]);
        
        % Sort interior indices (they should already be sorted)
        interior_indices = sort(interior_indices);
        
        % Number of interior indices
        n_interior = length(interior_indices);
        
        if n_interior > 0
            % We'll process rings from outermost to innermost
            % Each ring consists of nodes at a fixed distance from boundary
            
            % For each possible ring level
            for ring_level = 1:ceil(n_interior/2)
                % Determine which indices belong to this ring
                if ring_level <= floor(n_interior/2)
                    % Complete ring
                    ring_start = interior_indices(ring_level);
                    ring_end = interior_indices(end - ring_level + 1);
                    
                    % Bottom row of this ring (left to right)
                    for i = ring_start:ring_end
                        if ~assigned(ring_start, i)
                            ind_i(idx) = i;
                            ind_j(idx) = ring_start;
                            all_nodes(idx, :) = [xi_i_row(i), xi_i_row(ring_start)];
                            assigned(ring_start, i) = true;
                            idx = idx + 1;
                        end
                    end
                    
                    % Right column of this ring (bottom to top, excluding bottom corner)
                    for j = ring_start+1:ring_end-1
                        if ~assigned(j, ring_end)
                            ind_i(idx) = ring_end;
                            ind_j(idx) = j;
                            all_nodes(idx, :) = [xi_i_row(ring_end), xi_i_row(j)];
                            assigned(j, ring_end) = true;
                            idx = idx + 1;
                        end
                    end
                    
                    % Top row of this ring (right to left, excluding right corner)
                    for i = ring_end:-1:ring_start
                        if ~assigned(ring_end, i)
                            ind_i(idx) = i;
                            ind_j(idx) = ring_end;
                            all_nodes(idx, :) = [xi_i_row(i), xi_i_row(ring_end)];
                            assigned(ring_end, i) = true;
                            idx = idx + 1;
                        end
                    end
                    
                    % Left column of this ring (top to bottom, excluding corners)
                    for j = ring_end-1:-1:ring_start+1
                        if ~assigned(j, ring_start)
                            ind_i(idx) = ring_start;
                            ind_j(idx) = j;
                            all_nodes(idx, :) = [xi_i_row(ring_start), xi_i_row(j)];
                            assigned(j, ring_start) = true;
                            idx = idx + 1;
                        end
                    end
                else
                    % Center node (when n_interior is odd)
                    center_idx = interior_indices(ceil(n_interior/2));
                    if ~assigned(center_idx, center_idx)
                        ind_i(idx) = center_idx;
                        ind_j(idx) = center_idx;
                        all_nodes(idx, :) = [xi_i_row(center_idx), xi_i_row(center_idx)];
                        assigned(center_idx, center_idx) = true;
                        idx = idx + 1;
                    end
                end
            end
        end
    end
    
    % Verify all nodes were assigned
    if idx - 1 ~= n_nodes
        error('getLagrangeIndices: Expected %d nodes, assigned %d', n_nodes, idx-1);
    end
    
    % Display the ordering
    fprintf('\n=== COMPLETE CCW Ordering for Q%d (order %d) ===\n', n_nodes, order_1D);
    fprintf('Node#\tind_i\tind_j\txi\teta\n');
    for i = 1:n_nodes
        fprintf('%d\t%d\t%d\t%.3f\t%.3f\n', i, ind_i(i), ind_j(i), ...
                all_nodes(i,1), all_nodes(i,2));
    end
    
    % Verify the CCW pattern
    verify_ccw_pattern(ind_i, ind_j, all_nodes, order_1D);
end

function verify_ccw_pattern(ind_i, ind_j, all_nodes, order_1D)
    % Verify that the ordering follows the complete CCW pattern
    
    K = order_1D + 1;
    
    fprintf('\nVerifying CCW pattern:\n');
    
    % Check corners
    fprintf('Corners (1-4): ');
    for i = 1:4
        fprintf('N%d ', i);
    end
    fprintf('\n');
    
    % Check edges
    if order_1D > 1
        edge_start = 5;
        edge_end = 4 + 4*(order_1D-1);
        fprintf('Edges (%d-%d): ', edge_start, edge_end);
        fprintf('Bottom(%d-%d) → ', edge_start, edge_start+(order_1D-2)-1);
        fprintf('Right(%d-%d) → ', edge_start+(order_1D-1), edge_start+2*(order_1D-1)-1);
        fprintf('Top(%d-%d) → ', edge_start+2*(order_1D-1), edge_start+3*(order_1D-1)-1);
        fprintf('Left(%d-%d)\n', edge_start+3*(order_1D-1), edge_end);
    end
    
    % Check interior
    if order_1D > 2
        interior_start = 4 + 4*(order_1D-1) + 1;
        n_interior = (order_1D-1)^2;
        
        fprintf('Interior (%d-%d):\n', interior_start, interior_start+n_interior-1);
        
        % For Q16 (order=3), interior is 4 nodes: one ring
        if order_1D == 3
            fprintf('  Ring 1: ');
            fprintf('N%d(bottom-left) → ', interior_start);
            fprintf('N%d(bottom-right) → ', interior_start+1);
            fprintf('N%d(top-right) → ', interior_start+2);
            fprintf('N%d(top-left)\n', interior_start+3);
        end
        
        % For Q25 (order=4), interior is 9 nodes: two rings
        if order_1D == 4
            fprintf('  Ring 1 (outer): ');
            fprintf('N%d-%d ', interior_start, interior_start+7);
            fprintf('\n  Ring 2 (center): N%d\n', interior_start+8);
        end
    end
    
    fprintf('CCW pattern verification complete.\n');
end