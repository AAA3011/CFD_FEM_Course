function [quadN_conn, xCoordN, yCoordN] = quad4_to_quadHigh(quad4_conn, xCoord, yCoord, order)
% QUAD4_TO_QUADN Convert Quad4 mesh to Q_p (order) elements
%   [quadN_conn, xCoordN, yCoordN] = quad4_to_quadN(quad4_conn, xCoord, yCoord, order)
% ... (same header as before) ...

    % --- validate inputs -------------------------------------------------
    if nargin < 4
        error('Must provide quad4_conn, xCoord, yCoord, and order (p).');
    end
    p = double(order);
    if p < 1 || p ~= floor(p)
        error('order must be an integer >= 1');
    end

    nelem = size(quad4_conn,1);
    nnode_original = length(xCoord(:));

    % Number of nodes per Q_p element
    nodes_per_elem = (p+1)^2;

    % Pre-estimate max new nodes to preallocate (safe upper bound)
    max_new_nodes = nnode_original + nelem*( (p-1)^2 + 4*(p-1) + 1 );
    xCoordN = nan(max_new_nodes,1);
    yCoordN = nan(max_new_nodes,1);
    % copy originals
    xCoordN(1:nnode_original) = xCoord(:);
    yCoordN(1:nnode_original) = yCoord(:);
    current_node = nnode_original + 1;

    % initialize connectivity
    quadN_conn = zeros(nelem, nodes_per_elem);

    % edge-map to reuse edge nodes: key = 'min-max-k' where k indexes position along edge
    edge_map = containers.Map('KeyType','char','ValueType','double');

    fprintf('Converting %d Quad4 elements to Q_%d (%d nodes/elem)...\n', nelem, p, nodes_per_elem);

    % reference corner local coords (xi,eta) for corners in the same CCW order as input:
    % corner 1 -> (0,0), 2 -> (1,0), 3 -> (1,1), 4 -> (0,1)
    corner_ref = [0,0; 1,0; 1,1; 0,1];

    % loop elements
    for elem = 1:nelem
        corners = quad4_conn(elem, :);   % indices to original nodes (1x4), CCW

        % store corner node indices into connectivity first (positions 1..4)
        quadN_conn(elem,1:4) = corners;

        % convenience: corner coordinates
        x1 = xCoord(corners(1)); y1 = yCoord(corners(1));
        x2 = xCoord(corners(2)); y2 = yCoord(corners(2));
        x3 = xCoord(corners(3)); y3 = yCoord(corners(3));
        x4 = xCoord(corners(4)); y4 = yCoord(corners(4));

        % Build list of (xi,eta) and decide whether node is corner/edge/interior
        % We'll assign ordering as:
        %  1..4 corners (given)
        %  then edges: edge 1 (1->2) nodes at xi = k/p, eta=0 (k=1..p-1)
        %              edge 2 (2->3) xi=1, eta=k/p
        %              edge 3 (3->4) xi = 1 - k/p, eta=1  (to go 3->4 CCW)
        %              edge 4 (4->1) xi=0, eta = 1 - k/p
        %  then interior nodes: now ordered in CCW concentric rings.

        % Keep a local position pointer for connectivity positions
        pos = 5; % next position after corners (1..4)

        % --- edges -------------------------------------------------------
        % define edge endpoint indices in the element (index into 'corners')
        edge_pairs = [1,2; 2,3; 3,4; 4,1];

        for edge = 1:4
            idxA = edge_pairs(edge,1);
            idxB = edge_pairs(edge,2);
            nodeA_idx = corners(idxA);
            nodeB_idx = corners(idxB);

            for k = 1:(p-1)
                % local param t along the edge from nodeA to nodeB
                t = k / p;

                % compute reference (xi,eta) depending on edge, preserving CCW local direction
                switch edge
                    case 1 % 1->2 : eta=0, xi=t
                        xi = t; eta = 0;
                    case 2 % 2->3 : xi=1, eta=t
                        xi = 1; eta = t;
                    case 3 % 3->4 : eta=1, xi = 1 - t  (3->4)
                        xi = 1 - t; eta = 1;
                    case 4 % 4->1 : xi=0, eta = 1 - t  (4->1)
                        xi = 0; eta = 1 - t;
                end

                % Determine unique edge key so two adjacent elements share the same node.
                % Key uses sorted endpoints plus a k_key that counts from min_node->max_node.
                minnode = min(nodeA_idx, nodeB_idx);
                maxnode = max(nodeA_idx, nodeB_idx);
                if nodeA_idx < nodeB_idx
                    k_key = k; % same orientation as min->max
                else
                    k_key = p - k; % reversed orientation
                end
                key = sprintf('%d-%d-%d', minnode, maxnode, k_key);

                if isKey(edge_map, key)
                    nid = edge_map(key);
                else
                    % compute physical coordinates via bilinear mapping
                    x_phys = (1 - xi)*(1 - eta)*x1 + xi*(1 - eta)*x2 + xi*eta*x3 + (1 - xi)*eta*x4;
                    y_phys = (1 - xi)*(1 - eta)*y1 + xi*(1 - eta)*y2 + xi*eta*y3 + (1 - xi)*eta*y4;

                    nid = current_node;
                    xCoordN(nid) = x_phys;
                    yCoordN(nid) = y_phys;
                    edge_map(key) = nid;
                    current_node = current_node + 1;
                end

                quadN_conn(elem, pos) = nid;
                pos = pos + 1;
            end
        end

        % --- interior nodes (now created and ordered in CCW concentric rings) ----------
        m = p - 1;  % interior grid size (m x m)
        if m > 0
            % create interior node id grid: indices i=1..m (xi left->right), j=1..m (eta bottom->top)
            interior_ids = zeros(m,m);

            % First create all interior nodes and store their ids in the grid
            for j = 1:m         % eta = j/p (bottom->top)
                eta = j / p;
                for i = 1:m     % xi = i/p (left->right)
                    xi = i / p;
                    % bilinear map for interior node
                    x_phys = (1 - xi)*(1 - eta)*x1 + xi*(1 - eta)*x2 + xi*eta*x3 + (1 - xi)*eta*x4;
                    y_phys = (1 - xi)*(1 - eta)*y1 + xi*(1 - eta)*y2 + xi*eta*y3 + (1 - xi)*eta*y4;

                    nid = current_node;
                    xCoordN(nid) = x_phys;
                    yCoordN(nid) = y_phys;
                    current_node = current_node + 1;

                    interior_ids(i,j) = nid; % store by (i=xi, j=eta)
                end
            end

            % Now output the interior node ids in CCW concentric rings:
            rings = ceil(m/2);
            for r = 1:rings
                left = r;
                right = m - r + 1;
                bottom = r;
                top = m - r + 1;

                if left == right && bottom == top
                    % single center node
                    quadN_conn(elem, pos) = interior_ids(left, bottom);
                    pos = pos + 1;
                else
                    % 1) bottom row, left->right (i = left..right, j = bottom)
                    for i = left:right
                        quadN_conn(elem, pos) = interior_ids(i, bottom);
                        pos = pos + 1;
                    end
                    % 2) right column, bottom+1 .. top (j increases)
                    for j = bottom+1:top
                        quadN_conn(elem, pos) = interior_ids(right, j);
                        pos = pos + 1;
                    end
                    % 3) top row, right-1 .. left (i decreases) only if top != bottom
                    if top ~= bottom
                        for i = (right-1):-1:left
                            quadN_conn(elem, pos) = interior_ids(i, top);
                            pos = pos + 1;
                        end
                    end
                    % 4) left column, top-1 .. bottom+1 (j decreases) only if left ~= right
                    if left ~= right
                        for j = (top-1):-1:(bottom+1)
                            quadN_conn(elem, pos) = interior_ids(left, j);
                            pos = pos + 1;
                        end
                    end
                end
            end
        end

        if pos - 1 ~= nodes_per_elem
            error('Unexpected node count per element constructed for element %d', elem);
        end
    end

    % trim coordinate arrays
    last_node = current_node - 1;
    xCoordN = xCoordN(1:last_node);
    yCoordN = yCoordN(1:last_node);

    fprintf('Conversion complete!\n');
    fprintf('Original nodes: %d\n', nnode_original);
    fprintf('New nodes: %d\n', last_node);
end
