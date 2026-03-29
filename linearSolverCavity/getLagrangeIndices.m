function [ind_i, ind_j, all_nodes] = getLagrangeIndices(order_1D)
% Get Lagrange indices (ind_i, ind_j) for COMPLETE CCW ordering
%
% Ordering (UNCHANGED):
% 1. Corners: bottom-left, bottom-right, top-right, top-left
% 2. Edges: bottom → right → top → left (CCW)
% 3. Interior: CCW rings from outer to inner

    K = order_1D + 1;        % nodes per edge
    n_nodes = K^2;

    % --- 1D locations (UNCHANGED) ---
    xi_i_row = linspace(-1,1,K);
    xi_i_row = xi_i_row([1, K, 2:K-1]);   % keep your reordering

    % Preallocate
    ind_i     = zeros(1,n_nodes);
    ind_j     = zeros(1,n_nodes);
    all_nodes = zeros(n_nodes,2);
    assigned  = false(K,K);

    idx = 1;

    % ============================================================
    % 1. CORNERS (CCW) — UNCHANGED
    % ============================================================

    [~, actual_min_idx] = min(xi_i_row);
    [~, actual_max_idx] = max(xi_i_row);

    % BL
    ind_i(idx)=actual_min_idx; ind_j(idx)=actual_min_idx;
    all_nodes(idx,:)=[xi_i_row(actual_min_idx),xi_i_row(actual_min_idx)];
    assigned(actual_min_idx,actual_min_idx)=true; idx=idx+1;

    % BR
    ind_i(idx)=actual_max_idx; ind_j(idx)=actual_min_idx;
    all_nodes(idx,:)=[xi_i_row(actual_max_idx),xi_i_row(actual_min_idx)];
    assigned(actual_min_idx,actual_max_idx)=true; idx=idx+1;

    % TR
    ind_i(idx)=actual_max_idx; ind_j(idx)=actual_max_idx;
    all_nodes(idx,:)=[xi_i_row(actual_max_idx),xi_i_row(actual_max_idx)];
    assigned(actual_max_idx,actual_max_idx)=true; idx=idx+1;

    % TL
    ind_i(idx)=actual_min_idx; ind_j(idx)=actual_max_idx;
    all_nodes(idx,:)=[xi_i_row(actual_min_idx),xi_i_row(actual_max_idx)];
    assigned(actual_max_idx,actual_min_idx)=true; idx=idx+1;

    % ============================================================
    % 2. EDGES (CCW) — UNCHANGED
    % ============================================================

    if order_1D > 1
        all_indices = 1:K;

        % Bottom
        for i = all_indices
            if ~assigned(actual_min_idx,i) && i~=actual_min_idx && i~=actual_max_idx
                ind_i(idx)=i; ind_j(idx)=actual_min_idx;
                all_nodes(idx,:)=[xi_i_row(i),xi_i_row(actual_min_idx)];
                assigned(actual_min_idx,i)=true; idx=idx+1;
            end
        end

        % Right
        for j = all_indices
            if ~assigned(j,actual_max_idx) && j~=actual_min_idx && j~=actual_max_idx
                ind_i(idx)=actual_max_idx; ind_j(idx)=j;
                all_nodes(idx,:)=[xi_i_row(actual_max_idx),xi_i_row(j)];
                assigned(j,actual_max_idx)=true; idx=idx+1;
            end
        end

        % Top
        for i = all_indices(end:-1:1)
            if ~assigned(actual_max_idx,i) && i~=actual_min_idx && i~=actual_max_idx
                ind_i(idx)=i; ind_j(idx)=actual_max_idx;
                all_nodes(idx,:)=[xi_i_row(i),xi_i_row(actual_max_idx)];
                assigned(actual_max_idx,i)=true; idx=idx+1;
            end
        end

        % Left
        for j = all_indices(end:-1:1)
            if ~assigned(j,actual_min_idx) && j~=actual_min_idx && j~=actual_max_idx
                ind_i(idx)=actual_min_idx; ind_j(idx)=j;
                all_nodes(idx,:)=[xi_i_row(actual_min_idx),xi_i_row(j)];
                assigned(j,actual_min_idx)=true; idx=idx+1;
            end
        end
    end

    % ============================================================
    % 3. INTERIOR NODES (FIXED — SAME ORDERING)
    % ============================================================

    n_interior = (order_1D-1)^2;
    if n_interior > 0

        interior_indices = setdiff(1:K,[actual_min_idx,actual_max_idx]);
        interior_indices = sort(interior_indices);

        % -------- Q2 SPECIAL CASE (single center) --------
        if order_1D == 2
            c = interior_indices(1);
            ind_i(idx)=c; ind_j(idx)=c;
            all_nodes(idx,:)=[xi_i_row(c),xi_i_row(c)];
            assigned(c,c)=true; idx=idx+1;

        % -------- Order ≥ 3 (CCW rings, UNCHANGED LOGIC) --------
        else
            n_int = length(interior_indices);
            for ring_level = 1:ceil(n_int/2)

                if ring_level <= floor(n_int/2)
                    ring_start = interior_indices(ring_level);
                    ring_end   = interior_indices(end-ring_level+1);

                    % bottom
                    for i = ring_start:ring_end
                        if ~assigned(ring_start,i)
                            ind_i(idx)=i; ind_j(idx)=ring_start;
                            all_nodes(idx,:)=[xi_i_row(i),xi_i_row(ring_start)];
                            assigned(ring_start,i)=true; idx=idx+1;
                        end
                    end

                    % right
                    for j = ring_start+1:ring_end-1
                        if ~assigned(j,ring_end)
                            ind_i(idx)=ring_end; ind_j(idx)=j;
                            all_nodes(idx,:)=[xi_i_row(ring_end),xi_i_row(j)];
                            assigned(j,ring_end)=true; idx=idx+1;
                        end
                    end

                    % top
                    for i = ring_end:-1:ring_start
                        if ~assigned(ring_end,i)
                            ind_i(idx)=i; ind_j(idx)=ring_end;
                            all_nodes(idx,:)=[xi_i_row(i),xi_i_row(ring_end)];
                            assigned(ring_end,i)=true; idx=idx+1;
                        end
                    end

                    % left
                    for j = ring_end-1:-1:ring_start+1
                        if ~assigned(j,ring_start)
                            ind_i(idx)=ring_start; ind_j(idx)=j;
                            all_nodes(idx,:)=[xi_i_row(ring_start),xi_i_row(j)];
                            assigned(j,ring_start)=true; idx=idx+1;
                        end
                    end
                else
                    % center
                    c = interior_indices(ceil(n_int/2));
                    if ~assigned(c,c)
                        ind_i(idx)=c; ind_j(idx)=c;
                        all_nodes(idx,:)=[xi_i_row(c),xi_i_row(c)];
                        assigned(c,c)=true; idx=idx+1;
                    end
                end
            end
        end
    end

    % ============================================================
    % FINAL CHECK
    % ============================================================

    if idx-1 ~= n_nodes
        error('getLagrangeIndices: Expected %d nodes, assigned %d',n_nodes,idx-1);
    end
end
