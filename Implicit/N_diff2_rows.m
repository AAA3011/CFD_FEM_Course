% N_diff2_rows - Compute derivatives of Lagrange shape functions.
%
% FILE: N_diff2_rows.m
% DESCRIPTION:
% Compute second derivatives of 1D Lagrange basis functions evaluated at
% given reference coordinates. For linear elements this returns zeros.
%
% Inputs:
%   K (variable): Number of nodes per element side for the solution element (K = order+1).
%   xi_elev (variable): Natural coordinate evaluation points along the xi direction.
% Outputs:
%   none
function N_diff2_row_pages = N_diff2_rows(K,xi_elev)
    N_points = size(xi_elev,3); % number of points at which shape functions are to be computed
    xi_i_row = linspace(-1,1,K);
    xi_i_row = xi_i_row([1,K,2:K-1]); % reordering nodes, because node 2 is always at the end

    % Keep the same output layout as N_diff_rows: 1 x K x N_points
    N_diff2_row_pages = zeros(1,K,N_points);

    % Linear elements have zero second derivatives.
    if K <= 2
        return;
    end

    % For L_k(x) = prod_{j!=k} ((x-x_j)/(x_k-x_j)),
    % L''_k(x) = sum_{i!=k} sum_{m!=k,m!=i} [1/(x_k-x_i) 1/(x_k-x_m)
    %             * prod_{j!=k,j!=i,j!=m} ((x-x_j)/(x_k-x_j))].
    for k = 1:K
        idx_other = [1:k-1, k+1:K];

        for a = 1:numel(idx_other)
            i = idx_other(a);
            coeff_i = 1/(xi_i_row(k)-xi_i_row(i));

            idx_after_i = idx_other;
            idx_after_i(a) = [];

            for b = 1:numel(idx_after_i)
                m = idx_after_i(b);
                coeff = coeff_i/(xi_i_row(k)-xi_i_row(m));

                idx_rest = idx_after_i;
                idx_rest(b) = [];

                term = ones(1,1,N_points);
                for r = 1:numel(idx_rest)
                    j = idx_rest(r);
                    term = term .* ((xi_elev-xi_i_row(j))/(xi_i_row(k)-xi_i_row(j)));
                end

                N_diff2_row_pages(1,k,:) = N_diff2_row_pages(1,k,:) + coeff*term;
            end
        end
    end
end
