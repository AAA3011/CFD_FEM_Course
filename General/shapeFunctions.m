% shapeFunctions - Evaluate shape functions and natural derivatives.
%
%
% Supports multiple element types (Quad4, Quad8, Tri6, etc.) and evaluates
% basis functions and their natural derivatives at provided natural coords.
%
% Inputs:
%   n_ElementType (variable): Number of nodes in one element direction for the current element type.
%   coords_row_pages (variable):Gauss point natural coordinates arranged page-wise.
% Outputs:
%   N_row_points_pages : Shape function values evaluated at Gauss points, stored page-wise.
%   N_diff_rows_points_pages : Natural derivatives of shape functions at Gauss points.
function [N_row_points_pages,N_diff_rows_points_pages]=shapeFunctions(n_ElementType,coords_row_pages)
% Modified from the following
% Created in Jan 2021 by CUFE AER-2021 Seniors:
% Mahmoud Ahmed Moustafa, Mohamed Ahmed Ali, Ahmed Elrawy & Osama Mamdouh
% For our Graduation Project, under the instruction and supervision of Dr. Ahmed Rashed

N_points=size(coords_row_pages,3);

xi_elev=coords_row_pages(1,1,:);
if size(coords_row_pages,2)>1
    eta_elev=coords_row_pages(1,2,:);
end

if size(coords_row_pages,2)>2
    zeta_elev=coords_row_pages(1,3,:);
end

switch n_ElementType
    case 0   %Truss
        N_row_points_pages=[1-xi_elev, 1+xi_elev]/2;
        N_diff_rows_points_pages=repmat([-1, 1]/2,1,1,N_points);

    case 1   %Quad4
        N_row_points_pages=[(1+xi_elev).*(1-eta_elev) ...
                            (1+xi_elev).*(1+eta_elev) ...
                            (1-xi_elev).*(1+eta_elev) ...
                            (1-xi_elev).*(1-eta_elev)]/4;

        N_diff_rows_points_pages=[  1-eta_elev, 1+eta_elev, -(1+eta_elev), -(1-eta_elev)
                                  -(1+xi_elev),  1+xi_elev,     1-xi_elev,  -(1-xi_elev)]/4;

    case 2   %Quad8
        N_row_points_pages=[-1/4*(1+xi_elev).*(1-eta_elev).*(1-xi_elev+eta_elev)  ...
                             1/2*(1+xi_elev).*(1-eta_elev.*eta_elev) ...
                            -1/4*(1+xi_elev).*(1+eta_elev).*(1-xi_elev-eta_elev) ...
                             1/2*(1-xi_elev.*xi_elev).*(1+eta_elev) ...
                            -1/4*(1-xi_elev).*(1+eta_elev).*(1+xi_elev-eta_elev) ...
                             1/2*(1-xi_elev).*(1-eta_elev.*eta_elev) ...
                            -1/4*(1-xi_elev).*(1-eta_elev).*(1+xi_elev+eta_elev) ...
                             1/2*(1-xi_elev.*xi_elev).*(1-eta_elev)];

        N_diff_rows_points_pages=permute([ ((eta_elev-2.*xi_elev).*(eta_elev-1))/4,  ((xi_elev+1).*(2*eta_elev-xi_elev))/4
                                           1/2-eta_elev.^2/2,            -eta_elev.*(xi_elev+1)
                                          ((eta_elev+2*xi_elev).*(eta_elev+1))/4,  ((2*eta_elev+xi_elev).*(xi_elev+1))/4
                                           -xi_elev.*(eta_elev+1),             1/2-xi_elev.^2/2
                                           -((eta_elev-2*xi_elev).*(eta_elev+1))/4, -((xi_elev-1).*(2*eta_elev-xi_elev))/4
                                            eta_elev.^2/2-1/2,             eta_elev.*(xi_elev-1)
                                           -((eta_elev+2*xi_elev).*(eta_elev-1))/4, -((2*eta_elev+xi_elev).*(xi_elev-1))/4
                                            xi_elev.*(eta_elev-1),             xi_elev.^2/2-1/2],[2,1,3]);
        
    case 3   %Quad9
        N_row_points_pages=permute([0.25*xi_elev.*(xi_elev+1).*eta_elev.*(eta_elev-1)
                 -0.5*xi_elev.*(xi_elev+1).*(eta_elev.^2-1)
               0.25*xi_elev.*(xi_elev+1).*eta_elev.*(eta_elev+1)
                -0.5*(xi_elev.^2-1).*eta_elev.*(eta_elev+1)
               0.25*xi_elev.*(xi_elev-1).*eta_elev.*(eta_elev+1)
                 -0.5*xi_elev.*(xi_elev-1).*(eta_elev.^2-1)
               0.25*xi_elev.*(xi_elev-1).*eta_elev.*(eta_elev-1)
                -0.5*(xi_elev.^2-1).*eta_elev.*(eta_elev-1)
                       (xi_elev.^2-1).*(eta_elev.^2-1)],[2,1,3]);
                   
         N_diff_rows_points_pages=[(eta_elev.*(2*xi_elev + 1).*(eta_elev - 1))/4, -((eta_elev.^2 - 1).*(2*xi_elev + 1))/2, (eta_elev.*(2*xi_elev + 1).*(eta_elev + 1))/4,           -eta_elev.*xi_elev.*(eta_elev + 1), (eta_elev.*(2*xi_elev - 1).*(eta_elev + 1))/4, -((eta_elev.^2 - 1).*(2*xi_elev - 1))/2, (eta_elev.*(2*xi_elev - 1).*(eta_elev - 1))/4,           -eta_elev.*xi_elev.*(eta_elev - 1), 2*xi_elev.*(eta_elev.^2 - 1)
                                   (xi_elev.*(2*eta_elev - 1).*(xi_elev + 1))/4,            -eta_elev.*xi_elev.*(xi_elev + 1),  (xi_elev.*(2*eta_elev + 1).*(xi_elev + 1))/4, -((2*eta_elev + 1).*(xi_elev.^2 - 1))/2,  (xi_elev.*(2*eta_elev + 1).*(xi_elev - 1))/4,            -eta_elev.*xi_elev.*(xi_elev - 1),  (xi_elev.*(2*eta_elev - 1).*(xi_elev - 1))/4, -((2*eta_elev - 1).*(xi_elev.^2 - 1))/2, 2*eta_elev.*(xi_elev.^2 - 1)];
        
    case 4   %T3
        lamda_pages=1-xi_elev-eta_elev;
        N_row_points_pages=[xi_elev   eta_elev   lamda_pages];
        N_diff_rows_points_pages=repmat([1, 0, -1
                                  0, 1, -1],1,1,N_points);
    
    case 5   %Tri6
        lamda_pages=1-xi_elev-eta_elev;
        N_row_points_pages=[-xi_elev.*(1-2*xi_elev) ...
                            4*xi_elev.*eta_elev ...
                            -eta_elev.*(1-2*eta_elev) ...
                            4*eta_elev.*lamda_pages ...
                            -lamda_pages.*(1-2*lamda_pages) ...
                            4*xi_elev.*lamda_pages];

        N_diff_rows_points_pages=[-1+4*xi_elev, 4*eta_elev,        zeros(1,1,N_points),        -4*eta_elev, 1-4*lamda_pages, 4*(lamda_pages-xi_elev)
                                  zeros(1,1,N_points),  4*xi_elev, -1+4*eta_elev, 4*(lamda_pages-eta_elev), 1-4*lamda_pages,       -4*xi_elev];

    case 11    %Hex8
        N_row_points_pages=[(1-zeta_elev).*(1-eta_elev).*(1-xi_elev) ...
                            (1+zeta_elev).*(1-eta_elev).*(1-xi_elev) ...
                            (1+zeta_elev).*(1+eta_elev).*(1-xi_elev) ...
                            (1-zeta_elev).*(1+eta_elev).*(1-xi_elev) ...
                            (1-zeta_elev).*(1-eta_elev).*(1+xi_elev) ...
                            (1+zeta_elev).*(1-eta_elev).*(1+xi_elev) ...
                            (1+zeta_elev).*(1+eta_elev).*(1+xi_elev) ...
                            (1-zeta_elev).*(1+eta_elev).*(1+xi_elev)]/8;
        
        N_diff_rows_points_pages=permute([  -(eta_elev - 1).*(zeta_elev - 1), -(xi_elev - 1).*(zeta_elev - 1), -(eta_elev - 1).*(xi_elev - 1)
                                             (eta_elev + 1).*(zeta_elev - 1),  (xi_elev - 1).*(zeta_elev - 1),  (eta_elev + 1).*(xi_elev - 1)
                                            -(eta_elev + 1).*(zeta_elev + 1), -(xi_elev - 1).*(zeta_elev + 1), -(eta_elev + 1).*(xi_elev - 1)
                                             (eta_elev - 1).*(zeta_elev + 1),  (xi_elev - 1).*(zeta_elev + 1),  (eta_elev - 1).*(xi_elev - 1)
                                             (eta_elev - 1).*(zeta_elev - 1),  (xi_elev + 1).*(zeta_elev - 1),  (eta_elev - 1).*(xi_elev + 1)
                                            -(eta_elev + 1).*(zeta_elev - 1), -(xi_elev + 1).*(zeta_elev - 1), -(eta_elev + 1).*(xi_elev + 1)
                                             (eta_elev + 1).*(zeta_elev + 1),  (xi_elev + 1).*(zeta_elev + 1),  (eta_elev + 1).*(xi_elev + 1)
                                            -(eta_elev - 1).*(zeta_elev + 1), -(xi_elev + 1).*(zeta_elev + 1), -(eta_elev - 1).*(xi_elev + 1)]/8,[2,1,3]);

    case 22   %Hex20
        N_row_points_pages=[(xi_elev/8 - 1/8).*(eta_elev - 1).*(zeta_elev - 1).*(eta_elev + xi_elev + zeta_elev + 2) ...
                           -(xi_elev/8 + 1/8).*(eta_elev - 1).*(zeta_elev - 1).*(eta_elev - xi_elev + zeta_elev + 2) ...
                           -(xi_elev/8 + 1/8).*(eta_elev + 1).*(zeta_elev - 1).*(eta_elev + xi_elev - zeta_elev - 2) ...
                           -(xi_elev/8 - 1/8).*(eta_elev + 1).*(zeta_elev - 1).*(xi_elev - eta_elev + zeta_elev + 2) ...
                           -(xi_elev/8 - 1/8).*(eta_elev - 1).*(zeta_elev + 1).*(eta_elev + xi_elev - zeta_elev + 2) ...
                            (xi_elev/8 + 1/8).*(eta_elev - 1).*(zeta_elev + 1).*(eta_elev - xi_elev - zeta_elev + 2) ...
                            (xi_elev/8 + 1/8).*(eta_elev + 1).*(zeta_elev + 1).*(eta_elev + xi_elev + zeta_elev - 2) ...
                           -(xi_elev/8 - 1/8).*(eta_elev + 1).*(zeta_elev + 1).*(eta_elev - xi_elev + zeta_elev - 2) ...
                           -(xi_elev.^2/4 - 1/4).*(eta_elev - 1).*(zeta_elev - 1) ...
                            (eta_elev.^2 - 1).*(xi_elev/4 + 1/4).*(zeta_elev - 1) ...
                            (xi_elev.^2/4 - 1/4).*(eta_elev + 1).*(zeta_elev - 1) ...
                           -(eta_elev.^2 - 1).*(xi_elev/4 - 1/4).*(zeta_elev - 1) ...
                           -(xi_elev/4 - 1/4).*(zeta_elev.^2 - 1).*(eta_elev - 1) ...
                            (xi_elev/4 + 1/4).*(zeta_elev.^2 - 1).*(eta_elev - 1) ...
                           -(xi_elev/4 + 1/4).*(zeta_elev.^2 - 1).*(eta_elev + 1) ...
                            (xi_elev/4 - 1/4).*(zeta_elev.^2 - 1).*(eta_elev + 1) ...
                            (xi_elev.^2/4 - 1/4).*(eta_elev - 1).*(zeta_elev + 1) ...
                           -(eta_elev.^2 - 1).*(xi_elev/4 + 1/4).*(zeta_elev + 1) ...
                           -(xi_elev.^2/4 - 1/4).*(eta_elev + 1).*(zeta_elev + 1) ...
                            (eta_elev.^2 - 1).*(xi_elev/4 - 1/4).*(zeta_elev + 1)];
        
        N_diff_rows_points_pages=permute([ ((eta_elev - 1).*(zeta_elev - 1).*(eta_elev + 2*xi_elev + zeta_elev + 1))/8,  ((xi_elev - 1).*(zeta_elev - 1).*(2*eta_elev + xi_elev + zeta_elev + 1))/8,  ((eta_elev - 1).*(xi_elev - 1).*(eta_elev + xi_elev + 2*zeta_elev + 1))/8
                                 -((eta_elev - 1).*(zeta_elev - 1).*(eta_elev - 2*xi_elev + zeta_elev + 1))/8, -((xi_elev + 1).*(zeta_elev - 1).*(2*eta_elev - xi_elev + zeta_elev + 1))/8, -((eta_elev - 1).*(xi_elev + 1).*(eta_elev - xi_elev + 2*zeta_elev + 1))/8
                                 -((eta_elev + 1).*(zeta_elev - 1).*(eta_elev + 2*xi_elev - zeta_elev - 1))/8, -((xi_elev + 1).*(zeta_elev - 1).*(2*eta_elev + xi_elev - zeta_elev - 1))/8, -((eta_elev + 1).*(xi_elev + 1).*(eta_elev + xi_elev - 2*zeta_elev - 1))/8
                                 -((eta_elev + 1).*(zeta_elev - 1).*(2*xi_elev - eta_elev + zeta_elev + 1))/8, -((xi_elev - 1).*(zeta_elev - 1).*(xi_elev - 2*eta_elev + zeta_elev + 1))/8, -((eta_elev + 1).*(xi_elev - 1).*(xi_elev - eta_elev + 2*zeta_elev + 1))/8
                                 -((eta_elev - 1).*(zeta_elev + 1).*(eta_elev + 2*xi_elev - zeta_elev + 1))/8, -((xi_elev - 1).*(zeta_elev + 1).*(2*eta_elev + xi_elev - zeta_elev + 1))/8, -((eta_elev - 1).*(xi_elev - 1).*(eta_elev + xi_elev - 2*zeta_elev + 1))/8
                                  ((eta_elev - 1).*(zeta_elev + 1).*(eta_elev - 2*xi_elev - zeta_elev + 1))/8,  ((xi_elev + 1).*(zeta_elev + 1).*(2*eta_elev - xi_elev - zeta_elev + 1))/8,  ((eta_elev - 1).*(xi_elev + 1).*(eta_elev - xi_elev - 2*zeta_elev + 1))/8
                                  ((eta_elev + 1).*(zeta_elev + 1).*(eta_elev + 2*xi_elev + zeta_elev - 1))/8,  ((xi_elev + 1).*(zeta_elev + 1).*(2*eta_elev + xi_elev + zeta_elev - 1))/8,  ((eta_elev + 1).*(xi_elev + 1).*(eta_elev + xi_elev + 2*zeta_elev - 1))/8
                                 -((eta_elev + 1).*(zeta_elev + 1).*(eta_elev - 2*xi_elev + zeta_elev - 1))/8, -((xi_elev - 1).*(zeta_elev + 1).*(2*eta_elev - xi_elev + zeta_elev - 1))/8, -((eta_elev + 1).*(xi_elev - 1).*(eta_elev - xi_elev + 2*zeta_elev - 1))/8
                                                      -(xi_elev.*(eta_elev - 1).*(zeta_elev - 1))/2,                       -(xi_elev.^2/4 - 1/4).*(zeta_elev - 1),                       -(xi_elev.^2/4 - 1/4).*(eta_elev - 1)
                                                        ((eta_elev.^2 - 1).*(zeta_elev - 1))/4,                    2*eta_elev.*(xi_elev/4 + 1/4).*(zeta_elev - 1),                        (eta_elev.^2 - 1).*(xi_elev/4 + 1/4)
                                                       (xi_elev.*(eta_elev + 1).*(zeta_elev - 1))/2,                        (xi_elev.^2/4 - 1/4).*(zeta_elev - 1),                        (xi_elev.^2/4 - 1/4).*(eta_elev + 1)
                                                       -((eta_elev.^2 - 1).*(zeta_elev - 1))/4,                   -2*eta_elev.*(xi_elev/4 - 1/4).*(zeta_elev - 1),                       -(eta_elev.^2 - 1).*(xi_elev/4 - 1/4)
                                                       -((zeta_elev.^2 - 1).*(eta_elev - 1))/4,                       -(xi_elev/4 - 1/4).*(zeta_elev.^2 - 1),                  -2*zeta_elev.*(xi_elev/4 - 1/4).*(eta_elev - 1)
                                                        ((zeta_elev.^2 - 1).*(eta_elev - 1))/4,                        (xi_elev/4 + 1/4).*(zeta_elev.^2 - 1),                   2*zeta_elev.*(xi_elev/4 + 1/4).*(eta_elev - 1)
                                                       -((zeta_elev.^2 - 1).*(eta_elev + 1))/4,                       -(xi_elev/4 + 1/4).*(zeta_elev.^2 - 1),                  -2*zeta_elev.*(xi_elev/4 + 1/4).*(eta_elev + 1)
                                                        ((zeta_elev.^2 - 1).*(eta_elev + 1))/4,                        (xi_elev/4 - 1/4).*(zeta_elev.^2 - 1),                   2*zeta_elev.*(xi_elev/4 - 1/4).*(eta_elev + 1)
                                                       (xi_elev.*(eta_elev - 1).*(zeta_elev + 1))/2,                        (xi_elev.^2/4 - 1/4).*(zeta_elev + 1),                        (xi_elev.^2/4 - 1/4).*(eta_elev - 1)
                                                       -((eta_elev.^2 - 1).*(zeta_elev + 1))/4,                   -2*eta_elev.*(xi_elev/4 + 1/4).*(zeta_elev + 1),                       -(eta_elev.^2 - 1).*(xi_elev/4 + 1/4)
                                                      -(xi_elev.*(eta_elev + 1).*(zeta_elev + 1))/2,                       -(xi_elev.^2/4 - 1/4).*(zeta_elev + 1),                       -(xi_elev.^2/4 - 1/4).*(eta_elev + 1)
                                                        ((eta_elev.^2 - 1).*(zeta_elev + 1))/4,                    2*eta_elev.*(xi_elev/4 - 1/4).*(zeta_elev + 1),                        (eta_elev.^2 - 1).*(xi_elev/4 - 1/4)],[2,1,3]);
        
    case 33   %Tet4
       N_row_points_pages=[1-xi_elev-eta_elev-zeta_elev   xi_elev   eta_elev   zeta_elev];
       N_diff_rows_points_pages=repmat([-1 1 0 0
                                        -1 0 1 0
                                        -1 0 0 1],1,1,N_points);
%        N_row_points_pages=circshift(N_row_points_pages,-1,2);
%        N_diff_rows_points_pages=circshift(N_diff_rows_points_pages,-1,2);
    
    case 44   %Tet10
        N_row_points_pages=[-(1-xi_elev-eta_elev-zeta_elev).*(1-2*(1-xi_elev-eta_elev-zeta_elev)) ...
                            -xi_elev.*(1-2*xi_elev) ...
                            -eta_elev.*(1-2*eta_elev) ...
                            -zeta_elev.*(1-2*zeta_elev) ...
                            4*xi_elev.*(1-xi_elev-eta_elev-zeta_elev) ...
                            4*xi_elev.*eta_elev ...
                            4*eta_elev.*(1-xi_elev-eta_elev-zeta_elev) ...
                            4*zeta_elev.*(1-xi_elev-eta_elev-zeta_elev) ...
                            4*xi_elev.*zeta_elev ...
                            4*eta_elev.*zeta_elev];
        
        N_diff_rows_points_pages=permute([4*eta_elev + 4*xi_elev+4*zeta_elev-3, 4*eta_elev+4*xi_elev+4*zeta_elev-3, 4*eta_elev+4*xi_elev+4*zeta_elev-3
                                                4*xi_elev-1,                   zeros(1,1,N_points),                   zeros(1,1,N_points)
                                                     zeros(1,1,N_points),             4*eta_elev-1,                   zeros(1,1,N_points)
                                                     zeros(1,1,N_points),                   zeros(1,1,N_points),            4*zeta_elev-1
                                   4-8*xi_elev-4*zeta_elev-4*eta_elev,               -4*xi_elev,               -4*xi_elev
                                                 4*eta_elev,                4*xi_elev,                   zeros(1,1,N_points)
                                                -4*eta_elev, 4-4*xi_elev-4*zeta_elev-8*eta_elev,              -4*eta_elev
                                               -4*zeta_elev,             -4*zeta_elev, 4-4*xi_elev-8*zeta_elev-4*eta_elev
                                                4*zeta_elev,                   zeros(1,1,N_points),                4*xi_elev
                                                     zeros(1,1,N_points),              4*zeta_elev,               4*eta_elev],[2,1,3]);
        
end
