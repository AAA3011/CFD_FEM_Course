% readGeoAndBC - Read geometry and boundary condition node sets.
%
%
% Reads element and node lists and boundary node sets from text files used
% by the cylinder test case and returns connectivity and boolean node vectors.
%
% Outputs:
%   connectivityMatrix_mat : Element connectivity matrix (nElements x nodesPerElement)
%   xCoord_vec : Nodal x-coordinates of the mesh.
%   yCoord_vec : Nodal y-coordinates of the mesh.
%   top_vec : Node indices on the top boundary.
%   outlet_vec : Node indices on the outlet boundary.
%   inlet_vec : Node indices on the inlet boundary.
%   cylinderWall_vec : Node indices on the cylinder wall boundary.
%   bottom_vec : Node indices on the bottom boundary.
function [connectivityMatrix_mat, xCoord_vec, yCoord_vec, top_vec, outlet_vec, inlet_vec, cylinderWall_vec, bottom_vec] = readGeoAndBC()

Elements      = readtable('Elements11.txt',      'VariableNamingRule', 'preserve');
Nodes         = readtable('Nodes11.txt',         'VariableNamingRule', 'preserve');
top           = readtable('top11.txt',           'VariableNamingRule', 'preserve');
bottom        = readtable('bottom11.txt',        'VariableNamingRule', 'preserve');
inlet         = readtable('inlet11.txt',         'VariableNamingRule', 'preserve');
outlet        = readtable('outlet11.txt',        'VariableNamingRule', 'preserve');
cylinderWall  = readtable('cylinderWall11.txt',  'VariableNamingRule', 'preserve');

%% Geometry
connectivityMatrix_mat = Elements{:, 3:6};
xCoord_vec             = Nodes{:, 2}';
yCoord_vec             = Nodes{:, 3}';

%% Boundary Condition
top_vec          = top{:, 1}';
outlet_vec       = outlet{:, 1}';
inlet_vec        = inlet{:, 1}';
cylinderWall_vec = cylinderWall{:, 1}';
bottom_vec       = bottom{:, 1}';

end
