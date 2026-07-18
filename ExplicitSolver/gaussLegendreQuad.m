% gaussLegendreQuad - Return GaussLegendre points and weights.
%
% FILE: gaussLegendreQuad.m
% DESCRIPTION:
% Returns 1D Gauss-Legendre quadrature points and weights for integration
% on reference interval [-1, 1]. Used to construct 2D quadrature rules via
% tensor product for QUAD4 elements.
%
% Inputs:
%   numGaussPoints (variable): Number of Gauss points used per dimension
% Outputs:
%   gaussPoints_col : Column vector of Gauss-Legendre quadrature points on [-1,1].
%   gaussWeights_col : Column vector of Gauss-Legendre quadrature weights.
function [gaussPoints_col,gaussWeights_col] = gaussLegendreQuad(numGaussPoints)

    if numGaussPoints == 1
    
        gaussPoints_col  = 0;
        gaussWeights_col = 2;
    
    elseif numGaussPoints == 2
    
        gaussPoints_col  = [1/sqrt(3);-1/sqrt(3)];
        gaussWeights_col = [1;1];
    
    elseif numGaussPoints == 3
    
        gaussPoints_col  = [0;sqrt(3/5);-sqrt(3/5)];
        gaussWeights_col = [8/9;5/9;5/9];
    
    elseif numGaussPoints == 4
    
        gaussPoints_col  = [sqrt(3/7 - 2/7 * sqrt(6/5));-sqrt(3/7 - 2/7 * sqrt(6/5));sqrt(3/7 + 2/7 * sqrt(6/5));-sqrt(3/7 + 2/7 * sqrt(6/5))];
        gaussWeights_col = [(18+sqrt(30))/36;(18+sqrt(30))/36;(18-sqrt(30))/36;(18-sqrt(30))/36];
    
    elseif numGaussPoints == 5
    
        gaussPoints_col  = [0;1/3 * sqrt(5 - 2 * sqrt(10/7));-1/3 * sqrt(5 - 2 * sqrt(10/7));1/3 * sqrt(5 + 2 * sqrt(10/7));-1/3 * sqrt(5 + 2 * sqrt(10/7))];
        gaussWeights_col = [128/225; (322 + 13*sqrt(70))/900;(322 + 13*sqrt(70))/900;(322 - 13*sqrt(70))/900;(322 - 13*sqrt(70))/900];
    
    end

end
