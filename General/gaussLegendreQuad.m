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
    elseif numGaussPoints == 9
    
        gaussPoints_col = [ ...
            -0.968160239507626;
            -0.836031107326636;
            -0.613371432700590;
            -0.324253423403809;
            0;
            0.324253423403809;
            0.613371432700590;
            0.836031107326636;
            0.968160239507626 ];
    
        gaussWeights_col = [ ...
            0.081274388361574;
            0.180648160694857;
            0.260610696402935;
            0.312347077040003;
            0.330239355001260;
            0.312347077040003;
            0.260610696402935;
            0.180648160694857;
            0.081274388361574 ];
    elseif numGaussPoints == 10
    
        gaussPoints_col = [
            -0.973906528517172;
            -0.865063366688985;
            -0.679409568299024;
            -0.433395394129247;
            -0.148874338981631;
            0.148874338981631;
            0.433395394129247;
            0.679409568299024;
            0.865063366688985;
            0.973906528517172 ];
    
        gaussWeights_col = [
            0.066671344308688;
            0.149451349150581;
            0.219086362515982;
            0.269266719309996;
            0.295524224714753;
            0.295524224714753;
            0.269266719309996;
            0.219086362515982;
            0.149451349150581;
            0.066671344308688 ];
            elseif numGaussPoints == 20
    
        gaussPoints_col = [
            -0.9931285991850949;
            -0.9639719272779138;
            -0.9122344282513260;
            -0.8391169718222188;
            -0.7463319064601508;
            -0.6360536807265150;
            -0.5108670019508271;
            -0.3737060887154196;
            -0.2277858511416451;
            -0.07652652113349733;
             0.07652652113349733;
             0.2277858511416451;
             0.3737060887154196;
             0.5108670019508271;
             0.6360536807265150;
             0.7463319064601508;
             0.8391169718222188;
             0.9122344282513260;
             0.9639719272779138;
             0.9931285991850949 ];
    
        gaussWeights_col = [
            0.01761400713915212;
            0.04060142980038694;
            0.06267204833410906;
            0.08327674157670475;
            0.1019301198172404;
            0.1181945319615184;
            0.1316886384491766;
            0.1420961093183821;
            0.1491729864726037;
            0.1527533871307258;
            0.1527533871307258;
            0.1491729864726037;
            0.1420961093183821;
            0.1316886384491766;
            0.1181945319615184;
            0.1019301198172404;
            0.08327674157670475;
            0.06267204833410906;
            0.04060142980038694;
            0.01761400713915212 ];
    end
end
