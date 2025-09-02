function [gaussPoints_cvec,gaussWeights_cvec] = gaussLegendreQuad(numGaussPoints)

    if numGaussPoints == 1
    
        gaussPoints_cvec  = 0;
        gaussWeights_cvec = 2;
    
    elseif numGaussPoints == 2
    
        gaussPoints_cvec  = [1/sqrt(3);-1/sqrt(3)];
        gaussWeights_cvec = [1;1];
    
    elseif numGaussPoints == 3
    
        gaussPoints_cvec  = [0;sqrt(3/5);-sqrt(3/5)];
        gaussWeights_cvec = [8/9;5/9;5/9];
    
    elseif numGaussPoints == 4
    
        gaussPoints_cvec  = [sqrt(3/7 - 2/7 * sqrt(6/5));-sqrt(3/7 - 2/7 * sqrt(6/5));sqrt(3/7 + 2/7 * sqrt(6/5));-sqrt(3/7 + 2/7 * sqrt(6/5))];
        gaussWeights_cvec = [(18+sqrt(30))/36;(18+sqrt(30))/36;(18-sqrt(30))/36;(18-sqrt(30))/36];
    
    elseif numGaussPoints == 5
    
        gaussPoints_cvec  = [0;1/3 * sqrt(5 - 2 * sqrt(10/7));-1/3 * sqrt(5 - 2 * sqrt(10/7));1/3 * sqrt(5 + 2 * sqrt(10/7));-1/3 * sqrt(5 + 2 * sqrt(10/7))];
        gaussWeights_cvec = [128/225; (322 + 13*sqrt(70))/900;(322 + 13*sqrt(70))/900;(322 - 13*sqrt(70))/900;(322 - 13*sqrt(70))/900];
    
    end

end