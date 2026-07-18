function [x_eta1, x_eta2, y_eta1, y_eta2] = computeDerivatives(domainXStart, domainXEnd, domainYStart, domainYEnd, x, y, xCoord_mat, yCoord_mat, i, j, eta1_mat, eta2_mat)
    disp(i)
    disp(j)

    % Get points for eta1 direction (varying i index, fixed j)
    [x1_eta1, x2_eta1, eta11, eta12] = getEta1Points(domainXStart, domainXEnd, domainYStart, domainYEnd, x, y, xCoord_mat, yCoord_mat, i, j, eta1_mat, eta2_mat);
    
    % Get points for eta2 direction (varying j index, fixed i)  
    [y1_eta2, y2_eta2, eta21, eta22] = getEta2Points(domainXStart, domainXEnd, domainYStart, domainYEnd, x, y, xCoord_mat, yCoord_mat, i, j, eta1_mat, eta2_mat);
    
    % Get y values for eta1 direction and x values for eta2 direction
    y1_eta1 = yCoord_mat(i+1, j);
    y2_eta1 = yCoord_mat(i-1, j);
    x1_eta2 = xCoord_mat(i, j+1);
    x2_eta2 = xCoord_mat(i, j-1);

    % Compute derivatives
    x_eta1 = (x1_eta1 - x2_eta1) / (eta11 - eta12);
    y_eta1 = (y1_eta1 - y2_eta1) / (eta11 - eta12);
    x_eta2 = (x1_eta2 - x2_eta2) / (eta21 - eta22);
    y_eta2 = (y1_eta2 - y2_eta2) / (eta21 - eta22);
end