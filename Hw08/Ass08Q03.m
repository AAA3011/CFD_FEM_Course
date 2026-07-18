clc;clearvars;close all;

%% Solution Using Picard Method
range     = -100:0.001:100;
solutionp  = zeros(1,length(range));

for i = 1:length(range)

    xP_k      = range(i);
    diffP     = 1;
    tolP      = 1e-7;
    while(diffP > tolP)
    
        xP_kplus1 = (6-5*xP_k)/xP_k;
        diffP     = abs(xP_kplus1 - xP_k);
        xP_k      = xP_kplus1 ;
    
    end
    solutionp(i) = xP_k;
end

figure;
plot(range,solutionp)
title('variation of solution with initial guess Using picard Method','Interpreter','latex')
xlabel('initial guess','Interpreter','latex')
ylabel('solution','Interpreter','latex')
grid on

%% Solution Using Newton-Raphson Method

solutionN = zeros(1,length(range));

for j = 1:length(range)
    
    xN_k  = range(j);
    diffN = 1;
    tolN  = 1e-7;
    while(diffN > tolN)
    
        deltaX    = (6 - 5*xN_k - xN_k^2)/(2*xN_k + 5);
        xN_kPlus1 = xN_k + deltaX;
        diffN     = abs(xN_kPlus1 - xN_k); 
        xN_k      = xN_kPlus1;
    end
        solutionN(j) = xN_k;
end

figure;
plot(range,solutionN)
title('variation of solution with initial guess Using Newton Raphson Method','Interpreter','latex')
xlabel('initial guess','Interpreter','latex')
ylabel('solution','Interpreter','latex')
grid on
