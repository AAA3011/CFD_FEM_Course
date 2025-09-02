clc;clearvars;close all;

numElements_vec = [420 2625 5145 10500 15120 17745 20580];
meanVel_vec     = [1.7118 1.5790 1.5498 1.5335 1.5258 1.5229 1.5189];
FOS = 1.25;
% Grid convergence index
firstSolution  = meanVel_vec(end-2);
secondSolution = meanVel_vec(end-1);
thirdSolution  = meanVel_vec(end);

epslon_21 = (secondSolution - firstSolution)/secondSolution;
epslon_32 = (thirdSolution  - secondSolution)/thirdSolution;

gridRefinementRatio = mean([(numElements_vec(end)/numElements_vec(end-1))^(1/2),(numElements_vec(end-1)/numElements_vec(end-2))^(1/2)]);
orderOfConvergence = round(log((thirdSolution - secondSolution) / (secondSolution - firstSolution)) / log(gridRefinementRatio));

r_p = gridRefinementRatio^orderOfConvergence;

richardsonSolution = secondSolution + (secondSolution - firstSolution)/(r_p-1);
GCI_21 = FOS*abs(epslon_21)/(r_p-1)*100;
GCI_32 = FOS*abs(epslon_32)/(r_p-1)*100;

ratio = GCI_32/GCI_21;

asymptoticRange = ratio / r_p;

figure;
loglog(numElements_vec,meanVel_vec)
xlabel('number of Elements $$\times 10^4$','Interpreter','latex')
ylabel('mean velocity','Interpreter','latex')
title('Grid convergence study','Interpreter','latex')
grid on

hold on;
h = zeros(5, 1);
h(1) = plot(NaN, NaN, 'Visible', 'off');
h(2) = plot(NaN, NaN, 'Visible', 'off');
h(3) = plot(NaN, NaN, 'Visible', 'off');
h(4) = plot(NaN, NaN, 'Visible', 'off');
h(5) = plot(NaN, NaN, 'Visible', 'off');

legend_text = sprintf([ ...
    'Order (p): %.1f (Higher=Better)\n' ...
    'GCI_{21}: %.2f%% (Lower=Better)\n' ...
    'GCI_{32}: %.2f%% (Lower=Better)\n' ...
    'Asymp. Range: %.2f (Near 1=Good)\n' ...
    'Ref. Ratio: %.2f (>1.3=Good)'], ...
    orderOfConvergence, GCI_21, GCI_32, asymptoticRange, gridRefinementRatio);

legend(legend_text, 'Location', 'best', 'FontSize', 9, 'Interpreter', 'latex');
