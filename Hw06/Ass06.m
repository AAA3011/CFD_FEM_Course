clc;clearvars;close all;

%% Q2 Ass 06

x = linspace(0,1,100);
C = -3:0.2:3;

figure;
for i = 1:length(C)
    y = 2/3 *(x.^(3/2))+C(i);
    plot(x,y)
    hold on
    y = -2/3 * (x.^(3/2))+C(i);
    plot(x,y)
    hold on
    grid on
end
title(' The Euler–Tricomi characteristics','Interpreter','latex')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
set(gca,'YDir','reverse')
ylim([0 1])
xlim([0 1])


%% Q3 inviscid Burger's equation
% clc;clearvars;close all;

r1 = linspace(-2.5,0.5,11);
t1 = linspace(0,2,20);

figure;
for i = 1:length(r1)
    x = t1./(1+r1(i)^2) + r1(i);
    plot(x,t1,'DisplayName',['r = ' num2str(r1(i))])
    hold on
end
grid on
legend show 
xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
title('Curves for different r values','Interpreter','latex')

r = linspace(-10,10,300);
t = linspace(0,10,20);
figure;
for i = 1:length(t)

    u = 1./(1+r.^2);
    x = t(i)./(1+r.^2) + r;
    plot3(x, t(i)*ones(size(r)), u);
    hold on
    grid on
    set(gca,'YDir','reverse')

end
xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('u-axis','Interpreter','latex')
title('3D Stepping Plot of $$u $$','Interpreter','latex')

% surface plotting 
[R,T] = meshgrid(r,t);

W = 1./(1+R.^2);
X = T./(1+R.^2) + R;

figure;
surf(X, T, W)

shading interp              
set(gca,'YDir','reverse')   

xlabel('x-axis','Interpreter','latex')
ylabel('t-axis','Interpreter','latex')
zlabel('u-axis','Interpreter','latex')
title('3D Surface Plot of u','Interpreter','latex')


t_SW = 1.5396;

t_vec = [0 0.5*t_SW t_SW 1.5*t_SW];

for i = 1:length(t_vec)
    figure;
    u = 1./(1+r.^2);
    x = t_vec(i)./(1+r.^2) + r;
    plot(x,u);
    grid on
    title(['Inviscid Burgers equation at $t = ', [num2str(t_vec(i)/t_SW), '\times t_{sw}'], '$'], 'Interpreter', 'latex');
    xlabel('x','Interpreter','latex')
    ylabel('u','Interpreter','latex')
end