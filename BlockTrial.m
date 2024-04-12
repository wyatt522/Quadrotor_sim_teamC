%% Initializations
clc; clear; close all;
m = 10; % Mass of the block [kg]
b = 2; % Damping coefficient [Ns/m]
% k = 4; % Spring stiffness [N/m]
k = 50;
t0 = 0; % Simulation start time
tf = 15; % Simulation finish time
y0 = 0; % y0 = y(t0)
v0 = 1; % v0 = dy/dt(t0)
 % Reference (desired) input
%% System Dynamics
% Defining control as an anonymous function increases code modularity.
ts = 3;
r = 5;
rdot = @(x) 0-x;
phi = @(x) k*(rdot(x(1)) - x(2));
f = @(t, x, u)[x(2); u/m - b/m * x(2)];
%% Solving for The Solution
x0 = [y0; v0];
t = linspace(t0, tf, 100)';
% We can use an anonymous function to reduce the input numbers of the f
% function by fixing u as phi(x)
[~,x] = ode45(@(t,x)f(t,x,phi(x)), t, x0);
% for i=1:100
%     rdot_disp(i) = rdot(x(:,1));
% end
%% Display the results
ax = axes('NextPlot','add','TickLabelInterpreter','LaTeX','FontSize',20,...
'XLim',[t(1), t(end)],...
'Xgrid','on','Ygrid','on','Box','on');
line(t,x(:,1),'color',[0, 0.447, 0.741],'LineWidth',2,...
'Parent',ax,'DisplayName','$y(t)$');
line(t,x(:,2),'color',[0.85, 0.325, 0.098],'LineWidth',2,...
'Parent',ax,'DisplayName','$\dot{y}(t)$');
% line(t([1,end]),[r r],'color',[0 0 0],'LineWidth',2,...
% 'LineStyle','--','Parent',ax,'DisplayName','$r$');
% line(t,rdot_disp,'color',[0 0 0],'LineWidth',2,...
% 'LineStyle','--','Parent',ax,'DisplayName','$rdot$');
xlabel(ax,'$t$ [s]','Interpreter','LaTeX','FontSize',20);
legend(ax,'Interpreter','LaTeX','FontSize',20,...
'Position', [0.786, 0.666, 0.098, 0.135]);

% 'Ylim',[-0.1 r + 0.1]