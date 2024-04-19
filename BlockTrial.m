%% Initializations
clc; clear; close all;
m = 10; % Mass of the block [kg]
b = 2; % Damping coefficient [Ns/m]
% k = 4; % Spring stiffness [N/m]
k = 25;
t0 = 0; % Simulation start time
tf = 50; % Simulation finish time
N = (tf-t0) * 10;
y0 = 0; % y0 = y(t0)
v0 = 1; % v0 = dy/dt(t0)
 % Reference (desired) input
%% System Dynamics
% Defining control as an anonymous function increases code modularity.

r = @(t) 5;
rdot = @(t) 0;

error = @(x,t) r(t) - x;
errordot = @(x,t) rdot(t) - x;

phi = @(x,t) (errordot(x(2),t)) + k*(error(x(1),t)); % ctrl input with edot
% phi = @(x,t) k*error(x(1),t); % ctrl input without edot
f = @(t, x, u)[x(2); u/m - b/m * x(2)];
%% Solving for The Solution
x0 = [y0; v0];
t = linspace(t0, tf, N)';
for i=1:N
    refVec(i) = r(t(i));
    refDotVec(i) = rdot(t(i));
end
[~,x] = ode45(@(t,x)f(t,x,phi(x,t)), t, x0);
for i=1:N
    eVec(i) = error(x(i,1),t(i));
    edotVec(i) = errordot(x(i,2),t(i));
end
%% Display the results
SysResp = figure;
yyaxis left;
plot(t,x(:,1),'b'); hold on;
plot(t, refVec,'b--');
ylabel('Position, x [m]');
yyaxis right;
plot(t,x(:,2),'r');
plot(t, refDotVec,'r--');
ylabel('Velocity, v [m/s]');
legend('Block Position','Block Velocity','Reference Position','Reference Velocity');
title('System Response'); xlabel('Time, t [s]');

ErrorDyn = figure;
plot(t, eVec,'b'); hold on;
plot(t, edotVec,'r')
xlabel('Time, t [s]');
legend('Error','Change in Error');
title("Error of System");
