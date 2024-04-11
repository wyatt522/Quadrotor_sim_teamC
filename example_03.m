clc; clear; close all

addpath('./src');

% QUADROTOR

g = 9.81;  % The gravitational acceleration [m/s^2]
l = 0.2;  % Distance from the center of mass to each rotor [m]
m = 0.5;  % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48];  % Mass moment of inertia [kg m^2]
mu = 3.0;  % Maximum thrust of each rotor [N]
sigma = 0.01;  % The proportionality constant relating thrust to torque [m]

quad = quadrotor(g, l, m, diag(I), mu, sigma);

% INTRUDER
path = @(t) [cos(t); sin(t); 2];
dist = struct("r", @(t)0.1*[sin(t); sin(2*t); sin(4*t)],...
    "n", @(t)[0.1; 0.01; 0.1]);

intruder = uav(path, dist);

% CONTROLLER
ctrl = controller2([-0.5 0.5 3], 0, [3 1], quad);

% SIMULATION

sim = simulator(quad, ctrl, intruder);
sim.simtime = [0 5];
sim.timestep = 0.01;

z0 = zeros(12,1);

[t,z,u,d,y] = sim.simulate(z0);


%PLOT
u1i = 1:4:length(u);
u3i = 3:4:length(u);



figure;

% Plot for X Position
subplot(2, 1, 1); 
plot(t, z(:, 1), 'bo-', 'DisplayName', 'X Position');
hold on;
plot(t, z(:, 7), 'ro-', 'DisplayName', 'Velocity');
grid on;
xlabel('Time (s)');
ylabel('X Position (m)');
title('X Position over Time');
legend('Location', 'best');

% Plot for motors
subplot(2, 1, 2); 
grid on; % Turn on grid before plotting the lines
xlabel('Time (s)');
ylabel('Rotors');
title('Rotors over Time');
plot(t(1:500, :), ctrl.uvec(1,:), 'bo-', 'DisplayName', 'u1');
hold on;
plot(t(1:500, :), ctrl.uvec(3,:), 'ro-', 'DisplayName', 'u3');
legend('Location', 'best');
% ANIMATION
sim.animate(t, z, y);


