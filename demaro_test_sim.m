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
%path = @(t) [(-5 + t); (5 - t); 5];
%path = @(t) [-5 + 1.5*t; -3+t; t];
path = @(t) [4*cos(t); 4*sin(t); 4 + sin(t)];
% path = @(t) [-5+t; cos(t); 5];
dist = struct("r", @(t,z)0.1*[sin(t); sin(2*t); sin(4*t)],...
    "n", @(t,z) 0.1*[0.1; 0.01; 0.1]);

intruder = uav(path, dist);

ctrl = SAC(quad);

% SIMULATION

sim = simulator(quad, ctrl, intruder);
sim.simtime = [0 40];
sim.timestep = 0.01;
sim.epsilon = 0.1;

z0 = zeros(12,1);

[t,z,u,d,y] = sim.simulate(z0);

% ANIMATION
sim.animate(t, z, y);
% quad.plot(t,z)