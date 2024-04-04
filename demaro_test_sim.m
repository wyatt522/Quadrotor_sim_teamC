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

% tspan = 0:5:0.1;
% y = intruder.location(false,tspan,[]);
% zdes = zeros(12,1); zdes(1:3) = y(1:3);

% CONTROLLER

syms x y z xdot ydot zdot phi theta psi omega1 omega2 omega3 u1 u2 u3 u4 ...
    n1 n2 n3 r1 r2 r3

l = 0.2; m = 0.5; I1 = 1.24; I2 = 1.24; I3 = 2.48; g = 9.8; sigma = 0.01;

Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; % R3
Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]; %R2
Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)]; %R1

R_CE = Rz * Ry * Rx;
R_EC = transpose(Rx) * transpose(Ry) * transpose(Rz);

T = [1 0 -sin(theta); 0 cos(phi) sin(phi)*cos(theta); 0 -sin(phi) cos(phi)*cos(theta)]; % omega = T * alphadot

pos = [x; y; z];
alpha = [phi; theta; psi];
vel = [xdot; ydot; zdot];
omega = [omega1; omega2; omega3];

q = [pos; alpha; vel; omega];

I = diag([I1, I2, I3]);
u = [u1; u2; u3; u4];
n = [n1; n2; n3];
r = [r1; r2; r3];

pos_dot = [xdot; ydot; zdot];
alpha_dot = T\omega;
vel_dot = (-g * [0; 0; 1]) + (1/m * R_CE * (u(1) + u(2) + u(3) + u(4)) * [0; 0; 1]) + (1/m * R_CE * r);
omega_dot = I\(((u(2) - u(4)) * l * [1;0;0]) + ((u(3) - u(1)) * l * [0;1;0]) + ((u(1) - u(2) + u(3) - u(4)) * sigma * [0;0;1]) + n - cross(omega,I*omega));

qdot = [pos_dot; alpha_dot; vel_dot; omega_dot];

u0 = m*g/4 * ones(4,1);
qdes = zeros(12,1); qdes(3) = 1;
Ja = jacobian(qdot, q);
Ja_eval = subs(Ja, [q; u; r], [qdes; u0; zeros(size(r))]);
A = eval(Ja_eval);
Jb = jacobian(qdot, u);
Jb_eval = subs(Jb, [q; u; r], [qdes; u0; zeros(size(r))]);
B = eval(Jb_eval);

if rank(ctrb(A,B)) < ndims(A)
    disp("Not controllable")
end

Q = diag([1,1,1,0.5,0.5,0.5,0,0,0,2,2,2]);
R = eye(4);

[K,S,CLP] = lqr(A,B,Q,R);

ctrl = SAC(K, quad);

% SIMULATION

sim = simulator(quad, ctrl, intruder);
sim.simtime = [0 20];
sim.timestep = 0.01;

z0 = zeros(12,1);

[t,z,u,d,y] = sim.simulate(z0);

% ANIMATION
sim.animate(t, z, y);
% quad.plot(t,z)