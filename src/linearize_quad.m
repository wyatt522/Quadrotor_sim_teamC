function [A, B] = linearize_quad(quad, position)
    syms x y z xdot ydot zdot phi theta psiVar omega1 omega2 omega3 u1 u2 u3 u4 ...
        n1 n2 n3 r1 r2 r3
    l = quad.l; m = quad.m; I1 = quad.I(1, 1); I2 = quad.I(2, 2); I3 = quad.I(3,3); g = quad.g; sigma = quad.sigma; % all in quad object
    
    Rz = [cos(psiVar) -sin(psiVar) 0; sin(psiVar) cos(psiVar) 0; 0 0 1]; % R3
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]; %R2
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)]; %R1
    
    R_CE = Rz * Ry * Rx; 
    R_EC = transpose(Rx) * transpose(Ry) * transpose(Rz);
    
    T = [1 0 -sin(theta); 0 cos(phi) sin(phi)*cos(theta); 0 -sin(phi) cos(phi)*cos(theta)]; % omega = T * alphadot
    
    pos = [x; y; z];
    alpha = [phi; theta; psiVar];
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
    
    z0 = [position, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    u0 = m*g/4 * ones(4,1);
    Ja = jacobian(qdot, q);
    Ja_eval = subs(Ja, [q; u; r], [z0; u0; zeros(size(r))]);
    A = eval(Ja_eval);
    Jb = jacobian(qdot, u);
    Jb_eval = subs(Jb, [q; u; r], [z0; u0; zeros(size(r))]);
    B = eval(Jb_eval);
    if rank(ctrb(A,B)) ~= 12
        error("Controllability maxtrix is not full rank.")
    end
end 
function zdot = f(z, u)
    Rz = [cos(q(6)) -sin(q(6)) 0; sin(q(6)) cos(q(6)) 0; 0 0 1]; % R3
    Ry = [cos(q(5)) 0 sin(q(5)); 0 1 0; -sin(q(5)) 0 cos(q(5))]; %R2
    Rx = [1 0 0; 0 cos(q(4)) -sin(q(4)); 0 sin(q(4)) cos(q(4))]; %R1
    R_CE = Rz * Ry * Rx;

    T = [1 0 -sin(q(5)); 0 cos(q(4)) sin(q(4))*cos(q(5)); 0 -sin(q(4)) cos(q(4))*cos(q(5))]; % omega = T * alphadot
    I = diag([I1, I2, I3]);
    u = [u1; u2; u3; u4];
    n = [n1; n2; n3];
    r = [r1; r2; r3];
    omega = q(10:12,1);

    pos_dot = q(7:9,1);
    alpha_dot = T\omega;
    vel_dot = (-g * [0; 0; 1]) + (1/m * R_CE * (u(1) + u(2) + u(3) + u(4)) * [0; 0; 1]) + (1/m * R_CE * r);
    omega_dot = I\(((u(2) - u(4)) * l * [1;0;0]) + ((u(3) - u(1)) * l * [0;1;0]) + ((u(1) - u(2) + u(3) - u(4)) * sigma * [0;0;1]) + n - cross(omega,I*omega));
end