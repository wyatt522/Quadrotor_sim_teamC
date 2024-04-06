classdef lineartrial3 < handle
    properties(Access = public)
        zdes(12,1) double;
        K(4,12) double;
        u0(4,1) double;
        y_prev(3,1) double;
        y_prev_init(3,1) double;
        A(12,12) double;
        B(12,4) double;
    end


    methods(Access = public)
        function obj = lineartrial3(quadrotor)

            obj.zdes = zeros(12,1);
            obj.u0 = repmat(quadrotor.m*quadrotor.g/4, [4,1]);
            [obj.A, obj.B] = linearize(zeros(12,1), obj.u0, quadrotor);

            Q = diag([10 10 10 2.5 2.5 2.5 5 5 5 5 5 5]);
            R = 2.5*eye(4);
            obj.K = lqr(obj.A,obj.B,Q,R);

            obj.y_prev = [NaN; NaN; NaN];
            obj.y_prev_init = [NaN; NaN; NaN];
        end

        function projected_dist = calc_ref(obj, z, y)
            if obj.y_prev ~= obj.y_prev_init
                delta_y = y - obj.y_prev;
                normalized_direction = normalize(delta_y);
                dist_quad_uav = norm(z(1:3) - y(1:3));
                projected_dist = y(1:3) + dist_quad_uav * normalized_direction;
            else
                projected_dist = y(1:3);
            end
            obj.y_prev = y;
        end

        function u = output(obj, ~, z, y)
            ref = obj.calc_ref(z, y);
            % disp(ref);
            obj.zdes(1:3) = ref;
            e = obj.zdes - z;
            v = -obj.K * e;
            u = obj.u0 - v;
            % disp(u);
        end
    end
end

function [A, B] = linearize(q_eq, u_eq, quadrotor)
    syms x y z xdot ydot zdot alpha_phi alpha_theta alpha_psi omega1 omega2 omega3 u1 u2 u3 u4 ...
        n1 n2 n3 r1 r2 r3
    
    Rz = [cos(alpha_psi) -sin(alpha_psi) 0; sin(alpha_psi) cos(alpha_psi) 0; 0 0 1]; % R3
    Ry = [cos(alpha_theta) 0 sin(alpha_theta); 0 1 0; -sin(alpha_theta) 0 cos(alpha_theta)]; %R2
    Rx = [1 0 0; 0 cos(alpha_phi) -sin(alpha_phi); 0 sin(alpha_phi) cos(alpha_phi)]; %R1
    
    R_CE = Rz * Ry * Rx;
    R_EC = transpose(Rx) * transpose(Ry) * transpose(Rz);
    
    T = [1 0 -sin(alpha_theta); 0 cos(alpha_phi) sin(alpha_phi)*cos(alpha_theta); 0 -sin(alpha_phi) cos(alpha_phi)*cos(alpha_theta)]; % omega = T * alphadot
    
    pos = [x; y; z];
    alpha = [alpha_phi; alpha_theta; alpha_psi];
    vel = [xdot; ydot; zdot];
    omega = [omega1; omega2; omega3];
    
    q = [pos; alpha; vel; omega];
    
    u = [u1; u2; u3; u4];
    n = [n1; n2; n3];
    r = [r1; r2; r3];
    
    pos_dot = [xdot; ydot; zdot];
    alpha_dot = T\omega;
    vel_dot = (-quadrotor.g * [0; 0; 1]) + (1/quadrotor.m * R_CE * (u(1) + u(2) + u(3) + u(4)) * [0; 0; 1]) + (1/quadrotor.m * R_CE * r);
    omega_dot = quadrotor.I\(((u(2) - u(4)) * quadrotor.l * [1;0;0]) + ((u(3) - u(1)) * quadrotor.l * [0;1;0]) + ((u(1) - u(2) + u(3) - u(4)) * quadrotor.sigma * [0;0;1]) + n - cross(omega,quadrotor.I*omega));
    
    qdot = [pos_dot; alpha_dot; vel_dot; omega_dot];
    
    Ja = jacobian(qdot, q);
    Ja_eval = subs(Ja, [q; u; r], [q_eq; u_eq; zeros(size(r))]);
    A = eval(Ja_eval);
    Jb = jacobian(qdot, u);
    Jb_eval = subs(Jb, [q; u; r], [q_eq; u_eq; zeros(size(r))]);
    B = eval(Jb_eval);
end