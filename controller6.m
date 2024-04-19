%% Controller 6
% Attempting to properly cast reference positions and velocities
classdef controller6 < handle
    properties(Access = public)
        zdes(12,1) double;
        K_to_capture(4,12) double;
        K_to_return(4,12) double;
        u0(4,1) double;
        y_prev(3,1) double;
        y_prev_init(3,1) double;
        A(12,12) double;
        B(12,4) double;
        epsilon(1,1) double;
        uvec;
        projdist;
        prev_ref(3,1) double;
        refvelvec;
        ref;
        time_step;
        uav_vel;
        ref_vel;
        compVel;
        last_e;
        prev_coeffs;
    end


    methods(Access = public)
        function obj = controller6(quadrotor)

            obj.zdes = zeros(12,1);
            obj.u0 = repmat(quadrotor.m*quadrotor.g/4, [4,1]);
            [obj.A, obj.B] = linearize(zeros(12,1), obj.u0, quadrotor);
            % [obj.A, obj.B] = linearize_quad(quadrotor, [0, 0, 1]);

            Q = diag([1 1 1 1 1 1 1 1 1 1 1 1]);
            R = 2*eye(4);
     
            obj.K_to_capture = lqr(obj.A,obj.B,Q,R);


            Q = diag([5 5 5 10 10 10 2.5 2.5 2.5 7.5 7.5 7.5]);
            obj.K_to_return = lqr(obj.A,obj.B,Q,R);

            obj.y_prev = [NaN; NaN; NaN];
            obj.y_prev_init = [NaN; NaN; NaN];

            obj.epsilon = quadrotor.l/2;

            obj.uvec = [];
            obj.prev_ref = 0;
            obj.ref_vel = zeros(3,1);
            obj.uav_vel = zeros(3,1);
            obj.compVel = zeros(3,1);

            obj.last_e = zeros(12,1);
            obj.time_step = 0.01;
            obj.prev_coeffs = zeros(3, 6);
        end

        % function [ref_pos, ref_vel] = calc_ref(obj, z, y)
        %     if obj.y_prev ~= obj.y_prev_init % if after t0; have history of UAV
        %         delta_y = y - obj.y_prev; % change in UAV position
        %         uav_vel = delta_y / obj.time_step; % UAV velocity over past time-step
        %         normalized_direction = (delta_y)/norm(delta_y); % heading of UAV
        %         dist_quad_uav = norm(z(1:3) - y(1:3)); % scalar distance between Quad and UAV
        %         ref_pos = y(1:3) + dist_quad_uav * normalized_direction; % cast a reference position 
        %             % by projecting UAV heading by the scalar distance between Q & U
        %     else % if no history available, cast reference position as current UAV position
        %         ref_pos = y(1:3);
        %     end
        %     ref_pos = y(1:3);
        %     obj.y_prev = y;
        %     % obj.projdist = cat(2,obj.projdist,ref);
        %     if norm(y(1:3) - z(1:3)) == 0
        %         comp_vel = zeros(3,1)
        %     else
        %         if norm(uav_vel) == 0
        %             comp_vel = (y(1:3) - z(1:3))/norm(y(1:3) - z(1:3))
        %         else
        %             comp_vel = norm(uav_vel)*((y(1:3) - z(1:3))/norm(y(1:3) - z(1:3)))
        %         end
        %     end
        %     ref_vel = (uav_vel + comp_vel);
        %     obj.compVel = cat(2,obj.compVel,comp_vel);
        %     obj.uav_vel = cat(2,obj.uav_vel,uav_vel);
        %     obj.ref_vel = cat(2,obj.ref_vel,ref_vel);
        % end

        function points = castRefPoints(self, y)
            ydot = (y - self.prev_coeffs(:,6)) / self.timeStep;
            yddot = (ydot - self.prev_coeffs(:, 5)) / self.timeStep;
            y3dot = (yddot - self.prev_coeffs(:, 4)*2) / self.timeStep;
            y4dot = (y3dot - self.prev_coeffs(:, 3)*6) / self.timeStep;
            y5dot = (y4dot - self.prev_coeffs(:, 2)*24) / self.timeStep;
        
            a0 = y;
            a1 = ydot;
            a2 = yddot / 2;
            a3 = y3dot / 6;
            a4 = y4dot / 24;
            a5 = y5dot / 120;
            self.prev_coeffs = [a5, a4, a3, a2, a1, a0];
            coeffs = self.prev_coeffs;
        
            times = [self.time_step, self.time_step*2, self.time_step*3, self.time_step*4, self.time_step*5];
            for i=1:len(times)
                points(i) = coeffs*[times(i)^5; times(i)^4; times(i)^3; times(i)^2; times(i); 1];
            end
        end

        function [ref_pos, ref_vel] = calc_ref(obj, z, y)
            % if self.output_count > 15
            %     %predict intruder traj
            %     coeffs = self.solveCoeffs(y);
            %     r(1:3) = findIntersectionPoint(coeffs, z, self.maxVels);
            % else
            %     r(1:3) = y(1:3);
            %     self.output_count = self.output_count + 1;
            if any(isnan(obj.y_prev)) == 0
                delta_y = y - obj.y_prev; % change in UAV position
                uav_vel = delta_y / obj.time_step; % UAV velocity over past time-step
                uav_speed = norm(uav_vel);
                if norm(delta_y) == 0
                    uav_heading = zeros(3,1);
                else
                    uav_heading = (delta_y)/norm(delta_y); % heading of UAV
                end
            %     dist_quad_uav = norm(z(1:3) - y(1:3)); % scalar distance between Quad and UAV
            %     ref_pos = y(1:3) + dist_quad_uav * normalized_direction; % cast a reference position 
            %         % by projecting UAV heading by the scalar distance between Q & U
            else % if no history available, cast reference position as current UAV position
                delta_y = 0;
                ref_pos = y(1:3) + delta_y;
                uav_heading = zeros(3,1);
                uav_speed = 0;
                uav_vel = zeros(3,1);
            end
            ref_pos = y(1:3);
            obj.y_prev = y;
            % obj.projdist = cat(2,obj.projdist,ref);
            if norm(y(1:3) - z(1:3)) == 0
                comp_heading = zeros(3,1);
            else
                comp_heading = (y(1:3) - z(1:3))/norm(y(1:3) - z(1:3));
                % if norm(uav_vel) == 0
                %     comp_heading = (y(1:3) - z(1:3))/norm(y(1:3) - z(1:3))
                % else
                %     comp_vel = norm(uav_vel)*((y(1:3) - z(1:3))/norm(y(1:3) - z(1:3)))
                % end
            end
            % ref_heading = (uav_heading + comp_heading)/norm(uav_heading + comp_heading);
            ref_heading = (delta_y + (y - z(1:3))) / norm(delta_y + (y - z(1:3)));
            vel_Mag = uav_speed;
            % ref_vel = vel_Mag * ref_heading;
            ref_vel = uav_vel;
            obj.compVel = cat(2,obj.compVel,comp_heading/norm(uav_heading + comp_heading)*vel_Mag);
            obj.uav_vel = cat(2,obj.uav_vel,uav_heading/norm(uav_heading + comp_heading)*vel_Mag);
            obj.ref_vel = cat(2,obj.ref_vel,ref_vel);
        end

        function iscaptured = find_capture(obj, z, y)
            if norm(z(1:3)-y,2) < obj.epsilon
                iscaptured = true;
            else
                iscaptured = false;
            end
        end

        function u = output(obj, ~, z, y)
            iscaptured = find_capture(obj, z, y);
            if iscaptured == false
                [ref_pos, ref_vel] = obj.calc_ref(z, y);
                obj.zdes(1:3) = ref_pos;
                obj.zdes(7:9) = ref_vel;
                u = obj.u0 + obj.K_to_capture * (obj.zdes - z)
            elseif iscaptured == true
                obj.zdes = zeros(12,1);
                obj.zdes(3) = 3;
                u = obj.u0 + obj.K_to_return * (obj.zdes - z)
            end
            obj.uvec = cat(2,obj.uvec,u);
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