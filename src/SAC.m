classdef SAC < handle
    properties(Access = public)
        altitude(1,1) double;
        k(4, 12) double;
        u0(1,1) double;
        A(12,12) double;
        B(12,4) double;
        prev_y(3,2) double;
        timeStep(1,1) double;
        maxVels(3,1) double;
    end

    methods(Access = public)
        function obj = SAC(gains, quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.k = gains;
            obj.prev_y = zeros(3,2);
            obj.timeStep = 0.01;
            obj.maxVels = [0.75; 0.75; 2.25];
            position = [1,0,0];
            [obj.A, obj.B] = linearize_quad(position);
        end

        function [u, r] = output(self, isCaptured, z, y)
            if isCaptured == false
                %predict intruder traj
                coeffs = self.solveCoeffs(y);
                % disp([target, y])
                %find desired point
                r = zeros(12, 1);
                if ~(all(self.prev_y == 0))
                    r(1:3) = self.findIntersectionPoint(coeffs, z);
                else
                    r(1:3) = y(1:3);
                end
                u = repmat(self.u0, [4,1]) + self.k*(r - z);
                %solve for u
                %return u
            end
        end

        function coeffs = solveCoeffs(self, y)
            ydot = (y - self.prev_y(:,1)) / self.timeStep;
            yddot = (ydot - ((self.prev_y(:,1) - self.prev_y(:,2)) / self.timeStep)) / self.timeStep;
            c = y;
            b = ydot;
            a = yddot / 2;
            self.prev_y(:,2) = self.prev_y(:,1);
            self.prev_y(:,1) = y;
            coeffs = [a, b ,c];
        end

        function futurePose = solvePredict(self, coeffs, t)
            futurePose = coeffs*[t^2; t; 1];
        end

        function desired = findIntersectionPoint(self, coeffs, z)
            for t = 0:0.2:8
                uavPose = self.solvePredict(coeffs, t);
                min = z(1:3) - self.maxVels*t;
                max = z(1:3) + self.maxVels*t;
                if all((uavPose < max) & (uavPose > min))
                    desired = uavPose;
                    disp("found intersection")
                    return
                end
            end
            desired = self.solvePredict(coeffs, 0);
            disp("did not find intersection")
        end
    end
end

function [A, B] = linearize_quad(position)
            syms x y z xdot ydot zdot phi theta psiVar omega1 omega2 omega3 u1 u2 u3 u4 ...
                n1 n2 n3 r1 r2 r3
            l = 0.2; m = 0.5; I1 = 1.24; I2 = 1.24; I3 = 2.48; g = 9.8; sigma = 0.01; % all in quad object
            
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