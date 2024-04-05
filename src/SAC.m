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
        function obj = SAC(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.prev_y = zeros(3,2);
            obj.timeStep = 0.01;
            obj.maxVels = [0.75; 0.75; 2.25];
            position = [1,0,0];
            Q = diag([1,1,1,0.5,0.5,0.5,0,0,0,2,2,2]);
            R = eye(4);
            [A, B] = linearize_quad(quadrotor, position);
            [K,~, ~] = lqr(A,B,Q,R);
            obj.k = K;
        end

        function [u, r] = output(self, isCaptured, z, y)
            if isCaptured == false
                %predict intruder traj
                coeffs = self.solveCoeffs(y);
                % disp([target, y])
                %find desired point
                r = zeros(12, 1);
                if ~(all(self.prev_y == 0))
                    r(1:3) = findIntersectionPoint(coeffs, z, self.maxVels);
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
    end
end