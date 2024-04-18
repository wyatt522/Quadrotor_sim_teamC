classdef SAC < handle
    properties(Access = public)
        altitude(1,1) double;
        k(4, 12) double;
        u0(1,1) double;
        A(12,12) double;
        B(12,4) double;
        prev_coeffs(3, 6) double;
        timeStep(1,1) double;
        maxVels(3,1) double;
        output_count(1,1) double;

    end

    methods(Access = public)
        function obj = SAC(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.timeStep = 0.01;
            obj.maxVels = [0.25; 0.25; 1];
            position = [1,0,0];
            Q = diag([10,10,10,2.5,2.5,2.5,5,5,5,5,5,5]);
            R = 3*eye(4);
            [A, B] = linearize_quad(quadrotor, position);
            [K,~, ~] = lqr(A,B,Q,R);
            obj.k = K;
            obj.output_count = 0;
            obj.prev_coeffs = zeros(3, 6);
        end

        function [u, r] = output(self, isCaptured, z, y)
            if isCaptured == false
                % disp([target, y])
                %find desired point
                r = zeros(12, 1);
                if self.output_count > 15
                    %predict intruder traj
                    coeffs = self.solveCoeffs(y);
                    r(1:3) = findIntersectionPoint(coeffs, z, self.maxVels);
                else
                    r(1:3) = y(1:3);
                    self.output_count = self.output_count + 1;
                end
                u = repmat(self.u0, [4,1]) + self.k*(r - z);
            else
                home = zeros(12, 1);
                home(3) = 1;
                u = repmat(self.u0, [4,1]) + self.k*(home - z);
            end
        end

        function coeffs = solveCoeffs(self, y)
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
        end
    end
end