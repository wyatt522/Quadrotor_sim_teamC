classdef SAC < handle
    properties(Access = public)
        altitude(1,1) double;
        k(4, 12) double;
        u0(1,1) double;
        A(12,12) double;
        B(12,4) double;
        prev_coeffs(3, 8) double;
        timeStep(1,1) double;
        maxVels(3,1) double;
        output_count(1,1) double;
        target_time(1,1) double;
        jump_ahead_level(1,1) double;
        fast_k(4, 12) double;

    end

    methods(Access = public)
        function obj = SAC(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.timeStep = 0.01;
            obj.maxVels = [0.15; 0.15; 1];
            position = [1,0,0];
            Q = diag([4, 4,4,15, 15, 4,1.5,1.5,1.5,1.5,1.5,1.5]);
            fast_Q = diag([10,10,10, 8,8, 4,1.5,1.5,1.5,1.5,1.5,1.5]);
            R = 3*eye(4);
            [A, B] = linearize_quad(quadrotor, position);
            [K,~, ~] = lqr(A,B,Q,R);
            [fast_K,~, ~] = lqr(A,B,fast_Q,R);
            obj.k = K;
            obj.output_count = 0;
            obj.prev_coeffs = zeros(3, 8);
            obj.target_time = -1;
            obj.jump_ahead_level = 4.0;
            obj.fast_k = fast_K;

        end

        function [u, r] = output(self, isCaptured, z, y)
            if isCaptured == false
                % disp([target, y])
                %find desired point
                r = zeros(12, 1);
                % determine trajectory of uav
                [~, coeff_size] = size(self.prev_coeffs);
                coeffs = self.solveCoeffs(y, coeff_size - 1);
                % collect data first 100 timesteps
                
                error_vec = z(1:3) - y;
                error_mag = norm(error_vec);

                if (self.output_count > 100 && error_mag > 1)
                    % no target time, attempt to aquire one
                    if (self.target_time == -1)
                        % check to find target, go straight to UAV
                        self.target_time = self.jump_ahead_level;
                        self.jump_ahead_level = self.jump_ahead_level + 0.1;
                    else
                        % update for time passing, and go to the expected
                        % location
                        self.target_time = max(self.target_time - self.timeStep, 0);
                        if (self.target_time == 0)
                            self.target_time = -1;
                        end
                    end
                    r(1:3) = solvePoly(coeffs, self.target_time);
                else
                    % update data collection timer, go straight to UAV
                    self.output_count = self.output_count + 1;
                    r(1:3) = y(1:3);
                end
                if error_mag < 1
                    temp_k = self.fast_k;
                else
                    temp_k = self.k;
                end
                u = repmat(self.u0, [4,1]) + temp_k*(r - z);
            else
                home = zeros(12, 1);
                home(3) = 1;
                u = repmat(self.u0, [4,1]) + self.k*(home - z);
            end
        end

        function coeffs = solveCoeffs(self, y, degree)
            yder = zeros(3, degree + 1);
            yder(:, 1) = y;
            for i = 2:degree + 1
                yder(:, i) = (yder(:, i - 1) - factorial(i - 2)*self.prev_coeffs(:, (degree + 3 - i))) / self.timeStep;
            end
            coeffs = zeros(3, degree + 1);
            for i = 1:degree + 1
                coeffs(:, i) = yder(:, degree + 2 - i) / factorial(degree + 1 - i);
            end

            self.prev_coeffs = coeffs;
        end
    end
end