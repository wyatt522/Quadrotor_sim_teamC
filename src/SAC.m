classdef SAC < handle
    properties(Access = public)
        altitude(1,1) double;
        k(4, 12) double;
        u0(1,1) double;
        A(12,12) double;
        B(12,4) double;
        prev_coeffs; 
        timeStep(1,1) double;
        output_count(1,1) double;
        target_time(1,1) double;
        jump_ahead_level(1,1) double;
        fast_k(4, 12) double;
        home_k(4, 12) double;

    end

    methods(Access = public)
        function obj = SAC(quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.timeStep = 0.01;
            position = [0,0,1];
            Q = diag([4, 4,6,15, 15, 4,1.5,1.5,1.5,3,3,1.5]);
            R = 2*eye(4);

            fast_Q = diag([14, 14, 16, 12, 12, 4,1.0,1.0,1.0,1.5,1.5,1.5]);
            fast_R = 0.3*eye(4);

            home_Q = diag([4, 4, 6, 18, 18, 10,8.0,8.0,8.0,4.5,4.5,4.5]);
            home_R = 3*eye(4);

            [A, B] = linearize_quad(quadrotor, position);
            [K,~, ~] = lqr(A,B,Q,R);
            [fast_K,~, ~] = lqr(A,B,fast_Q,fast_R);
            [home_K,~, ~] = lqr(A,B,home_Q,home_R);

            obj.k = K;
            obj.fast_k = fast_K;
            obj.home_k = home_K;
            obj.output_count = 0;
            obj.prev_coeffs = zeros(3, 8);
            obj.target_time = -1;
            obj.jump_ahead_level = 3.6;

        end

        function u = output(self, isCaptured, z, y)
            if isCaptured == false
                %find desired point
                r = zeros(12, 1);
                % determine trajectory of uav
                [~, coeff_size] = size(self.prev_coeffs);
                coeffs = self.solveCoeffs(y, coeff_size - 1);
                % collect data first 15 timesteps
                
                error_vec = z(1:3) - y;
                error_mag = norm(error_vec);
                if (self.output_count > 10)
                    % no target time, attempt to aquire one
                    if (self.target_time == -1)
                        if (error_mag < 1)
                            r(1:3) = solvePoly(coeffs, 0.03);
                            disp("in kill mode");
                            disp(error_mag);
                        else
                            self.target_time = self.jump_ahead_level;
                            self.jump_ahead_level = self.jump_ahead_level + 0.1;
                            r(1:3) = solvePoly(coeffs, self.target_time);
                        end
                    else
                        % update for time passing, and go to the expected
                        % location
                        self.target_time = max(self.target_time - self.timeStep, 0);
                        if ((self.target_time == 0) || (error_mag < 1))
                            self.target_time = -1;
                        end
                        r(1:3) = solvePoly(coeffs, self.target_time);

                    end
                    
                else
                    % update data collection timer, go straight to UAV
                    self.output_count = self.output_count + 1;
                    r(1:3) = y(1:3);
                end
                if error_mag < 1
                    temp_k = self.fast_k; % When within range of UAV
                else
                    temp_k = self.k;
                end

                %bound target location
                r(1:2) = min(5, max(-5, r(1:2)));
                r(3) = min(10, max(0, r(3)));

                u = repmat(self.u0, [4,1]) + temp_k*(r - z);
            else
                home = zeros(12, 1);
                if not((abs(z(1)) < 0.5) && (abs(z(2)) < 0.5))
                    home(3) = z(3);
                else
                    home(3) = 0.2;
                end
                u = repmat(self.u0, [4,1]) + self.home_k*(home - z);
            end
        end

        function coeffs = solveCoeffs(self, y, degree)
            threshold = 1e-11; % threshold for rounding errors
            yder = zeros(3, degree + 1);
            yder(:, 1) = y;
            for i = 2:degree + 1
                yder(:, i) = (yder(:, i - 1) - factorial(i - 2)*self.prev_coeffs(:, (degree + 3 - i))) / self.timeStep;
                for j = 1:3
                    if (abs(yder(j, i)) < threshold)
                        yder(j, i) = 0;
                    end
                end
            end
            coeffs = zeros(3, degree + 1);
            for i = 1:degree + 1
                coeffs(:, i) = yder(:, degree + 2 - i) / factorial(degree + 1 - i);
            end

            self.prev_coeffs = coeffs;
        end
    end
end