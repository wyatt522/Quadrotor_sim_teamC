classdef lineartrial2 < handle
    properties(Access = public)
        zdes(12,1) double;
        K(4,12) double;
        u0(4,1) double;
        prev_y(3,2) double;
        timeStep(1,1) double;
        maxVels(3,1) double;
    end


    methods(Access = public)
        function obj = lineartrial2(gains, quadrotor)
            obj.zdes = zeros(12,1);
            obj.u0 = repmat(quadrotor.m*quadrotor.g/4, [4,1]);
            obj.K = gains;
        end

        % function [u, ref] = output(obj, isCaptured, z, y)
        %     if isCaptured == false
        %         %predict intruder traj
        %         coeffs = obj.solveCoeffs(y);
        %         target = obj.solvePredict(coeffs, 0.01);
        %         % disp([target, y])
        %         %find desired point
        %         if ~(all(obj.prev_y == 0))
        %             ref = obj.findIntersectionPoint(coeffs, z);
        %             % disp([ref, y])
        %         else
        %             ref = y(1:3);
        %         end
        % 
        %         %solve for u
        %         %return u
        %     end
        %     obj.zdes(1:3) = ref(1:3);
        %     % disp(['ref = ', ref'])
        %     disp(ref)
        %     e = obj.zdes - z;
        %     v = -obj.K * e;
        %     u = obj.u0 - v;
        %     % disp(e);
        % end
        % 
        % function coeffs = solveCoeffs(obj, y)
        %     ydot = (y - obj.prev_y(:,1)) / obj.timeStep;
        %     yddot = (ydot - ((obj.prev_y(:,1) - obj.prev_y(:,2)) / obj.timeStep)) / obj.timeStep;
        %     c = y;
        %     b = ydot;
        %     a = yddot / 2;
        %     obj.prev_y(:,2) = obj.prev_y(:,1);
        %     obj.prev_y(:,1) = y;
        %     coeffs = [a, b ,c];
        % end
        % 
        % function futurePose = solvePredict(obj, coeffs, t)
        %     futurePose = coeffs*[t^2; t; 1];
        % end
        % 
        % function desired = findIntersectionPoint(obj, coeffs, z)
        %     for t = 0:0.1:2
        %         uavPose = obj.solvePredict(coeffs, t);
        %         min = z(1:3) - obj.maxVels*t;
        %         max = z(1:3) + obj.maxVels*t;
        %         if all((uavPose < max) & (uavPose > min))
        %             desired = uavPose;
        %             return
        %         end
        %         desired = z(1:3);
        %     end
        % end

        function u = output(obj, ~, z, y)
            obj.zdes(1:3) = y(1:3);
            e = obj.zdes - z;
            v = -obj.K * e;
            u = obj.u0 - v;
            disp(u);
        end
    end


end