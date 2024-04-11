classdef controller2 < handle
    properties(Access = public)
        xyz (1,3) double;
        altitudeDot (1,1) double;
        k(1,2) double;
        u0(1,1) double;
        prev_error double;
        int_error double;
        uvec;
        error;
        d;
    end


    methods(Access = public)
        function obj = controller2(xyz, altitudeDot, gains, quadrotor)
            obj.u0 = quadrotor.m*quadrotor.g/4;
            obj.altitudeDot = altitudeDot;
            obj.xyz = xyz;
            obj.k = gains;
            obj.prev_error = 0;
            obj.error = 0;
            obj.d = 0;
        end

        function u = output(obj, ~, z, ~)
          
            u = obj.u0 + repmat(obj.k*[(obj.xyz(3) - z(3)); obj.altitudeDot-z(9)],[4,1]);
            ex = obj.xyz(1) - z(1); % pos = needs to move right, neg = needs to move left
            roll = z(4); %rotors 2 and 4
            %pitch = z(5); %rotors 1 and 3
            %b = 100;
            %pu13 = obj.pit(u(1), u(3), ex, 0.15, 0.25);
            %u(1) = pu13(1);
            %u(3) = pu13(2);
            %if((-0.15 <ex)&&(ex>0.15)) %stabilize
             %   u(3) = u(3) + b*abs(z(5));
             %   u(1) = u(1) - b*abs(z(5));
            %end
            %desPitch = 1*ex^2;
            %dPdot = 2*ex;
            %u(1) = u(1) - 2*ex^3 + 1*desPitch;
            %u(3) = u(3) + 2*ex^3 + 1*desPitch;

            i = trapz(0.1, [obj.prev_error, ex]);
            obj.d = gradient(obj.error);

            
            u(1) = u(1) - 0.1*ex - 10*i - 10*obj.d(end);
            u(3) = u(3) + 0.1*ex + 10*i + 10*obj.d(end);

            obj.prev_error = ex;

            disp("ex: "+ ex)
            disp("roll: " + z(4))
            disp("pitch: " + z(5))
            disp("velocity: "+ z(9))
            disp("u: " + u)
            disp("------------------")
            obj.uvec = cat(2, obj.uvec, u);
            obj.error = cat(2, obj.error,  ex);

   
        end
    end
    
    methods(Access = private) % Define pit as a private method
        function pu13 = pit(~, u1, u3, ex, threshold, gain)
            if (ex > threshold) % move right by rotors 3 
                u3 = u3 + gain * abs(ex); 
            end
            if (ex < -threshold) % move left by rotor 1
                u1 = u1 + gain * abs(ex); 
            end 
            pu13 = [u1 u3];
        end
    end


end

