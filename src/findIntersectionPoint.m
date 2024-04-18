

function desired = findIntersectionPoint(coeffs, z, maxVels)
    for t = 0:0.2:8
        uavPose = solvePredictQuintic(coeffs, t);
        min = z(1:3) - maxVels*t;
        max = z(1:3) + maxVels*t;
        if all((uavPose < max) & (uavPose > min))
            desired = uavPose;
            disp("found intersection")
            disp(desired)
            return
        end
    end
    desired = solvePredictQuintic(coeffs, 0);
    disp("did not find intersection")
end

function futurePose = solvePredictQuadratic(coeffs, t)
    futurePose = coeffs*[t^2; t; 1];
end

function futurePose = solvePredictQuintic(coeffs, t)
    futurePose = coeffs*[t^5; t^4; t^3; t^2; t; 1];
end

function futurePose = solvePredict10thPower(coeffs, t)
    futurePose = coeffs*[t^10; t^9; t^8; t^7; t^6; t^5; t^4; t^3; t^2; t; 1];
end