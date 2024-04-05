

function desired = findIntersectionPoint(coeffs, z, maxVels)
    for t = 0:0.2:8
        uavPose = solvePredict(coeffs, t);
        min = z(1:3) - maxVels*t;
        max = z(1:3) + maxVels*t;
        if all((uavPose < max) & (uavPose > min))
            desired = uavPose;
            disp("found intersection")
            return
        end
    end
    desired = solvePredict(coeffs, 0);
    disp("did not find intersection")
end

function futurePose = solvePredict(coeffs, t)
    futurePose = coeffs*[t^2; t; 1];
end