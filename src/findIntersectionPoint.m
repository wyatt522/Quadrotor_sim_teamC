function t = findIntersectionPoint(coeffs, z, maxVels)
    for t = 0:0.2:8
        uavPose = solvePoly(coeffs, t);
        min = z(1:3) - maxVels*t;
        max = z(1:3) + maxVels*t;
        if all((uavPose < max) & (uavPose > min))
            desired = uavPose;
            disp("found intersection")
            disp(desired)
            disp(t)
            return
        end
    end
    disp("did not find intersection")
    t = -1;
end