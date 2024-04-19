function futurePose = solvePoly(coeffs, t)
    [~, col] = size(coeffs);
    time = [];
    for i = (col-1):-1:0
        time = [time; (t^(i))];
    end
    futurePose = coeffs*time;
end