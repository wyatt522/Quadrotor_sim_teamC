function futurePose = solveQuintic(coeffs, t)
    futurePose = coeffs*[t^5; t^4; t^3; t^2; t; 1];
end