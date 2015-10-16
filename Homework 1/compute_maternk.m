function MaternKernel = compute_maternk(xi, xj, v, l)
    % v = kernelScaleFactor, l = kernelLengthScale
    r = abs(xi - xj);
    if r == 0
        MaternKernel = 1;
    else
        MaternKernel = ((2.^(1-v))/gamma(v)) * ((sqrt(2*v) * r) / l).^v * besselk(v, ((sqrt(2*v) * r) / l));
    end
end