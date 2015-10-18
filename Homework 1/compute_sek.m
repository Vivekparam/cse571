function SEKernel = compute_sek(xi, xj, kernelScaleFactor, kernelLengthScale)
   % SEKernel = kernelScaleFactor.^2 * exp(-0.5 * (1/(kernelLengthScale.^2)) * (xi-xj).^2);
     SEKernel = kernelScaleFactor.^2 * exp(-0.5 * dot((xi - xj) .* (kernelLengthScale.^-2), transpose(xi - xj)));
%   m = kernelLengthScale .^ (-2);
%   SEKernel = (kernelScaleFactor .^ 2) *  exp(-(dot((xi - xj) .* m, transpose(xi - xj))) / 2);
end