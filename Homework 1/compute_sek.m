function SEKernel = compute_sek(xi, xj, kernelScaleFactor, kernelLengthScale)
     SEKernel = kernelScaleFactor.^2 * exp(-0.5 * dot((xi - xj) .* (kernelLengthScale.^-2), transpose(xi - xj)));
end