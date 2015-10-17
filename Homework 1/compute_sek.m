function SEKernel = compute_sek(xi, xj, kernelScaleFactor, kernelLengthScale)
   % SEKernel = kernelScaleFactor.^2 * exp(-0.5 * (1/(kernelLengthScale.^2)) * (xi-xj).^2);
    SEKernel = kernelScaleFactor.^5 * exp(-0.5 * dot(transpose(xi - xj), (1/(kernelLengthScale.^2)) * (xi - xj)));
end