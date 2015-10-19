function variances = noiseFromMotion(motion, alphas)
variances(1)=alphas(1)*motion(1)^2+alphas(2)*motion(2)^2;
variances(2)=alphas(3)*motion(2)^2+alphas(4)*(motion(1)^2+motion(3)^2);
variances(3)=alphas(1)*motion(3)^2+alphas(2)*motion(2)^2;