%-------------------------------------------------------
% noisy version of prediction
% predicts the new state given the current state and motion
% motion in form of [drot1,dtrans,drot2]
%-------------------------------------------------------
function state = sampleOdometry(motion, state, alphas)

variances = noiseFromMotion(motion, alphas);

noisymotion(1)=sample(motion(1),variances(1),1);
noisymotion(2)=sample(motion(2),variances(2),1);
noisymotion(3)=sample(motion(3),variances(3),1);

state(3)=state(3)+noisymotion(1);
state(1)=state(1)+noisymotion(2)*cos(state(3));
state(2)=state(2)+noisymotion(2)*sin(state(3));
state(3)=state(3)+noisymotion(3);
state(3)=minimizedAngle(state(3));
