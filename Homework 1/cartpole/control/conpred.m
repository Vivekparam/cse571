%% conpred.m 
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-01-24
%
% Modified by Arunkumar Byravan for CSE 571: Probabilistic Robotics 
% (Autumn 2015) at  University of Washington

function [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdp, dSdp, dCdp] ...
  = conpred(policy, m, s)
%% Code

% 1. Extract policy parameters
policy.hyp = policy.p.hyp;
policy.inputs = policy.p.inputs;
policy.targets = policy.p.targets;

% fix policy signal and the noise variance 
% (avoids some potential numerical problems)
policy.hyp(end-1,:) = log(1);                  % set signal variance to 1
policy.hyp(end,:) = log(0.01);                 % set noise standard dev to 0.01

if nargout < 4                                 % if no derivatives are required
  [M, S, C] = pred(policy, m, s);
else                                             % else compute derivatives too
  [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdi, dSdi, dCdi, dMdt, ...
    dSdt, dCdt, dMdh, dSdh, dCdh] = predd(policy, m, s);
  
  % 3. Set derivatives of non-free parameters to zero: signal and noise variance
  d = size(policy.inputs,2);            
  d2 = size(policy.hyp,1); dimU = size(policy.targets,2);
  sidx = bsxfun(@plus,(d+1:d2)',(0:dimU-1)*d2);
  dMdh(:,sidx(:)) = 0; dSdh(:,sidx(:)) = 0; dCdh(:,sidx(:)) = 0;
  
  % 4. Merge derivatives
  dMdp = [dMdh dMdi dMdt]; dSdp = [dSdh dSdi dSdt]; dCdp = [dCdh dCdi dCdt];
end

function [M, S, V] = pred(model, m, s)
%% Code
persistent iK oldX oldIn oldOut beta oldn;
D = size(model.inputs,2);    
[n, E] = size(model.targets);      

input = model.inputs;  target = model.targets; X = model.hyp;

if numel(X) ~= numel(oldX) || isempty(iK) ||  n ~= oldn || ...
    sum(any(X ~= oldX)) || sum(any(oldIn ~= input)) || ...
    sum(any(oldOut ~= target))
  oldX = X; oldIn = input; oldOut = target; oldn = n;
  K = zeros(n,n,E); iK = K; beta = zeros(n,E);
  
  for i=1:E                                              
    inp = bsxfun(@rdivide,model.inputs,exp(X(1:D,i)'));
    dst = bsxfun(@plus,sum(inp.*inp,2),sum(inp.*inp,2)')-2*(inp*inp');
    K(:,:,i) = exp(2*X(D+1,i)-dst/2);
    if isfield(model,'nigp')
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n) + diag(model.nigp(:,i)))';
    else
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n))';
    end
    iK(:,:,i) = L'\(L\eye(n));
    beta(:,i) = L'\(L\model.targets(:,i));
  end
end

k = zeros(n,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E);

inp = bsxfun(@minus,model.inputs,m');                    

for i=1:E
  iL = diag(exp(-X(1:D,i))); 
  in = inp*iL;
  B = iL*s*iL+eye(D);
  
  t = in/B;
  l = exp(-sum(in.*t,2)/2); lb = l.*beta(:,i);
  tL = t*iL;
  c = exp(2*X(D+1,i))/sqrt(det(B));
  
  M(i) = sum(lb)*c;                                            
  V(:,i) = tL'*lb*c;                   
  k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2;
end

if nargout > 1

  for i=1:E
    ii = bsxfun(@rdivide,inp,exp(2*X(1:D,i)'));

    for j=1:i
      R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D);
      t = 1/sqrt(det(R));
      ij = bsxfun(@rdivide,inp,exp(2*X(1:D,j)'));
      Q = R\s/2;
      aQ = ii*Q; Kn = bsxfun(@plus,sum(aQ.*ii,2),sum((-ij)*Q.*(-ij),2)')+2*aQ*ij';
      L = exp(bsxfun(@plus,k(:,i),k(:,j)')+Kn);
      S(i,j) = t*beta(:,i)'*L*beta(:,j); S(j,i) = S(i,j);
    end

    S(i,i) = S(i,i) + 1e-6;          

  end

  S = S - M*M';

end
