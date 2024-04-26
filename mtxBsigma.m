function Bsigma = mtxBsigma(S,dom,zpars,sigma,eps,Q)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
%   Arguments:

lambda = zpars(1);
m0 = debyem0(sigma,lambda);

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals';
m0vals = surfacefun_to_array(m0,dom,S);
m0vals = m0vals';

% single layer potential evals
S0m01 = helm3d.dirichlet.eval(S,zpars,m0vals(1,:),eps,S,Q);
S0m02 = helm3d.dirichlet.eval(S,zpars,m0vals(2,:),eps,S,Q);
S0m03 = helm3d.dirichlet.eval(S,zpars,m0vals(3,:),eps,S,Q);
S0m0 = [S0m01 S0m02 S0m03];
S0m0 = array_to_surfacefun(S0m0,dom,S);

% need all arguments
[curlS0m0,gradS0sigma] = virtualcasing.evalgradcurlS0(S,m0vals,sigmavals,eps,S,Q);
curlS0m0 = array_to_surfacefun(curlS0m0',dom,S); % note transpose 
gradS0sigma = array_to_surfacefun(gradS0sigma',dom,S); % note transpose 

% derivatives of layer potentials
n = normal(dom);
nSkm0 = dot(n,S0m0);
ncurlSkm0 = dot(n,curlS0m0); 
ngradSksigma = dot(n,gradS0sigma);

% construct Bsigma
Bsigma = -sigma/2 - ngradSksigma + 1i*(lambda*nSkm0 + ncurlSkm0);

end