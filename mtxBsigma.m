function Bsigma = mtxBsigma(S,dom,zpars,sigma,eps,Q)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
%   Arguments:

lambda = zpars(1);
m0 = debyem0(sigma,lambda);

sigmavals = surfacefun_to_array(sigma,dom,S);
m0vals = surfacefun_to_array(m0,dom,S);

% single layer potential evals
Skm01 = helm3d.dirichlet.eval(S,zpars,m0vals(:,1),eps,S,Q);
Skm02 = helm3d.dirichlet.eval(S,zpars,m0vals(:,2),eps,S,Q);
Skm03 = helm3d.dirichlet.eval(S,zpars,m0vals(:,3),eps,S,Q);
Skm0 = [Skm01 Skm02 Skm03];
Skm0 = array_to_surfacefun(Skm0,dom,S);

Sksigma = helm3d.dirichlet.eval(S,zpars,sigmavals,eps,S,Q);
Sksigma = array_to_surfacefun(Sksigma,dom,S);

% derivatives of layer potentials
n = normal(dom);
% curlSkm0 = curl(Skm0); WRONG! this is not surface curl
nSkm0 = dot(n,Skm0);
ncurlSkm0 = dot(n,curlSkm0); 

gradSksigma = grad(Sksigma);
ngradSksigma = dot(n,gradSksigma);

% construct Bsigma
Bsigma = -sigma/2 - ngradSksigma + 1i*(lambda*nSkm0 + ncurlSkm0);

end