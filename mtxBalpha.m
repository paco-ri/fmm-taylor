function Balpha = mtxBalpha(S,dom,zpars,mH,eps,Q)
%MTXBALPHA compute alpha coefficient for surface magnetic field
%   Detailed explanation goes here

mHvals = surfacefun_to_array(mH,dom,S);

% single layer potential evals
SkmH1 = helm3d.dirichlet.eval(S,zpars,mHvals(:,1),eps,S,Q);
SkmH2 = helm3d.dirichlet.eval(S,zpars,mHvals(:,2),eps,S,Q);
SkmH3 = helm3d.dirichlet.eval(S,zpars,mHvals(:,3),eps,S,Q);
SkmH = [SkmH1 SkmH2 SkmH3];
SkmH = array_to_surfacefun(SkmH,dom,S);

n = normal(dom);
nSkmH = dot(n,SkmH);

% derivatives of layer potentials 

curlSkmH = curl(SkmH);
ncurlSkmH = dot(n,curlSkmH);

lambda = zpars(1);

Balpha = 1i*(lambda*nSkmH + ncurlSkmH);

end