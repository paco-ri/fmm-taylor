function Bsigma = mtxBsigma(S,dom,sigma,eps,Q)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
%   Arguments:

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals';

% need all arguments
dummyvals = repmat(sigmavals,3,1);
[~,gradS0sigma] = virtualcasing.evalgradcurlS0(S,dummyvals,sigmavals,eps);
gradS0sigma = array_to_surfacefun(gradS0sigma',dom,S); % note transpose 

% derivatives of layer potentials
n = normal(dom);
ngradS0sigma = dot(n,gradS0sigma);

% construct Bsigma
Bsigma = -sigma/2 - ngradS0sigma;

end