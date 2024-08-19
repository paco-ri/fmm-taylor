function m0 = debyem0(sigma,lambda,varargin)
%DEBYEM0 compute part of the Debye current associated with sigma
%   
%    Arguments
%      sigma [surfacefun] 
%      lambda [double] 
%      np [int>0] polynomial degree by which to oversample

oversample = false;
if nargin > 2
    oversample = true;
    np = varargin{1};
end

pdo = [];
pdo.lap = 1;

% resample
n = size(sigma.vals{1,1},1);
if oversample
    sigma = resample(sigma,n+np);
end

dom = sigma.domain;

% solve lap(u) = sigma
L = surfaceop(dom, pdo, sigma);
L.rankdef = true;
u = L.solve();

vn = normal(dom);
m0 = 1i.*lambda.*(grad(u) + 1i.*cross(vn, grad(u)));

% m0err = div(m0) - 1i*lambda.*sigma;
% fprintf('upsampled m0err = %f\n', norm(m0err))

if oversample
    % resample
    m0 = resample(m0,n);
end

end