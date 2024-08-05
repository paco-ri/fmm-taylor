function m0 = debyem0(sigma,lambda)
%DEBYEM0 compute part of the Debye current associated with sigma
%   Detailed explanation goes here

pdo = [];
pdo.lap = 1;

% resample
n = size(sigma.vals{1,1},1);
sigma = resample(sigma,n+1);
dom = sigma.domain;

% solve lap(u) = sigma
L = surfaceop(dom, pdo, sigma);
L.rankdef = true;
u = L.solve();

vn = normal(dom);
m0 = 1i.*lambda.*(grad(u) + 1i.*cross(vn, grad(u)));

% m0err = div(m0) - 1i*lambda.*sigma;
% fprintf('upsampled m0err = %f\n', norm(m0err))

% resample
m0 = resample(m0,n);

end