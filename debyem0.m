function m0 = debyem0(sigma,lambda)
%DEBYEM0 compute part of the Debye current associated with sigma
%   Detailed explanation goes here

pdo = [];
pdo.lap = 1;

% resample
% n = size(sigma.vals{1,1},1);
% sigma = resample(sigma,2*n);

dom = sigma.domain;
% S = surfer.surfacemesh_to_surfer(dom);

% solve lap(u) = sigma
L = surfaceop(dom, pdo, sigma);
L.rankdef = true;
u = L.solve();

Lr = surfaceop(dom, pdo, real(sigma));
Lr.rankdef = true;
ur = Lr.solve();
Li = surfaceop(dom, pdo, imag(sigma));
Li.rankdef = true;
ui = Li.solve();
u = ur + 1i.*ui;

vn = normal(dom);
m0 = 1i.*lambda.*(grad(u) + 1i.*cross(vn, grad(u)));
% uvals = surfacefun_to_array(u,dom,S);
% graduvals = get_surface_grad(S,uvals);
% gradu = array_to_surfacefun(graduvals.',dom,S);
% m0 = 1i.*lambda.*(gradu + 1i.*cross(vn, gradu));

% m0err = div(m0) - 1i*lambda.*sigma;
% fprintf('upsampled m0err = %f\n', norm(m0err))

% resample
% m0 = resample(m0,n);

end