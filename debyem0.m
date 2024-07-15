function m0 = debyem0(sigma,lambda)
%DEBYEM0 compute part of the Debye current associated with sigma
%   Detailed explanation goes here

pdo = [];
pdo.lap = 1;
dom = sigma.domain;
n = normal(dom);

% solve lap(u) = sigma
L = surfaceop(dom, pdo, sigma);
L.rankdef = true;
u = solve(L);

m0 = 1i.*lambda.*(grad(u) + 1i.*cross(n, grad(u)));
% m0 = surfacefunv(dom);
% gradu = grad(u);
% nxgradu = cross(n, gradu);
% m0.components{1} = lambda.*(1i.*gradu.components{1}-nxgradu.components{1});
% m0.components{2} = lambda.*(1i.*gradu.components{2}-nxgradu.components{2});
% m0.components{3} = lambda.*(1i.*gradu.components{3}-nxgradu.components{3});

end